"""Measurements behind docs/variation_design.md. Not part of VirPlot.

Reads a BAM and the FASTA it was mapped to (default: the shipped GLRaV-13
example) and prints, for each reference, two kinds of statistic, kept apart
on purpose:

* within-sample variation: how the reads at one site disagree with each
  other (minor-allele frequency, pi, the allele-frequency spectrum, indels
  carried by a minority, linkage);
* divergence from the reference: how the reads disagree with the sequence
  they were mapped to (sites where the read majority is not the reference,
  mismatches per read, soft-clip pileups, mate pairs, how coverage meets
  the contig ends).

    python3 docs/design/variation_measurements.py
    python3 docs/design/variation_measurements.py --bam S.bam --fasta ref.fa --around RNA2:1546

--around [REF:]POS adds a profile of both kinds of signal in windows either
side of POS. Soft clips appear only with a local aligner (bwa mem, bowtie2
--local); with end-to-end alignment the contig-end and mate checks carry the
weight. Standard library and numpy only. Bases count at base quality >= 30;
sites count with >= 100 such bases.
"""

import argparse
import collections
import gzip
import itertools
import os
import struct

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.join(HERE, "..", "..", "examples")
MIN_BQ, MIN_DEPTH, ISNV, INDEL = 30, 100, 0.02, 0.05
SKIP = 0x4 | 0x100 | 0x200 | 0x400 | 0x800       # unmapped, secondary, QC fail, duplicate, supplementary

SEQ = "=ACMGRSVTWYHKDBN"
CODE = np.full(16, -1, np.int8)
for i, b in enumerate("ACGT"):
    CODE[SEQ.index(b)] = i
OPS = "MIDNSHP=X"


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: np.array(["ACGT".find(c) for c in "".join(v)], np.int8) for k, v in seqs.items()}


def bam_open(path):
    """(references [(name, length)], iterator over raw alignment records)."""
    fp = gzip.open(path, "rb")
    assert fp.read(4) == b"BAM\x01", f"{path} is not a BAM file"
    (lt,) = struct.unpack("<i", fp.read(4))
    fp.read(lt)
    (nref,) = struct.unpack("<i", fp.read(4))
    refs = []
    for _ in range(nref):
        (ln,) = struct.unpack("<i", fp.read(4))
        name = fp.read(ln)[:-1].decode()
        (length,) = struct.unpack("<i", fp.read(4))
        refs.append((name, length))

    def records():
        while True:
            head = fp.read(4)
            if len(head) < 4:
                return
            (size,) = struct.unpack("<i", head)
            yield fp.read(size)
    return refs, records()


def aux_sa(b, off):
    """The SA tag's (reference, 1-based position) partners, if any."""
    out, i = [], off
    while i < len(b):
        tag, t = b[i:i + 2].decode(), chr(b[i + 2])
        i += 3
        if t == "Z" or t == "H":
            j = b.index(b"\0", i)
            if tag == "SA":
                for s in b[i:j].decode().rstrip(";").split(";"):
                    f = s.split(",")
                    out.append((f[0], int(f[1])))
            i = j + 1
        elif t in "cCA":
            i += 1
        elif t in "sS":
            i += 2
        elif t in "iIf":
            i += 4
        elif t == "B":
            sub, (cnt,) = chr(b[i]), struct.unpack_from("<i", b, i + 1)
            i += 5 + cnt * {"c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}[sub]
        else:
            break
    return out


class Pile:
    """Everything counted for one reference sequence."""

    def __init__(self, name, ref):
        L = len(ref)
        self.name, self.ref, self.L = name, ref, L
        self.counts = np.zeros((2, L, 4), np.int64)   # [fragment orientation, site, base] at Q >= MIN_BQ
        self.cov = np.zeros(L + 1, np.int64)           # aligned bases, any quality (difference array)
        self.dels = np.zeros(L + 1, np.int64)          # reads with a deletion over the site (difference array)
        self.ins = np.zeros(L, np.int64)               # reads with an insertion after the site
        self.del_len, self.ins_len = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
        self.starts, self.ends = np.zeros(L + 1, np.int64), np.zeros(L + 1, np.int64)   # aligned first / last base
        self.clips, self.sa = np.zeros(L + 1, np.int64), collections.defaultdict(collections.Counter)
        self.mate_unmapped, self.mate_elsewhere, self.improper = (np.zeros(L, np.int64) for _ in range(3))
        self.paired = 0
        self.reads = []                                 # (start, end, sites, bases, mismatches, indels)

    def add(self, b, refs):
        _rid, pos, l_rn, _mq, _bin, n_cig, flag, l_seq, mate_rid = struct.unpack_from("<iiBBHHHii", b, 0)
        off = 32 + l_rn
        cig = struct.unpack_from(f"<{n_cig}I", b, off)
        off += 4 * n_cig
        packed = np.frombuffer(b, np.uint8, (l_seq + 1) // 2, off)
        off += (l_seq + 1) // 2
        nib = np.empty(2 * len(packed), np.uint8)
        nib[0::2], nib[1::2] = packed >> 4, packed & 15
        base, qual = CODE[nib[:l_seq]], np.frombuffer(b, np.uint8, l_seq, off)
        off += l_seq
        # fragment orientation: a read 2 of a pair speaks for its read 1, flipped
        s = int(bool(flag & 0x10) != bool(flag & 0x1 and flag & 0x80))
        ref = self.ref
        rp, qp, sites, bases, mism, indels = pos, 0, [], [], 0, []
        lead = cig and OPS[cig[0] & 15] in "SH" and cig[0] >> 4 >= 5
        for c in cig:
            n, op = c >> 4, OPS[c & 15]
            if op in "M=X":
                bb, qq = base[qp:qp + n], qual[qp:qp + n]
                idx = np.arange(rp, rp + n)
                mism += int(((bb >= 0) & (bb != ref[idx])).sum())
                ok = (bb >= 0) & (qq >= MIN_BQ)
                np.add.at(self.counts[s], (idx[ok], bb[ok]), 1)
                sites.append(idx[ok])
                bases.append(bb[ok])
                self.cov[rp] += 1
                self.cov[rp + n] -= 1
                rp += n
                qp += n
            elif op == "I":
                indels.append((rp, n))
                if 0 < rp <= self.L:
                    self.ins[rp - 1] += 1
                    self.ins_len[rp - 1][n] += 1
                qp += n
                mism += n
            elif op == "S":
                qp += n
            elif op == "D":
                indels.append((rp, -n))
                self.dels[rp] += 1
                self.dels[rp + n] -= 1
                self.del_len[rp][n] += 1
                rp += n
                mism += n
            elif op == "N":
                rp += n
        trail = cig and OPS[cig[-1] & 15] in "SH" and cig[-1] >> 4 >= 5
        partners = aux_sa(b, off) if (lead or trail) else []
        for site, clipped in ((pos, lead), (rp, trail)):
            if clipped:
                self.clips[site] += 1
                for p in partners:
                    self.sa[site][p] += 1
        self.starts[pos] += 1
        self.ends[rp - 1] += 1
        mid = min((pos + rp) // 2, self.L - 1)
        if flag & 0x1:
            self.paired += 1
            if flag & 0x8:
                self.mate_unmapped[mid] += 1
            elif mate_rid != refs.index(self.name):
                self.mate_elsewhere[mid] += 1
            elif not flag & 0x2:
                self.improper[mid] += 1
        self.reads.append((pos, rp, np.concatenate(sites) if sites else np.empty(0, int),
                           np.concatenate(bases) if bases else np.empty(0, np.int8), mism, indels))


def report(p, around, window):
    L, ref, counts = p.L, p.ref, p.counts
    tot = counts.sum(0)
    n = tot.sum(1)
    ok = n >= MIN_DEPTH
    cov = np.cumsum(p.cov[:-1])
    dels = np.cumsum(p.dels[:-1])
    major = np.argmax(tot, 1)
    second = np.argsort(tot, 1)[:, -2]
    with np.errstate(invalid="ignore", divide="ignore"):
        maf = tot[np.arange(L), second] / n
        nonmajor = 1 - tot[np.arange(L), major] / n
        nonref = 1 - tot[np.arange(L), np.maximum(ref, 0)] / n
        pi = (n * (n - 1) - (tot * (tot - 1)).sum(1)) / (n * (n - 1))
        del_frac = dels / (cov + dels)
        ins_frac = p.ins / np.maximum(cov, 1)
    reads = p.reads
    mism = np.array([r[4] for r in reads])
    starts = np.array([r[0] for r in reads])
    print(f"== {p.name}: {len(reads)} reads, {L} sites, {ok.sum()} with >= {MIN_DEPTH} Q{MIN_BQ} bases ==\n")
    rl = int(np.median([r[1] - r[0] for r in reads]))
    half = max(rl // 2, 1)
    edge = {}
    for side, at, inner in (("5'", half, slice(2 * rl, 3 * rl)), ("3'", L - 1 - half, slice(L - 3 * rl, L - 2 * rl))):
        edge[side] = (cov[at] / max(cov[inner].mean(), 1), cov[inner].mean(), cov[at])
    # calls within one read length of an end the reads run past are alignment artefacts of the overhang
    near_bad_end = np.zeros(L, bool)
    if 0.25 <= edge["5'"][0] < 0.8:
        near_bad_end[:rl] = True
    if 0.25 <= edge["3'"][0] < 0.8:
        near_bad_end[L - rl:] = True
    mark = lambda i: "*" if near_bad_end[i] else ""

    print("WITHIN-SAMPLE VARIATION (reads against each other)")
    print(f"  pi                         {pi[ok].mean():.5f}"
          f"   (share from sites with MAF < 0.5%: {pi[ok & (maf < 0.005)].sum() / pi[ok].sum():.0%})")
    print(f"  mean MAF                   {maf[ok].mean():.5f}")
    print(f"  mean non-major fraction    {nonmajor[ok].mean():.5f}")
    for t in (0.01, 0.02, 0.05, 0.10, 0.20):
        k = (ok & (maf >= t)).sum()
        print(f"  sites with MAF >= {t:4.0%}     {k:4d}  ({1000 * k / max(ok.sum(), 1):.2f} per kb)")
    edges = [0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
    h, _ = np.histogram(maf[ok & (maf >= 0.01)], bins=edges)
    print("  folded spectrum            " + ", ".join(f"{a:.0%}-{b:.0%}: {x}" for a, b, x in zip(edges, edges[1:], h)))
    isnv = np.where(ok & (maf >= ISNV))[0]
    shown = isnv[np.argsort(-maf[isnv])][:30]
    print(f"  variant sites >= {ISNV:.0%}: {len(isnv)}" + (" (30 highest)" if len(isnv) > 30 else "") + ": "
          + ", ".join(f"{i + 1}{mark(i)} ({maf[i]:.1%})" for i in sorted(shown)))
    if near_bad_end[isnv].any():
        print("    * within one read length of an end the reads run past: alignment artefacts of the overhang, not variation")
    any_alt = sum(bool((bb != major[ss]).any()) for _, _, ss, bb, _, _ in reads if len(ss))
    print(f"  reads with >= 1 non-major Q{MIN_BQ} base: {any_alt / len(reads):.1%}"
          f"  (0.1% error per base over {rl} bp alone: {1 - 0.999 ** rl:.1%})")
    print("  linkage, variant pairs one read can span:")
    for i, j in itertools.combinations(isnv, 2):
        if j - i >= rl:
            continue
        t = np.zeros((2, 2), int)
        for _, _, ss, bb, _, _ in reads:
            if not len(ss) or ss[0] > i or ss[-1] < j:
                continue
            hi, hj = np.searchsorted(ss, i), np.searchsorted(ss, j)
            if hi >= len(ss) or hj >= len(ss) or ss[hi] != i or ss[hj] != j:
                continue
            ai = 0 if bb[hi] == major[i] else 1 if bb[hi] == second[i] else -1
            aj = 0 if bb[hj] == major[j] else 1 if bb[hj] == second[j] else -1
            if ai >= 0 and aj >= 0:
                t[ai, aj] += 1
        m = t.sum()
        if not m:
            continue
        pa, pb = t[1].sum() / m, t[:, 1].sum() / m
        if 0 < pa < 1 and 0 < pb < 1:
            r2 = f"{(t[1, 1] / m - pa * pb) ** 2 / (pa * (1 - pa) * pb * (1 - pb)):.2f}"
        else:
            r2 = "n/a (minor allele absent from spanning reads)"
        print(f"    {i + 1}-{j + 1} ({j - i} nt): {m} spanning reads, r^2 {r2}")
    fwd, rev = counts[0].sum(1), counts[1].sum(1)
    print(f"  fragment orientation: {rev[ok].sum() / n[ok].sum():.1%} of Q{MIN_BQ} bases from reverse fragments")
    print("  minor allele by fragment orientation (read 1's strand; a stranded library's RNA strand):")
    for i in sorted(shown)[:25]:
        f = counts[0][i, second[i]] / max(fwd[i], 1)
        r = counts[1][i, second[i]] / max(rev[i], 1)
        print(f"    {i + 1:>6}: forward {f:6.1%} of {fwd[i]:>6}, reverse {r:6.1%} of {rev[i]:>6}")
    # Indels: a read that ends inside a repeat aligns without the indel, so per-site counts
    # understate it. Judge each cluster of indel sites by the reads spanning all of it.
    ind = np.where(((del_frac >= INDEL) | (ins_frac >= INDEL)) & (cov >= MIN_DEPTH))[0]
    clusters = []
    for i in ind:
        if clusters and i - clusters[-1][1] <= 10:
            clusters[-1][1] = i
        else:
            clusters.append([i, i])
    seq = "".join("ACGTN"[x] for x in ref)

    def repeat_bounds(a, b):
        # widen [a, b) across a tandem repeat: an edge 8-mer that recurs within 30 nt is still inside it
        w0, w1 = max(0, a - 25), min(L, b + 26)
        while w0 > 8 and w0 > a - 125 and seq[w0 - 8:w0] in seq[w0:w0 + 30]:
            w0 -= 1
        while w1 < L - 8 and w1 < b + 125 and seq[w1:w1 + 8] in seq[max(0, w1 - 30):w1]:
            w1 += 1
        return max(0, w0 - 10), min(L, w1 + 10)

    merged = []
    for a, b in clusters:
        w0, w1 = repeat_bounds(a, b + 1)
        if merged and w0 <= merged[-1][3]:
            merged[-1] = [merged[-1][0], b, min(merged[-1][2], w0), max(merged[-1][3], w1)]
        else:
            merged.append([a, b, w0, w1])
    indel_calls = []
    for a, b, w0, w1 in merged:
        net = collections.Counter()
        for r in reads:
            if r[0] <= w0 and r[1] >= w1:
                net[sum(n_ for q, n_ in r[5] if w0 <= q <= w1)] += 1
        m = sum(net.values())
        per_site = max(np.nanmax(del_frac[a:b + 1]), np.nanmax(ins_frac[a:b + 1]))
        indel_calls.append((a, b, w0, w1, net, m, per_site))
    print(f"  indel clusters (sites where >= {INDEL:.0%} of reads carry an indel), judged by reads spanning the cluster and any tandem repeat around it:"
          f" {len(indel_calls)}")
    for a, b, w0, w1, net, m, per_site in indel_calls[:15]:
        alleles = ", ".join(f"{d:+d} nt {c / m:.1%}" for d, c in net.most_common(3)) if m else "no spanning reads"
        print(f"    {a + 1}-{b + 1}{mark(a)}: per site up to {per_site:.1%}; {m} reads span {w0 + 1}-{w1}: {alleles}")

    print(f"\nDIVERGENCE FROM THE REFERENCE (reads against {p.name})")
    snv = np.where(ok & (major != ref))[0]
    maj_indel = [c for c in indel_calls if c[5] >= 20 and c[4].get(0, 0) / c[5] < 0.5]
    print(f"  sites where the read majority is not the reference: {len(snv)} base, {len(maj_indel)} indel")
    for i in snv[:20]:
        print(f"    {i + 1} {'ACGT'[ref[i]]}->{'ACGT'[major[i]]} ({1 - nonmajor[i]:.0%} of {n[i]})")
    for a, b, w0, w1, net, m, _ in maj_indel:
        d, c = net.most_common(1)[0]
        print(f"    {a + 1}-{b + 1}{mark(a)}: {c / m:.0%} of {m} spanning reads are {d:+d} nt against the reference"
              f" (reference length in {net.get(0, 0) / m:.1%})")
    print(f"  mean non-reference fraction {nonref[ok].mean():.5f}"
          f"  (equals the non-major fraction wherever majority = reference)")
    print(f"  mismatches + indel bases per read: mean {mism.mean():.2f}, reads with >= 3: {(mism >= 3).mean():.1%}")
    w = max(100, int(round(L / 9, -2)))
    prof = [mism[(starts >= a) & (starts < a + w)].mean() for a in range(0, L, w)]
    print(f"  mean per read, by {w}-bp window of read start: " + " ".join(f"{x:.2f}" for x in prof))
    with np.errstate(invalid="ignore", divide="ignore"):
        sf, ef = p.starts[:L] / np.maximum(cov, 1), p.ends[:L] / np.maximum(cov, 1)
    for label, arr, frac in (("start", p.starts[:L], sf), ("end", p.ends[:L], ef)):
        hot = [i for i in np.argsort(arr)[::-1][:40] if arr[i] >= 50 and frac[i] >= 0.10][:5]
        if hot:
            print(f"  read-{label} pileups (>= 10% of the depth there {label}s at one base): "
                  + ", ".join(f"{i + 1} ({arr[i]} reads, {frac[i]:.0%})" for i in sorted(hot)))
    clipped = p.clips.sum()
    if clipped:
        with np.errstate(invalid="ignore", divide="ignore"):
            cf = p.clips[:L] / np.maximum(n, 1)
        top = np.argsort(cf * ok)[::-1][:5]
        print(f"  reads soft-clipped >= 5 nt: {clipped} ({clipped / len(reads):.1%});"
              f" sites where >= 5% of reads clip: {(ok & (cf >= 0.05)).sum()}")
        for i in top:
            partners = ", ".join(f"{q} ({c})" for (r, q), c in p.sa[i].most_common(3))
            print(f"    {i + 1}: {p.clips[i]} reads ({cf[i]:.1%})" + (f"; split partners {partners}" if partners else ""))
    else:
        print("  no soft clips (end-to-end alignment): contig ends and mates below carry the weight")
    if p.paired:
        print(f"  pairs: {p.paired} paired reads; mate unmapped {p.mate_unmapped.sum()},"
              f" mate on another reference {p.mate_elsewhere.sum()}, improper {p.improper.sum()}")
        for name, arr in (("mate unmapped", p.mate_unmapped), ("mate elsewhere", p.mate_elsewhere),
                          ("improper", p.improper)):
            if arr.sum():
                top = np.argsort(arr)[::-1][:4]
                bins = np.add.reduceat(arr, np.arange(0, L, 100))
                hot = np.argsort(bins)[::-1][:4]
                print(f"    {name}: busiest 100-nt windows (by read midpoint) "
                      + ", ".join(f"{100 * k + 1}-{min(100 * k + 100, L)} ({bins[k]})" for k in sorted(hot)))
    # With end-to-end alignment a read that would overhang the contig is not aligned, so if the
    # molecule runs past the end, coverage ramps down over one read length (about half the
    # interior depth half a read length in). If the molecule ends there, fragments end there
    # too and coverage holds up to the edge. With a local aligner, the overhang is clipped.
    print(f"  contig ends ({rl}-nt reads; ~50%: the molecule runs past the end; ~100%: it ends here;"
          f" near 0: the reads barely reach it):")
    for side, clip in (("5'", p.clips[0:5].sum()), ("3'", p.clips[L - 5:L + 1].sum())):
        ratio, interior, at = edge[side]
        print(f"    {side}: interior {interior:.0f}x; {half} nt in {at:.0f}x ({ratio:.0%}); reads clipped at the end {int(clip)}")

    if around and (around[0] in (None, p.name) or p.name.endswith(around[0] or "")):
        c, w = around[1], window
        print(f"\n  AROUND {c} ({w}-bp windows; depth = aligned bases; reads = reads starting in the window)")
        print("    window          depth  rev-frag  non-major  maj!=ref  deleted  read starts  read ends  mate-unmapped  mism/read")
        for a in range(max(0, c - 4 * w), min(L, c + 4 * w), w):
            sl = slice(a, min(a + w, L))
            sel = (starts >= a) & (starts < a + w)
            print(f"    {a + 1:>6}-{min(a + w, L):<6} {cov[sl].mean():7.0f}  {rev[sl].sum() / max(n[sl].sum(), 1):7.1%}"
                  f"  {np.nanmean(nonmajor[sl]):9.5f}  {int((ok[sl] & (major[sl] != ref[sl])).sum()):8d}"
                  f"  {np.nanmax(del_frac[sl]):7.1%}  {int(p.starts[sl].sum()):11d}  {int(p.ends[sl].sum()):9d}"
                  f"  {int(p.mate_unmapped[sl].sum()):13d}  {mism[sel].mean() if sel.any() else float('nan'):9.2f}")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", default=os.path.join(EXAMPLES, "glrav13.bam"))
    ap.add_argument("--fasta", default=os.path.join(EXAMPLES, "glrav13.fa"))
    ap.add_argument("--around", metavar="[REF:]POS",
                    help="profile depth, variation and divergence in windows either side of POS")
    ap.add_argument("--window", type=int, default=250, help="window size for --around [%(default)s]")
    args = ap.parse_args()
    around = None
    if args.around:
        r, _, q = args.around.rpartition(":")
        around = (r or None, int(q))
    fasta = read_fasta(args.fasta)
    refs, records = bam_open(args.bam)
    names = [r for r, _ in refs]
    piles = {}
    for name, length in refs:
        if name not in fasta:
            raise SystemExit(f"{name} is in the BAM header but not in {args.fasta}")
        if len(fasta[name]) != length:
            raise SystemExit(f"{name}: {len(fasta[name])} nt in the FASTA, {length} in the BAM header")
        piles[name] = Pile(name, fasta[name])
    for b in records:
        (rid,) = struct.unpack_from("<i", b, 0)
        (flag,) = struct.unpack_from("<H", b, 14)
        if flag & SKIP or rid < 0:
            continue
        piles[names[rid]].add(b, names)
    print(f"{os.path.basename(args.bam)}\n")
    for p in piles.values():
        if p.reads:
            report(p, around, args.window)


if __name__ == "__main__":
    main()
