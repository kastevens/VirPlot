"""Measurements behind docs/variation_design.md. Not part of VirPlot.

Reads a BAM and the FASTA it was mapped to (default: the shipped GLRaV-13
example) and prints two kinds of statistic, kept apart on purpose:

* within-sample variation: how the reads at one site disagree with each
  other (minor-allele frequency, pi, the allele-frequency spectrum, linkage);
* divergence from the reference: how the reads disagree with the sequence
  they were mapped to (sites where the read majority is not the reference
  base, soft-clip pileups, mismatches per read).

    python3 docs/design/variation_measurements.py
    python3 docs/design/variation_measurements.py --bam S.bam --fasta ref.fa --around 12239

--around POS adds a profile of both kinds of signal in windows either side of
POS, for asking what a depth step there is made of. Single-reference BAMs
only. Soft clips appear only with a local aligner (bwa mem, bowtie2 --local).
Standard library and numpy only. Bases count at base quality >= 30; sites
count with >= 100 such bases.
"""

import argparse
import gzip
import itertools
import os
import struct

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.join(HERE, "..", "..", "examples")
MIN_BQ, MIN_DEPTH, ISNV = 30, 100, 0.02
SKIP = 0x4 | 0x100 | 0x200 | 0x400 | 0x800       # unmapped, secondary, QC fail, duplicate, supplementary

SEQ = "=ACMGRSVTWYHKDBN"
CODE = np.full(16, -1, np.int8)
for i, b in enumerate("ACGT"):
    CODE[SEQ.index(b)] = i
OPS = "MIDNSHP=X"


def read_fasta(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    return np.array(["ACGT".find(c) for c in seq], np.int8)


def bam_records(path):
    fp = gzip.open(path, "rb")
    assert fp.read(4) == b"BAM\x01"
    (lt,) = struct.unpack("<i", fp.read(4))
    fp.read(lt)
    (nref,) = struct.unpack("<i", fp.read(4))
    for _ in range(nref):
        (ln,) = struct.unpack("<i", fp.read(4))
        fp.read(ln + 4)
    while True:
        head = fp.read(4)
        if len(head) < 4:
            return
        (size,) = struct.unpack("<i", head)
        yield fp.read(size)


def pileup(bam, ref):
    """Allele counts by read orientation, plus per-read records and clip pileups."""
    L = len(ref)
    counts = np.zeros((2, L, 4), np.int64)        # [orientation 0 fwd / 1 rev, site, base]
    clips = np.zeros(L + 1, np.int64)              # reads soft-clipped (>= 5 nt) at each site
    reads = []                                      # (start, end, sites, bases, mismatches)
    for b in bam_records(bam):
        _ref_id, pos, l_rn, _mapq, _bin, n_cig, flag, l_seq = struct.unpack_from("<iiBBHHHi", b, 0)
        if flag & SKIP:
            continue
        off = 32 + l_rn
        cig = struct.unpack_from(f"<{n_cig}I", b, off)
        off += 4 * n_cig
        packed = np.frombuffer(b, np.uint8, (l_seq + 1) // 2, off)
        off += (l_seq + 1) // 2
        nib = np.empty(2 * len(packed), np.uint8)
        nib[0::2], nib[1::2] = packed >> 4, packed & 15
        base, qual = CODE[nib[:l_seq]], np.frombuffer(b, np.uint8, l_seq, off)
        s = 1 if flag & 0x10 else 0
        if cig and OPS[cig[0] & 15] == "S" and cig[0] >> 4 >= 5:
            clips[pos] += 1
        rp, qp, sites, bases, mism = pos, 0, [], [], 0
        for c in cig:
            n, op = c >> 4, OPS[c & 15]
            if op in "M=X":
                bb, qq = base[qp:qp + n], qual[qp:qp + n]
                idx = np.arange(rp, rp + n)
                mism += int(((bb >= 0) & (bb != ref[idx])).sum())
                ok = (bb >= 0) & (qq >= MIN_BQ)
                np.add.at(counts[s], (idx[ok], bb[ok]), 1)
                sites.append(idx[ok])
                bases.append(bb[ok])
                rp += n
                qp += n
            elif op in "IS":
                qp += n
                mism += n if op == "I" else 0
            elif op in "DN":
                rp += n
                mism += n if op == "D" else 0
        if cig and OPS[cig[-1] & 15] == "S" and cig[-1] >> 4 >= 5:
            clips[rp] += 1
        reads.append((pos, rp, np.concatenate(sites) if sites else np.empty(0, int),
                      np.concatenate(bases) if bases else np.empty(0, np.int8), mism))
    return counts, clips[:L], reads


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", default=os.path.join(EXAMPLES, "glrav13.bam"))
    ap.add_argument("--fasta", default=os.path.join(EXAMPLES, "glrav13.fa"))
    ap.add_argument("--around", type=int, metavar="POS",
                    help="profile depth, variation and divergence in windows either side of POS")
    ap.add_argument("--window", type=int, default=250, help="window size for --around [%(default)s]")
    args = ap.parse_args()
    ref = read_fasta(args.fasta)
    counts, clips, reads = pileup(args.bam, ref)
    tot = counts.sum(0)
    n = tot.sum(1)
    ok = n >= MIN_DEPTH
    L = len(ref)
    major = np.argmax(tot, 1)
    order = np.argsort(tot, 1)
    second = order[:, -2]
    with np.errstate(invalid="ignore", divide="ignore"):
        maf = tot[np.arange(L), second] / n                      # second allele
        nonmajor = 1 - tot[np.arange(L), major] / n              # within-sample
        nonref = 1 - tot[np.arange(L), ref] / n                  # against the reference
        pi = (n * (n - 1) - (tot * (tot - 1)).sum(1)) / (n * (n - 1))
    print(f"{os.path.basename(args.bam)}: {len(reads)} reads, {L} sites, "
          f"{ok.sum()} with >= {MIN_DEPTH} Q{MIN_BQ} bases\n")

    print("WITHIN-SAMPLE VARIATION (reads against each other)")
    print(f"  pi                         {pi[ok].mean():.5f}"
          f"   (share from sites with MAF < 0.5%: {pi[ok & (maf < 0.005)].sum() / pi[ok].sum():.0%})")
    print(f"  mean MAF                   {maf[ok].mean():.5f}")
    print(f"  mean non-major fraction    {nonmajor[ok].mean():.5f}")
    for t in (0.01, 0.02, 0.05, 0.10, 0.20):
        k = (ok & (maf >= t)).sum()
        print(f"  sites with MAF >= {t:4.0%}     {k:4d}  ({1000 * k / ok.sum():.2f} per kb)")
    edges = [0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
    h, _ = np.histogram(maf[ok & (maf >= 0.01)], bins=edges)
    print("  folded spectrum            " + ", ".join(f"{a:.0%}-{b:.0%}: {x}" for a, b, x in zip(edges, edges[1:], h)))
    isnv = np.where(ok & (maf >= ISNV))[0]
    print(f"  variant sites >= {ISNV:.0%}:     " + ", ".join(f"{i + 1} ({maf[i]:.1%})" for i in isnv))
    any_alt = sum(bool((bb != major[ss]).any()) for _, _, ss, bb, _ in reads if len(ss))
    print(f"  reads with >= 1 non-major Q{MIN_BQ} base: {any_alt / len(reads):.1%}"
          f"  (0.1% error per base over 100 bp alone: {1 - 0.999 ** 100:.1%})")
    print("  linkage, variant pairs one read can span:")
    for i, j in itertools.combinations(isnv, 2):
        if j - i >= 100:
            continue
        t = np.zeros((2, 2), int)
        for _, _, ss, bb, _ in reads:
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
        pa, pb = t[1].sum() / m, t[:, 1].sum() / m
        if 0 < pa < 1 and 0 < pb < 1:
            r2 = f"{(t[1, 1] / m - pa * pb) ** 2 / (pa * (1 - pa) * pb * (1 - pb)):.2f}"
        else:
            r2 = "n/a (minor allele absent from spanning reads)"
        print(f"    {i + 1}-{j + 1} ({j - i} nt): {m} spanning reads, r^2 {r2}")
    fwd, rev = counts[0].sum(1), counts[1].sum(1)
    print("  minor allele by read orientation (reverse = + strand RNA in this dUTP library):")
    for i in isnv:
        f = counts[0][i, second[i]] / max(fwd[i], 1)
        r = counts[1][i, second[i]] / max(rev[i], 1)
        print(f"    {i + 1:>6}: forward reads {f:6.1%} of {fwd[i]:>4}, reverse reads {r:6.1%} of {rev[i]:>5}")

    print("\nDIVERGENCE FROM THE REFERENCE (reads against glrav13.fa)")
    disagree = np.where(ok & (major != ref))[0]
    print(f"  sites where the read majority is not the reference base: {len(disagree)}"
          + (": " + ", ".join(f"{i + 1} ({'ACGT'[ref[i]]}->{'ACGT'[major[i]]}, "
                               f"{1 - nonmajor[i]:.0%})" for i in disagree) if len(disagree) else ""))
    print(f"  mean non-reference fraction {nonref[ok].mean():.5f}"
          f"  (equals the non-major fraction wherever majority = reference)")
    with np.errstate(invalid="ignore", divide="ignore"):
        clip_frac = clips / np.maximum(n, 1)
    top = np.argsort(clip_frac * ok)[::-1][:5]
    print(f"  reads soft-clipped >= 5 nt: {clips.sum()} ({clips.sum() / len(reads):.1%} of reads);"
          f" sites where >= 5% of reads clip: {(ok & (clip_frac >= 0.05)).sum()}")
    print("  highest clip pileups: " + ", ".join(f"{i + 1} ({clips[i]} reads, {clip_frac[i]:.1%})" for i in top))
    mism = np.array([r[4] for r in reads])
    print(f"  mismatches + indel bases per read: mean {mism.mean():.2f}, reads with >= 3: {(mism >= 3).mean():.1%}")
    starts = np.array([r[0] for r in reads])
    w = 2000
    prof = [mism[(starts >= a) & (starts < a + w)].mean() for a in range(0, L, w)]
    print(f"  mean per read, by {w}-bp window of read start: " + " ".join(f"{p:.2f}" for p in prof))

    if args.around:
        c, w = args.around, args.window
        print(f"\nAROUND {c} ({w}-bp windows; depth = Q{MIN_BQ} bases; reads = reads starting in the window)")
        print("  window          depth  rev-read  non-major  maj!=ref  clipped  mism/read")
        for a in range(max(0, c - 4 * w), min(L, c + 4 * w), w):
            sl = slice(a, min(a + w, L))
            sel = (starts >= a) & (starts < a + w)
            print(f"  {a + 1:>6}-{min(a + w, L):<6} {n[sl].mean():7.0f}  {rev[sl].sum() / max(n[sl].sum(), 1):7.1%}"
                  f"  {np.nanmean(nonmajor[sl]):9.5f}  {int((ok[sl] & (major[sl] != ref[sl])).sum()):8d}"
                  f"  {int(clips[sl].sum()):7d}  {mism[sel].mean() if sel.any() else float('nan'):9.2f}")


if __name__ == "__main__":
    main()
