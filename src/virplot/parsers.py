"""GFF3 and depth-file parsers, plus the dispatcher that turns any depth
source (samtools depth text, SAM, BAM) into a per-base array."""

from __future__ import annotations

import logging
import os
import re
from urllib.parse import unquote

import numpy as np

from virplot.alignments import depth_from_alignments, is_alignment_file, resolve_reference
from virplot.models import RNA, Domain, Feature, Noncoding
from virplot.settings import classify_function

log = logging.getLogger(__name__)


class ParseError(ValueError):
    """A GFF3 or depth file VirPlot cannot use; the message says where and why."""


# --------------------------------------------------------------------------
# GFF3
# --------------------------------------------------------------------------

# Row types read as polyprotein domains: NCBI's GFF3 spelling of GenBank's
# mat_peptide, plus the GenBank name itself and the SO synonym, which other
# converters write.
_DOMAIN_TYPES = {"mature_protein_region_of_CDS", "mat_peptide", "mature_protein_region"}

# Non-coding landmarks (ICTV: "IR, UTR, stem-loop — grey or black", drawn on
# the genome line). Two tiers: types that are non-coding by definition are
# always read; catch-all types (GenBank misc_feature, which NCBI's GFF3 calls
# sequence_feature, regulatory, repeat_region) are read only when their name
# or note says intergenic / common region / UTR / stem-loop etc., because
# RefSeq also uses them for motifs inside ORFs that would clutter the line.
_NONCODING_ALWAYS = {
    "five_prime_UTR": ("region", "5\u2032 UTR"),
    "three_prime_UTR": ("region", "3\u2032 UTR"),
    "UTR": ("region", "UTR"),
    "intergenic_region": ("region", "IR"),
    "origin_of_replication": ("region", "ori"),
    "stem_loop": ("stem_loop", ""),
}
_NONCODING_IF_NAMED = {"sequence_feature", "misc_feature", "regulatory_region", "regulatory",
                       "repeat_region", "sequence_secondary_structure", "misc_structure"}
_STEM_LOOP_WORDS = re.compile(r"stem[- ]?loop|hairpin", re.I)


def parse_gff_rnas(gff_path: str) -> list[RNA]:
    """Parse a GFF3 file into one ``RNA`` per ``region`` line.

    RNAs are returned in the order their ``region`` lines appear. CDS features
    are attached to the RNA whose sequence id (column 1) they carry; a CDS on
    a sequence id with no ``region`` line is skipped with a warning. Mature
    protein rows become ``Domain`` segments of the CDS they belong to.
    """
    rnas: dict[str, RNA] = {}
    rows: dict[str, list[tuple[int, int, str, dict]]] = {}
    domain_rows: dict[str, list[tuple[int, int, str, dict]]] = {}
    noncoding: dict[str, list[Noncoding]] = {}

    # Collect everything first, attach afterwards: GFF3 does not require the
    # region line to precede the features of its sequence.
    with open(gff_path) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            seqid, _, feature_type, start, end, _, strand, _, attributes = parts

            if feature_type == "region":
                if seqid not in rnas:
                    attrs = _attributes(attributes)
                    circular = attrs.get("Is_circular", "").lower() == "true"
                    rnas[seqid] = RNA(name=seqid, seqid=seqid, length=int(end),
                                      circular=circular)
                continue

            is_nc = feature_type in _NONCODING_ALWAYS or feature_type in _NONCODING_IF_NAMED
            if feature_type != "CDS" and feature_type not in _DOMAIN_TYPES and not is_nc:
                continue

            info = _attributes(attributes)
            if is_nc:
                nc = _noncoding(feature_type, int(start), int(end), strand, info)
                if nc is not None:
                    noncoding.setdefault(seqid, []).append(nc)
                continue
            target = rows if feature_type == "CDS" else domain_rows
            target.setdefault(seqid, []).append((int(start), int(end), strand, info))

    if not rnas:
        raise ParseError(f"No 'region' feature found in GFF: {gff_path}")

    for seqid, rna in rnas.items():
        rna.features.extend(_features_from_rows(rows.get(seqid, []),
                                                domain_rows.get(seqid, [])))
        rna.noncoding.extend(noncoding.get(seqid, []))
        if not rna.circular:
            wrapped = [f for f in rna.features if rna.wraps(f)]
            if wrapped:
                log.warning("%s is not marked circular but %d feature(s) cross the origin "
                            "(%s); add Is_circular=true to its region line or --topology circular",
                            seqid, len(wrapped), ", ".join(f.product for f in wrapped[:3]))

    for seqid in sorted(set(rows) - set(rnas)):
        log.warning("CDS features on %r ignored: no 'region' line for that sequence", seqid)

    return list(rnas.values())


def _attributes(column: str) -> dict[str, str]:
    """GFF3 column 9 as a dict, with ``%XX`` escapes decoded.

    GFF3 reserves ``;``, ``=``, ``&`` and ``,`` inside values and writes them
    percent-encoded; NCBI does so routinely (``product=polyprotein%2C putative``,
    ``transl_except=(pos:3417..3419%2Caa:OTHER)``). Decoding here means labels
    read naturally and the function-word classifier sees real text.
    """
    out: dict[str, str] = {}
    for kv in column.split(";"):
        if "=" in kv:
            key, value = kv.split("=", 1)
            out[key.strip()] = unquote(value)
    return out


# --- expression mechanisms ---------------------------------------------------
#
# GFF3 has no field for "this ORF is reached by a frameshift / by reading
# through a stop". RefSeq encodes a slipped CDS as several rows sharing one ID
# with ``exception=ribosomal slippage``, and a readthrough product as a CDS
# carrying ``transl_except=``. VirPlot reads both, plus its own ``Note=``
# convention (``Note=+1 frameshift``, ``Note=readthrough``) for hand-written
# files. See docs/GFF_GUIDE.md section 8.

_SHIFT_IN_NOTE = re.compile(r"([+-])\s*([12])\s*(?:ribosomal\s+)?frameshift", re.I)


def _features_from_rows(rows: list[tuple[int, int, str, dict]],
                        domain_rows: list[tuple[int, int, str, dict]] = ()) -> list[Feature]:
    """Turn one molecule's CDS rows into Features, resolving joins and mechanisms,
    then hang any mature-protein rows off the CDS they belong to."""
    # 1. group rows that share an ID (a RefSeq join); keep first-seen order
    groups: dict[str, list[tuple[int, int, str, dict]]] = {}
    order: list[str] = []
    for i, row in enumerate(rows):
        key = row[3].get("ID") or f"__row{i}"
        if key not in groups:
            groups[key] = []; order.append(key)
        groups[key].append(row)

    feats: list[Feature] = []
    joined: set[int] = set()         # first segments of slipped joins
    for key in order:
        segs = groups[key]
        if len(segs) == 1:
            feats.append(_feature(*segs[0]))
            continue
        info = segs[0][3]
        slipped = ("slippage" in info.get("exception", "").lower()
                   or "frameshift" in info.get("Note", "").lower())
        if not slipped:
            log.info("CDS %s has %d rows without a ribosomal-slippage note; drawn as %d boxes",
                     key, len(segs), len(segs))
            feats.extend(_feature(*seg) for seg in segs)
            continue
        # translation order: ascending on +, descending on -
        forward = segs[0][2] != "-"
        segs = sorted(segs, key=lambda r: r[0], reverse=not forward)
        first = _feature(*segs[0])
        feats.append(first)
        joined.add(id(first))        # its Note describes the join, not itself
        prev = first
        for seg in segs[1:]:
            shift = _shift_between(prev, seg[0], seg[1], forward)
            cont = _feature(*seg, mechanism="frameshift", shift=shift,
                            show_label=seg[3].get("product") != segs[0][3].get("product"))
            feats.append(cont)
            prev = cont

    # 2. single rows that declare a mechanism in their own attributes
    resolved: list[Feature] = []
    for f in feats:
        if f.mechanism is not None or id(f) in joined:
            resolved.append(f); continue
        info = f._attrs                                    # stashed by _feature
        note = info.get("Note", "").lower()
        if "readthrough" in note or "read-through" in note or "transl_except" in info:
            resolved.append(_as_readthrough(f, feats))
        elif "frameshift" in note:
            m = _SHIFT_IN_NOTE.search(info.get("Note", ""))
            shift = int(m.group(1) + m.group(2)) if m else _shift_from_neighbour(f, feats)
            resolved.append(_with(f, mechanism="frameshift", shift=shift))
        else:
            resolved.append(f)

    # 3. polyprotein domains
    resolved = _attach_domains(resolved, domain_rows)
    return [_strip_attrs(f) for f in resolved]


# --- polyprotein domains -------------------------------------------------------
#
# RefSeq writes the mature proteins of a polyprotein as
# ``mature_protein_region_of_CDS`` rows (GenBank ``mat_peptide``) with
# ``Parent=cds-…`` naming the CDS and ``product=`` naming the protein. The ICTV
# figures draw them as one box split by thin lines, so they become ``Domain``
# segments of the CDS's Feature rather than Features of their own.

def _attach_domains(feats: list[Feature],
                    domain_rows: list[tuple[int, int, str, dict]]) -> list[Feature]:
    """Give each Feature the domains that belong to it; warn about any that fit nowhere.

    Matching: ``Parent=`` against the CDS ``ID`` first; failing that, the
    smallest CDS on the same strand that contains the domain (for files
    written without ``Parent``). A joined CDS is several Features sharing
    one ID; the domain goes to the segment it overlaps most.

    LIMIT: a domain that straddles a frameshift junction (coronavirus nsp12
    begins in ORF1a and ends in ORF1b) is attached to one segment only and
    drawn clipped to it, so its far end carries no divider.
    """
    if not domain_rows:
        return feats
    by_id: dict[str, list[int]] = {}
    for i, f in enumerate(feats):
        fid = f._attrs.get("ID")
        if fid:
            by_id.setdefault(fid, []).append(i)

    attached: dict[int, list[Domain]] = {}
    for start, end, strand, info in domain_rows:
        product = info.get("product", "unknown")
        candidates: list[int] = []
        for parent in info.get("Parent", "").split(","):
            candidates.extend(by_id.get(parent.strip(), []))
        if not candidates:
            mid = (start + end) / 2
            candidates = [i for i, f in enumerate(feats)
                          if f.strand == strand and f.start <= mid <= f.end]
            candidates.sort(key=lambda i: feats[i].length)       # smallest container
            candidates = candidates[:1]
        if not candidates:
            log.warning("Mature protein %r at %d..%d matches no CDS; ignored", product, start, end)
            continue
        best = max(candidates, key=lambda i: min(end, feats[i].end) - max(start, feats[i].start))
        attached.setdefault(best, []).append(Domain(start, end, product))

    return [_with(f, domains=tuple(sorted(attached[i], key=lambda d: d.start)))
            if i in attached else f
            for i, f in enumerate(feats)]


def _feature(start: int, end: int, strand: str, info: dict, **extra) -> Feature:
    f = Feature(start=start, end=end, strand=strand,
                product=info.get("product", "unknown"),
                gene=info.get("gene") or None, **extra)
    object.__setattr__(f, "_attrs", info)                # frozen dataclass: side channel
    return f


def _with(f: Feature, **changes) -> Feature:
    g = Feature(**{**{k: getattr(f, k) for k in Feature.__dataclass_fields__}, **changes})
    object.__setattr__(g, "_attrs", getattr(f, "_attrs", {}))
    return g


def _strip_attrs(f: Feature) -> Feature:
    if hasattr(f, "_attrs"):
        object.__delattr__(f, "_attrs")
    return f


def _shift_between(prev: Feature, start: int, end: int, forward: bool) -> int | None:
    """Frameshift sign from the join geometry: skipping 1 nt is +1, re-reading 1 nt is -1.
    None (drawn as a bare ``FS``) when the segments are in frame with each other.

    LIMIT: this is inference, not annotation. It is exactly right for a RefSeq
    join, whose segment boundaries are the slippage site, but a hand-written
    pair whose coordinates do not reflect the slip (say, ORF1b written from
    its first full codon) will get a wrong or zero sign. A ``Note=`` with an
    explicit ``+1``/``-1`` overrides this (see ``_features_from_rows``);
    docs/GFF_GUIDE.md section 8 tells authors to write one in that case.
    """
    gap = (start - prev.end - 1) if forward else (prev.start - end - 1)
    return {1: 1, 2: -1}.get(gap % 3, None)


def _shift_from_neighbour(f: Feature, feats: list[Feature]) -> int | None:
    """For a Note=frameshift row with no sign given: infer it from the nearest upstream ORF.

    Same LIMIT as ``_shift_between``: correct only when the two ORFs' coordinates
    meet at the slippage site. Returns None (drawn as a bare ``FS``) when no
    upstream ORF is found or the gap is a multiple of three.
    """
    same = [g for g in feats if g is not f and g.strand == f.strand]
    if f.forward:
        ups = [g for g in same if g.end <= f.start + 2]
        if not ups:
            return None
        prev = max(ups, key=lambda g: g.end)
        gap = f.start - prev.end - 1
    else:
        ups = [g for g in same if g.start >= f.end - 2]
        if not ups:
            return None
        prev = min(ups, key=lambda g: g.start)
        gap = prev.start - f.end - 1
    return {1: 1, 2: -1}.get(gap % 3, None)


def _as_readthrough(f: Feature, feats: list[Feature]) -> Feature:
    """A readthrough product spans the ORF it extends; draw only the extension.

    RefSeq annotates e.g. TMV 183K as 69..4919 beside 126K at 69..3419. The
    ICTV figures show the readthrough as a second box abutting the first, so
    the shared part is trimmed off and the bar marks the read-through stop.
    """
    for g in feats:
        if g is f or g.strand != f.strand or g.mechanism is not None:
            continue
        if f.forward and g.start == f.start and g.end < f.end:
            return _with(f, start=g.end + 1, mechanism="readthrough")
        if not f.forward and g.end == f.end and g.start > f.start:
            return _with(f, end=g.start - 1, mechanism="readthrough")
    return _with(f, mechanism="readthrough")


# --- non-coding landmarks --------------------------------------------------------

def _noncoding(ftype: str, start: int, end: int, strand: str, info: dict) -> Noncoding | None:
    """Build a ``Noncoding`` from a row, or None if a catch-all row is not one.

    The label is the first of ``Name``, ``product``, ``gene``,
    ``standard_name``, ``regulatory_class``, then ``Note``; for the dedicated
    types a default (``5′ UTR``, ``IR``, ``ori``) stands in when none is
    given; a ``stem_loop`` row is labelled only by its own ``Name`` /
    ``standard_name`` (not a ``gene=`` RefSeq may attach to it). A short
    region whose words say stem-loop / hairpin is drawn as one.
    """
    named = next((info[k] for k in ("Name", "product", "gene", "standard_name",
                                    "regulatory_class") if info.get(k)), "")
    note = info.get("Note", "")
    text = named or note
    if ftype in _NONCODING_ALWAYS:
        kind, default = _NONCODING_ALWAYS[ftype]
        if kind == "stem_loop":
            # a stem-loop's Note is a description, and RefSeq hangs the
            # neighbouring ORF's gene= on it (FBNYV DNA-U4): only an explicit
            # Name / standard_name is the hairpin's own name
            label = next((info[k] for k in ("Name", "standard_name") if info.get(k)), "")
        else:
            label = text or default
    else:
        if classify_function(f"{text} {note}") != "noncoding":
            return None                      # a motif, promoter, repeat inside an ORF…
        kind, label = "region", text
    if kind == "region" and _STEM_LOOP_WORDS.search(f"{text} {note}") and (end - start) < 100:
        kind = "stem_loop"
    return Noncoding(start, end, strand, label, kind)


def parse_gff(gff_path: str) -> tuple[int, list[Feature]]:
    """Parse a GFF3 file, returning (sequence_length, features) of its first RNA."""
    first = parse_gff_rnas(gff_path)[0]
    return first.length, first.features


def gff_seqid(gff_path: str) -> str | None:
    """Return the sequence id (column 1) of the GFF's first ``region`` line, if any."""
    with open(gff_path) as fp:
        for line in fp:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 9 and parts[2] == "region":
                return parts[0]
    return None


# --------------------------------------------------------------------------
# samtools depth text
# --------------------------------------------------------------------------

def parse_depth(depth_path: str, seq_len: int, *, seqid: str | None = None,
                ref: str | None = None, circular: bool = False) -> tuple[np.ndarray, int]:
    """Parse a samtools depth file directly into a numpy array.

    With ``seqid``/``ref`` the file may hold several sequences and the
    matching one is selected (see ``resolve_reference``); without either,
    every line is used, as before. With ``circular``, positions past the end
    (as a padded reference produces) wrap round instead of being dropped.
    Returns (depth_array, n_entries).
    """
    by_seq: dict[str, list[tuple[int, int]]] = {}
    with open(depth_path) as fp:
        for lineno, line in enumerate(fp, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ParseError(f"Malformed depth line in {depth_path}:{lineno} — expected 3 "
                                 f"tab-separated columns, got {len(fields)}")
            try:
                pos, cov = int(fields[1]), int(fields[2])
            except ValueError:
                raise ParseError(f"Non-integer value in {depth_path}:{lineno} — "
                                 f"{line.rstrip()!r}") from None
            by_seq.setdefault(fields[0], []).append((pos, cov))

    if seqid is None and ref is None:
        rows = [r for rows in by_seq.values() for r in rows]
    elif not by_seq:
        rows = []
    else:
        target = resolve_reference(by_seq, seqid=seqid, ref=ref, path=depth_path)
        rows = by_seq.get(target, [])

    y = np.zeros(seq_len, dtype=int)
    for pos, cov in rows:
        if circular and pos >= 1:
            y[(pos - 1) % seq_len] += cov
        elif 1 <= pos <= seq_len:
            y[pos - 1] = cov
    return y, len(rows)


# --------------------------------------------------------------------------
# dispatcher
# --------------------------------------------------------------------------

def load_depth(
    path: str,
    seq_len: int,
    *,
    seqid: str | None = None,
    ref: str | None = None,
    min_mapq: int = 0,
    circular: bool = False,
) -> tuple[np.ndarray, int, str]:
    """Load one depth track from a depth file or a SAM/BAM file.

    Returns ``(depth_array, n, kind)`` where ``n`` is the number of depth
    entries read (text) or reads counted (alignments) and ``kind`` is
    ``"depth"`` or ``"alignments"``.
    """
    if is_alignment_file(path):
        y, n = depth_from_alignments(path, seq_len, seqid=seqid, ref=ref,
                                     min_mapq=min_mapq, circular=circular)
        return y, n, "alignments"
    y, n = parse_depth(path, seq_len, seqid=seqid, ref=ref, circular=circular)
    return y, n, "depth"


def default_label(path: str) -> str:
    """Label for a depth source: its basename without extension."""
    return os.path.splitext(os.path.basename(path))[0]
