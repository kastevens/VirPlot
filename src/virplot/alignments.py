"""Per-base depth straight from SAM/BAM alignments.

This replaces the ``samtools depth -a`` step: VirPlot walks each alignment's
CIGAR string and counts, for every reference position, how many reads cover
it.  Defaults reproduce ``samtools depth -a`` as closely as practical:

* reads flagged unmapped, secondary, QC-fail or duplicate are skipped
  (samtools' default ``-G UNMAP,SECONDARY,QCFAIL,DUP``);
* only aligned bases count (CIGAR ``M``, ``=``, ``X``); deletions (``D``) and
  reference skips (``N``) do not, matching samtools without ``-J``;
* no mapping-quality or base-quality threshold unless ``min_mapq`` is set.

Both SAM (text) and BAM (BGZF-compressed binary) are read with the standard
library only — BGZF is a valid multi-member gzip stream, so :mod:`gzip`
decodes it.  CRAM is not supported.
"""

from __future__ import annotations

import gzip
import logging
import os
import re
import struct
from dataclasses import dataclass
from typing import Iterator

import numpy as np

log = logging.getLogger(__name__)

# SAM flag bits skipped by default (samtools depth: UNMAP,SECONDARY,QCFAIL,DUP)
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_BAM_CIGAR_OPS = "MIDNSHP=X"
_REF_CONSUMING = {"M", "D", "N", "=", "X"}
_DEPTH_OPS = {"M", "=", "X"}


class AlignmentError(ValueError):
    """Raised for malformed input or an unresolvable reference name."""


@dataclass(frozen=True)
class Alignment:
    """The subset of an alignment record that depth counting needs."""

    ref: str            # reference name (already resolved for BAM)
    pos: int            # 1-based leftmost mapping position
    flag: int
    mapq: int
    cigar: tuple[tuple[int, str], ...]   # ((length, op), ...)


@dataclass
class AlignmentFile:
    """Header references plus a lazy iterator over alignment records."""

    path: str
    references: dict[str, int]      # name -> length from @SQ
    records: Iterator[Alignment]


# --------------------------------------------------------------------------
# public API
# --------------------------------------------------------------------------

def is_alignment_file(path: str) -> bool:
    """True if ``path`` looks like SAM or BAM (by extension, then content)."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".sam", ".bam"):
        return True
    if ext in (".dep", ".depth", ".txt", ".tsv"):
        return False
    with open(path, "rb") as fp:
        head = fp.read(4)
    if head[:2] == b"\x1f\x8b":            # gzip magic -> assume BAM
        return True
    return head[:1] == b"@"                # SAM header line


def open_alignments(path: str) -> AlignmentFile:
    """Open a SAM or BAM file; the record iterator is consumed lazily."""
    with open(path, "rb") as fp:
        magic = fp.read(2)
    if magic == b"\x1f\x8b":
        return _open_bam(path)
    return _open_sam(path)


def depth_from_alignments(
    path: str,
    seq_len: int,
    *,
    seqid: str | None = None,
    ref: str | None = None,
    min_mapq: int = 0,
    circular: bool = False,
) -> tuple[np.ndarray, int]:
    """Compute per-base depth for one reference from a SAM/BAM file.

    ``seqid`` is the GFF sequence id used to choose among ``@SQ`` entries when
    the file has several; ``ref`` overrides that choice explicitly.

    With ``circular``, alignments that run off the end of the reference wrap
    round to position 1 instead of being clipped, and positions past the end
    are taken modulo the length. Aligners have no notion of a circular
    reference, so this is what makes depth continuous across the origin for a
    BAM produced by padding the reference and wrapping the coordinates back.

    Returns ``(depth_array of length seq_len, n_reads_counted)``.
    """
    af = open_alignments(path)
    target = resolve_reference(af.references, seqid=seqid, ref=ref, path=path)

    ref_len = af.references.get(target)
    if ref_len is not None and ref_len != seq_len and not (circular and ref_len > seq_len):
        log.warning(
            "Reference %s is %d bp in %s but %d bp in the GFF; positions beyond "
            "the GFF length are ignored", target, ref_len, os.path.basename(path), seq_len,
        )

    # difference-array trick: +1 at start, -1 after end, cumulative sum at the end
    diff = np.zeros(seq_len + 1, dtype=np.int64)
    counted = 0
    for aln in af.records:
        if aln.ref != target or aln.flag & SKIP_FLAGS or aln.mapq < min_mapq:
            continue
        counted += 1
        pos = aln.pos - 1                  # 0-based
        for length, op in aln.cigar:
            if op in _DEPTH_OPS:
                if circular:
                    _add_wrapped(diff, pos, length, seq_len)
                else:
                    s, e = pos, pos + length
                    if s < seq_len and e > 0:
                        diff[max(s, 0)] += 1
                        diff[min(e, seq_len)] -= 1
            if op in _REF_CONSUMING:
                pos += length

    depth = np.cumsum(diff[:-1]).astype(int)
    return depth, counted


def _add_wrapped(diff: np.ndarray, start: int, length: int, seq_len: int) -> None:
    """Mark ``length`` covered bases from ``start`` on a circular reference.

    ``start`` is 0-based and may sit past the end (an unwrapped padded
    alignment); the block may also run off the end, in which case it
    continues from position 1.
    """
    if length <= 0 or seq_len <= 0:
        return
    if length >= seq_len:                  # covers the whole circle
        diff[0] += 1
        diff[seq_len] -= 1
        return
    s = start % seq_len
    e = s + length
    if e <= seq_len:
        diff[s] += 1
        diff[e] -= 1
    else:                                  # split across the origin
        diff[s] += 1
        diff[seq_len] -= 1
        diff[0] += 1
        diff[e - seq_len] -= 1


# --------------------------------------------------------------------------
# reference selection
# --------------------------------------------------------------------------

def resolve_reference(references: dict, *, seqid: str | None,
                      ref: str | None, path: str) -> str:
    """Pick which sequence to use from a file's reference names.

    Priority: explicit ``ref``; the GFF ``seqid`` when present in the file (or
    when the file lists no references at all); the only reference if there is
    exactly one. Anything else is ambiguous and raises ``AlignmentError``.
    ``references`` may be any mapping or set of names.
    """
    if ref is not None:
        if references and ref not in references:
            raise AlignmentError(
                f"--ref {ref!r} not found in {os.path.basename(path)}; "
                f"header references: {', '.join(references) or '(none)'}"
            )
        return ref
    if seqid is not None and (seqid in references or not references):
        return seqid
    if len(references) == 1:
        (only,) = references
        if seqid is not None and only != seqid:
            log.warning("GFF sequence id %r not in %s; using its only reference %r",
                        seqid, os.path.basename(path), only)
        return only
    raise AlignmentError(
        f"Cannot decide which reference to use from {os.path.basename(path)} "
        f"(GFF id {seqid!r}; header references: {', '.join(references) or '(none)'}). "
        "Pass --ref NAME."
    )


# --------------------------------------------------------------------------
# SAM (text)
# --------------------------------------------------------------------------

def _parse_cigar(text: str) -> tuple[tuple[int, str], ...]:
    if text == "*":
        return ()
    parts = _CIGAR_RE.findall(text)
    if sum(len(n) + 1 for n, _ in parts) != len(text):
        raise AlignmentError(f"Malformed CIGAR string {text!r}")
    return tuple((int(n), op) for n, op in parts)


def _open_sam(path: str) -> AlignmentFile:
    fp = open(path, "rt", encoding="utf-8", errors="replace")
    references: dict[str, int] = {}
    first_record: str | None = None

    for line in fp:
        if line.startswith("@"):
            if line.startswith("@SQ"):
                fields = dict(f.split(":", 1) for f in line.rstrip("\n").split("\t")[1:] if ":" in f)
                if "SN" in fields:
                    references[fields["SN"]] = int(fields.get("LN", 0))
            continue
        first_record = line
        break

    def records() -> Iterator[Alignment]:
        with fp:
            lines = [first_record] if first_record is not None else []
            for lineno, line in enumerate(_chain(lines, fp), 1):
                if not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 11:
                    raise AlignmentError(
                        f"Malformed SAM record in {os.path.basename(path)}:{lineno} — "
                        f"expected ≥11 tab-separated columns, got {len(f)}"
                    )
                try:
                    flag, pos, mapq = int(f[1]), int(f[3]), int(f[4])
                except ValueError:
                    raise AlignmentError(
                        f"Non-integer FLAG/POS/MAPQ in {os.path.basename(path)}:{lineno}"
                    ) from None
                yield Alignment(ref=f[2], pos=pos, flag=flag, mapq=mapq,
                                cigar=_parse_cigar(f[5]))

    return AlignmentFile(path=path, references=references, records=records())


def _chain(first, rest):
    yield from first
    yield from rest


# --------------------------------------------------------------------------
# BAM (binary, BGZF)
# --------------------------------------------------------------------------

def _open_bam(path: str) -> AlignmentFile:
    fp = gzip.open(path, "rb")
    if fp.read(4) != b"BAM\x01":
        fp.close()
        raise AlignmentError(f"{os.path.basename(path)} is gzip-compressed but not BAM")

    (l_text,) = struct.unpack("<i", fp.read(4))
    fp.read(l_text)                                   # header text (ignored)
    (n_ref,) = struct.unpack("<i", fp.read(4))
    names: list[str] = []
    references: dict[str, int] = {}
    for _ in range(n_ref):
        (l_name,) = struct.unpack("<i", fp.read(4))
        name = fp.read(l_name)[:-1].decode()          # strip NUL
        (l_ref,) = struct.unpack("<i", fp.read(4))
        names.append(name)
        references[name] = l_ref

    def records() -> Iterator[Alignment]:
        with fp:
            while True:
                head = fp.read(4)
                if len(head) < 4:
                    return
                (block_size,) = struct.unpack("<i", head)
                block = fp.read(block_size)
                if len(block) < block_size:
                    raise AlignmentError(f"Truncated BAM record in {os.path.basename(path)}")
                (ref_id, pos, l_read_name, mapq, _bin, n_cigar, flag, _l_seq,
                 _next_ref, _next_pos, _tlen) = struct.unpack_from("<iiBBHHHiiii", block, 0)
                off = 32 + l_read_name
                raw = struct.unpack_from(f"<{n_cigar}I", block, off)
                cigar = tuple((c >> 4, _BAM_CIGAR_OPS[c & 0xF]) for c in raw)
                ref = names[ref_id] if 0 <= ref_id < len(names) else "*"
                yield Alignment(ref=ref, pos=pos + 1, flag=flag, mapq=mapq, cigar=cigar)

    return AlignmentFile(path=path, references=references, records=records())
