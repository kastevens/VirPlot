"""GFF3 and depth-file parsers, plus the dispatcher that turns any depth
source (samtools depth text, SAM, BAM) into a per-base array."""

from __future__ import annotations

import logging
import os
import sys

import numpy as np

from virplot.alignments import depth_from_alignments, is_alignment_file, resolve_reference
from virplot.models import RNA, Feature

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# GFF3
# --------------------------------------------------------------------------

def parse_gff_rnas(gff_path: str) -> list[RNA]:
    """Parse a GFF3 file into one ``RNA`` per ``region`` line.

    RNAs are returned in the order their ``region`` lines appear. CDS features
    are attached to the RNA whose sequence id (column 1) they carry; a CDS on
    a sequence id with no ``region`` line is skipped with a warning.
    """
    rnas: dict[str, RNA] = {}
    orphans: set[str] = set()

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
                    attrs = dict(
                        kv.split("=", 1) for kv in attributes.split(";") if "=" in kv
                    )
                    circular = attrs.get("Is_circular", "").lower() == "true"
                    rnas[seqid] = RNA(name=seqid, seqid=seqid, length=int(end),
                                      circular=circular)
                continue

            if feature_type != "CDS":
                continue

            info = dict(
                kv.split("=", 1) for kv in attributes.split(";") if "=" in kv
            )
            feat = Feature(
                start=int(start),
                end=int(end),
                strand=strand,
                product=info.get("product", "unknown"),
                gene=info.get("gene") or None,
            )
            if seqid in rnas:
                rnas[seqid].features.append(feat)
            else:
                orphans.add(seqid)

    if not rnas:
        log.error("No 'region' feature found in GFF: %s", gff_path)
        sys.exit(1)

    for seqid in sorted(orphans):
        log.warning("CDS features on %r ignored: no 'region' line for that sequence", seqid)

    return list(rnas.values())


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
                log.error("Malformed depth line in %s:%d — expected 3 "
                          "tab-separated columns, got %d", depth_path, lineno, len(fields))
                sys.exit(1)
            try:
                pos, cov = int(fields[1]), int(fields[2])
            except ValueError:
                log.error("Non-integer value in %s:%d — %r", depth_path, lineno, line.rstrip())
                sys.exit(1)
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
