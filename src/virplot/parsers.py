"""GFF3 and depth file parsers."""

from __future__ import annotations

import logging
import sys

import numpy as np

from virplot.models import Feature

log = logging.getLogger(__name__)


def parse_gff(gff_path: str) -> tuple[int, list[Feature]]:
    """Parse a GFF3 file, returning (sequence_length, features)."""
    features: list[Feature] = []
    sequence_length: int | None = None

    with open(gff_path) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            _, _, feature_type, start, end, _, strand, _, attributes = parts

            if feature_type == "region" and sequence_length is None:
                sequence_length = int(end)
                continue

            if feature_type != "CDS":
                continue

            info = dict(
                kv.split("=", 1) for kv in attributes.split(";") if "=" in kv
            )

            features.append(Feature(
                start=int(start),
                end=int(end),
                strand=strand,
                product=info.get("product", "unknown"),
            ))

    if sequence_length is None:
        log.error("No 'region' feature found in GFF: %s", gff_path)
        sys.exit(1)

    return sequence_length, features


def parse_depth(depth_path: str, seq_len: int) -> tuple[np.ndarray, int]:
    """Parse a samtools depth file directly into a numpy array.

    Returns (depth_array, n_entries).
    """
    y = np.zeros(seq_len, dtype=int)
    n = 0
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
            if 1 <= pos <= seq_len:
                y[pos - 1] = cov
            n += 1
    return y, n
