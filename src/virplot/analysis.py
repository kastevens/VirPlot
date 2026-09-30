"""Coverage interval/gap analysis, CSV reporting, and read mapping with bowtie2."""

from __future__ import annotations

import csv
import logging
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field

import numpy as np

log = logging.getLogger(__name__)


def smooth_depth(y_vals: np.ndarray, window_size: int = 10) -> np.ndarray:
    """Smooth a depth array using a moving average kernel."""
    kernel = np.ones(window_size) / window_size
    return np.convolve(y_vals, kernel, mode="same")


def call_blocks(y: np.ndarray, threshold: int) -> tuple[list[dict], list[dict], float]:
    """Identify contiguous intervals >= threshold and gaps < threshold.

    Returns (intervals, gaps, pct_covered).
    """
    covered = y >= threshold

    edges = np.diff(np.r_[False, covered, False].astype(int))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]

    intervals = []
    for s, e in zip(starts, ends):
        seg = y[s:e]
        intervals.append({
            "start_bp": int(s + 1),
            "end_bp": int(e),
            "length_bp": int(e - s),
            "min_depth": int(seg.min()) if seg.size else 0,
            "mean_depth": float(seg.mean()) if seg.size else 0.0,
        })

    not_covered = ~covered
    edges2 = np.diff(np.r_[False, not_covered, False].astype(int))
    g_starts = np.where(edges2 == 1)[0]
    g_ends = np.where(edges2 == -1)[0]

    gaps = [
        {"start_bp": int(s + 1), "end_bp": int(e), "length_bp": int(e - s)}
        for s, e in zip(g_starts, g_ends)
    ]

    covered_bp = sum(iv["length_bp"] for iv in intervals)
    pct_covered = 100.0 * covered_bp / len(y) if len(y) else 0.0

    return intervals, gaps, pct_covered


def write_csvs(intervals: list[dict], gaps: list[dict],
               outdir: str, base: str, threshold: int) -> None:
    """Write interval and gap CSV reports for a given threshold."""
    iv_path = os.path.join(outdir, f"{base}.intervals_ge{threshold}.csv")
    gp_path = os.path.join(outdir, f"{base}.gaps_lt{threshold}.csv")

    with open(iv_path, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=["start_bp", "end_bp", "length_bp", "min_depth", "mean_depth"])
        w.writeheader()
        w.writerows(intervals)

    with open(gp_path, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=["start_bp", "end_bp", "length_bp"])
        w.writeheader()
        w.writerows(gaps)

    log.info("Reports written: %s, %s", os.path.basename(iv_path), os.path.basename(gp_path))


# ==========================================================================
# Read mapping with bowtie2
# ==========================================================================
#
# VirPlot can start from FASTQ: reads are mapped to the reference with
# bowtie2 and the resulting SAM goes straight into the depth pipeline (no
# samtools needed — VirPlot reads SAM itself). bowtie2 and bowtie2-build must
# be on PATH; nothing here imports them, so the rest of the package works
# without them installed.

BOWTIE2_TOOLS = ("bowtie2", "bowtie2-build")
_FASTQ_SUFFIX = re.compile(r"(\.f(ast)?q)?(\.gz)?$", re.I)
_MATE_SUFFIX = re.compile(r"[._-](R?[12]|read[12])$", re.I)


class MappingError(RuntimeError):
    """A mapping step could not run: tool missing, bad inputs, or bowtie2 failed."""


@dataclass
class ReadSet:
    """One sample's reads, in bowtie2's terms: unpaired files, or mate-1/mate-2 files.

    ``unpaired``/``mate1``/``mate2`` are lists because bowtie2 accepts several
    comma-separated files per option; they are joined with commas again when
    the command is built.
    """

    label: str
    unpaired: list[str] = field(default_factory=list)
    mate1: list[str] = field(default_factory=list)
    mate2: list[str] = field(default_factory=list)

    @property
    def paired(self) -> bool:
        return bool(self.mate1)

    @property
    def files(self) -> list[str]:
        return self.unpaired + self.mate1 + self.mate2

    def validate(self) -> None:
        if self.unpaired and (self.mate1 or self.mate2):
            raise MappingError(f"{self.label}: give either -U or -1/-2, not both")
        if len(self.mate1) != len(self.mate2):
            raise MappingError(
                f"{self.label}: -1 has {len(self.mate1)} file(s) but -2 has {len(self.mate2)}")
        if not self.files:
            raise MappingError(f"{self.label}: no read files")
        for f in self.files:
            if not os.path.isfile(f):
                raise MappingError(f"{self.label}: read file not found: {f}")


def default_read_label(files: list[str]) -> str:
    """Sample name from the first read file: strip .fastq/.fq/.gz and a mate suffix."""
    base = os.path.basename(files[0])
    base = _FASTQ_SUFFIX.sub("", base)
    base = _FASTQ_SUFFIX.sub("", base)          # twice: name.fastq.gz -> name.fastq -> name
    return _MATE_SUFFIX.sub("", base) or base


def missing_bowtie2_tools() -> list[str]:
    """Names of the bowtie2 executables not found on PATH (empty when ready)."""
    return [t for t in BOWTIE2_TOOLS if shutil.which(t) is None]


def fasta_lengths(fasta: str) -> dict[str, int]:
    """Sequence name (first word of the header) -> length, for every record."""
    lengths: dict[str, int] = {}
    name = None
    with open(fasta) as fp:
        for line in fp:
            if line.startswith(">"):
                name = line[1:].split()[0] if line[1:].strip() else ""
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def ensure_bowtie2_index(fasta: str, outdir: str, threads: int = 1,
                         runner=subprocess.run) -> str:
    """Return a bowtie2 index basename for ``fasta``, building it if needed.

    An index already sitting next to the FASTA (``<fasta>.1.bt2`` …) is used
    as is. Otherwise one is built under ``<outdir>/<fasta stem>.bt2/`` and
    reused on later runs while it is newer than the FASTA.
    """
    if not os.path.isfile(fasta):
        raise MappingError(f"Reference FASTA not found: {fasta}")

    def complete(base: str) -> bool:
        parts = [f"{base}.{i}.bt2" for i in (1, 2, 3, 4)] + [f"{base}.rev.{i}.bt2" for i in (1, 2)]
        return all(os.path.isfile(p) for p in parts)

    if complete(fasta):
        log.info("Using existing bowtie2 index beside %s", os.path.basename(fasta))
        return fasta

    stem = _FASTQ_SUFFIX.sub("", os.path.splitext(os.path.basename(fasta))[0])
    index_dir = os.path.join(outdir, f"{stem}.bt2")
    base = os.path.join(index_dir, stem)
    if complete(base) and os.path.getmtime(f"{base}.1.bt2") >= os.path.getmtime(fasta):
        log.info("Using cached bowtie2 index %s", index_dir)
        return base

    os.makedirs(index_dir, exist_ok=True)
    cmd = ["bowtie2-build", "--threads", str(threads), "-q", fasta, base]
    log.info("Building bowtie2 index: %s", " ".join(cmd))
    _run(cmd, runner, "bowtie2-build")
    return base


def map_reads(reads: ReadSet, index_base: str, out_sam: str, threads: int = 1,
              extra_args: list[str] | None = None, runner=subprocess.run) -> str:
    """Map one sample with bowtie2 and return the SAM path it wrote."""
    reads.validate()
    cmd = ["bowtie2", "-p", str(threads), "-x", index_base]
    if reads.paired:
        cmd += ["-1", ",".join(reads.mate1), "-2", ",".join(reads.mate2)]
    else:
        cmd += ["-U", ",".join(reads.unpaired)]
    cmd += list(extra_args or []) + ["-S", out_sam]
    log.info("Mapping %s: %s", reads.label, " ".join(cmd))
    result = _run(cmd, runner, "bowtie2")
    # bowtie2 prints its alignment summary on stderr; keep it, it is the QC
    for line in (result.stderr or "").strip().splitlines():
        log.info("  bowtie2 | %s", line.rstrip())
    return out_sam


def _run(cmd: list[str], runner, tool: str):
    if shutil.which(cmd[0]) is None:
        raise MappingError(f"{tool} not found on PATH (install bowtie2, e.g. conda install -c bioconda bowtie2)")
    try:
        result = runner(cmd, check=False, capture_output=True, text=True)
    except OSError as exc:
        raise MappingError(f"could not start {tool}: {exc}") from exc
    if result.returncode != 0:
        tail = (result.stderr or "").strip().splitlines()[-5:]
        raise MappingError(f"{tool} exited with status {result.returncode}:\n  " + "\n  ".join(tail))
    return result
