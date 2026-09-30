"""CLI entry point for VirPlot."""

from __future__ import annotations

import argparse
import logging
import os
import sys

from virplot import __version__
from virplot.analysis import call_blocks, write_csvs
from virplot.models import RNA
from virplot.parsers import parse_gff, parse_depth
from virplot.plotting import LinearPlotter
from virplot.settings import load_settings

log = logging.getLogger(__name__)


class _Formatter(logging.Formatter):
    """Show log level only for non-INFO messages."""

    def format(self, record: logging.LogRecord) -> str:
        if record.levelno == logging.INFO:
            self._style._fmt = "[VirPlot] %(message)s"
        else:
            self._style._fmt = "[VirPlot] %(levelname)s: %(message)s"
        return super().format(record)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Virus genome feature annotation and depth plotting tool",
    )
    p.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("-g", "--gff", required=True,
                   help="Path to GFF3 annotation file")
    p.add_argument("-d", "--depth", nargs="+", required=True,
                   help="One or more depth files (stacked if multiple)")
    p.add_argument("-l", "--labels", nargs="+",
                   help="Label(s) for depth line (same order as --depth)")
    p.add_argument("-y", "--yaml", required=True,
                   help="YAML file for color mapping and other specs")
    p.add_argument("-o", "--outdir", default=".",
                   help="Output directory for the plot")
    p.add_argument("-n", "--normalize", action="store_true",
                   help="Normalize depth values to max=1")
    p.add_argument("--grid", action="store_true",
                   help="Enable background grid on depth plot")
    p.add_argument("--smooth", action="store_true",
                   help="Smooth depth plot using moving average")
    p.add_argument("--yscale", choices=["linear", "symlog"], default="linear",
                   help="Y-axis scale method for depth plot [%(default)s]")
    p.add_argument("--linthresh", type=float, default=10.0,
                   help="Symlog linear threshold around 0 [%(default)s]")
    p.add_argument("--name", default="virplot",
                   help="Base name for output file [%(default)s]")
    p.add_argument("--no-label", action="store_true",
                   help="Do not label feature names in feature rectangles")
    p.add_argument("--no-border", action="store_true",
                   help="Do not draw borders around feature rectangles")
    p.add_argument("-t", "--thresholds", nargs="+", type=int, default=[1, 5],
                   help="Coverage thresholds to call intervals and breaks [%(default)s]")
    p.add_argument("-r", "--report", action="store_true",
                   help="Write CSV reports of intervals and breaks for each threshold")
    p.add_argument("--shade-breaks", action="store_true",
                   help="Shade coverage gaps (<T) on the depth plot")
    p.add_argument("--legend", action="store_true",
                   help="Show depth plot legend")
    p.add_argument("--title", action="store_true",
                   help="Show title specified in YAML")
    p.add_argument("-f", "--format", choices=["svg", "pdf", "png"], default="svg",
                   help="Output format [%(default)s]")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Enable debug-level logging")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    handler = logging.StreamHandler()
    handler.setFormatter(_Formatter())
    logging.basicConfig(
        handlers=[handler],
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    # --- validate inputs ---
    for path, label in [(args.gff, "GFF"), (args.yaml, "YAML")]:
        if not os.path.isfile(path):
            log.error("%s file not found: %s", label, path)
            sys.exit(1)
    for f in args.depth:
        if not os.path.isfile(f):
            log.error("Depth file not found: %s", f)
            sys.exit(1)
    if args.labels and len(args.labels) != len(args.depth):
        log.error("Mismatching number of labels [%d] to depth files [%d]",
                  len(args.labels), len(args.depth))
        sys.exit(1)

    # --- parse ---
    sequence_length, features = parse_gff(args.gff)
    log.info("Parsed %d features from GFF", len(features))
    rna = RNA(name=args.name, length=sequence_length, features=features)

    labels = args.labels or [
        os.path.splitext(os.path.basename(f))[0] for f in args.depth
    ]
    counts: set[int] = set()
    for label, df in zip(labels, args.depth):
        y, n = parse_depth(df, rna.length)
        log.info("Parsed %d depth entries from %s", n, df)
        rna.add_depth(label, y)
        counts.add(n)

    if len(counts) != 1:
        log.error("Mismatching position count across depth files: %s", counts)
        sys.exit(1)

    settings = load_settings(args.yaml)
    log.info("Loaded settings from %s", args.yaml)

    # --- thresholds ---
    y_sum = rna.total_depth().astype(int)
    threshold_results: list[tuple] = []
    for T in args.thresholds:
        intervals, gaps, pct = call_blocks(y_sum, T)
        threshold_results.append((T, intervals, gaps, pct))
        log.info("T=%d: %d intervals, %d breaks, %.2f%% genome >=%dx",
                 T, len(intervals), len(gaps), pct, T)

    if args.report:
        for T, intervals, gaps, _ in threshold_results:
            write_csvs(intervals, gaps, args.outdir, args.name, T)

    # --- plot ---
    LinearPlotter(settings, args).render(rna, threshold_results)
