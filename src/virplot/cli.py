"""CLI entry point for VirPlot."""

from __future__ import annotations

import argparse
import logging
import os
import sys

from virplot import __version__
from virplot.analysis import call_blocks, write_csvs
from virplot.models import RNA
from virplot.alignments import AlignmentError
from virplot.parsers import default_label, load_depth, parse_gff_rnas
from virplot.plotting import CircularPlotter, LinearPlotter
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
                   help="One or more depth sources: samtools depth output, or "
                        "SAM/BAM alignments (stacked if multiple)")
    p.add_argument("-l", "--labels", nargs="+",
                   help="Label(s) for depth line (same order as --depth)")
    p.add_argument("--ref",
                   help="Reference name to use from SAM/BAM/depth files when the GFF "
                        "sequence id does not match (single-RNA GFF only)")
    p.add_argument("--rnas", nargs="+", metavar="SEQID",
                   help="Plot only these GFF sequence ids, in this order "
                        "(default: every 'region' in the GFF)")
    p.add_argument("--topology", choices=["auto", "circular", "linear"], default="auto",
                   help="Treat the genome(s) as circular or linear; 'auto' follows the "
                        "GFF region line's Is_circular attribute [%(default)s]")
    p.add_argument("--layout", choices=["linear", "circular", "auto"], default="linear",
                   help="Figure layout: a linear track, a circular (polar) plot, or "
                        "'auto' to draw circular genomes as circles [%(default)s]")
    p.add_argument("--free-y", action="store_true",
                   help="With several RNAs, let each figure pick its own depth y-limit "
                        "instead of sharing one")
    p.add_argument("--equal-width", action="store_true",
                   help="With several RNAs, draw every figure at full width instead of "
                        "proportional to RNA length")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="Skip SAM/BAM reads with MAPQ below this [%(default)s]")
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
    p.add_argument("--border", action="store_true",
                   help="Outline feature glyphs in black (default: no outline)")
    p.add_argument("--no-border", action="store_true",
                   help=argparse.SUPPRESS)          # former default; kept so old commands run
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

    # --- parse annotations ---
    rnas = parse_gff_rnas(args.gff)
    if args.rnas:
        by_id = {r.seqid: r for r in rnas}
        missing = [x for x in args.rnas if x not in by_id]
        if missing:
            log.error("--rnas not found in GFF: %s (available: %s)",
                      ", ".join(missing), ", ".join(by_id))
            sys.exit(1)
        rnas = [by_id[x] for x in args.rnas]
    if args.ref and len(rnas) > 1:
        log.error("--ref applies to a single RNA; use --rnas to select one, "
                  "or drop --ref so each RNA is matched by its GFF sequence id")
        sys.exit(1)
    if args.topology != "auto":
        for rna in rnas:
            rna.circular = args.topology == "circular"

    multi = len(rnas) > 1
    for rna in rnas:
        shape = "circular" if rna.circular else "linear"
        if multi:
            log.info("%s: %d bp, %d features, %s", rna.seqid, rna.length,
                     len(rna.features), shape)
        else:
            log.info("Parsed %d features from GFF (%d bp, %s)",
                     len(rna.features), rna.length, shape)

    # --- depth tracks ---
    labels = args.labels or [default_label(f) for f in args.depth]
    for rna in rnas:
        counts: set[int] = set()
        for label, df in zip(labels, args.depth):
            try:
                y, n, kind = load_depth(df, rna.length, seqid=rna.seqid,
                                        ref=args.ref, min_mapq=args.min_mapq,
                                        circular=rna.circular)
            except AlignmentError as exc:
                log.error("%s", exc)
                sys.exit(1)
            where = f"{rna.seqid} in {df}" if multi else df
            if kind == "alignments":
                log.info("Counted %d aligned reads for %s", n, where)
            else:
                log.info("Parsed %d depth entries for %s", n, where)
                counts.add(n)
            rna.add_depth(label, y)
        if len(counts) > 1:
            log.error("Mismatching position count across depth files for %s: %s",
                      rna.seqid, counts)
            sys.exit(1)

    settings = load_settings(args.yaml)
    log.info("Loaded settings from %s", args.yaml)

    os.makedirs(args.outdir, exist_ok=True)
    # One prepared plotter per layout, so cross-RNA scaling is shared; --layout
    # 'auto' can mix the two when a GFF holds both circular and linear RNAs.
    linear_plotter = LinearPlotter(settings, args)
    linear_plotter.prepare(rnas)
    circular_plotter = CircularPlotter(settings, args)
    circular_plotter.prepare(rnas)

    for rna in rnas:
        out_base = f"{args.name}.{rna.seqid}" if multi else args.name
        prefix = f"{rna.seqid}: " if multi else ""

        # --- thresholds ---
        y_sum = rna.total_depth().astype(int)
        threshold_results: list[tuple] = []
        for T in args.thresholds:
            intervals, gaps, pct = call_blocks(y_sum, T)
            threshold_results.append((T, intervals, gaps, pct))
            log.info("%sT=%d: %d intervals, %d breaks, %.2f%% genome >=%dx",
                     prefix, T, len(intervals), len(gaps), pct, T)

        if args.report:
            for T, intervals, gaps, _ in threshold_results:
                write_csvs(intervals, gaps, args.outdir, out_base, T)

        # --- plot ---
        circular_layout = args.layout == "circular" or (
            args.layout == "auto" and rna.circular)
        (circular_plotter if circular_layout else linear_plotter).render(
            rna, threshold_results, out_base)
