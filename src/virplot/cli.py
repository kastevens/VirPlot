"""CLI entry point for VirPlot."""

from __future__ import annotations

import argparse
import logging
import os
import sys

from virplot import __version__
from virplot.analysis import (MappingError, ReadSet, call_blocks, default_read_label,
                              ensure_bowtie2_index, fasta_lengths, map_reads,
                              missing_bowtie2_tools, write_csvs)
from virplot.models import RNA
from virplot.alignments import AlignmentError
from virplot.parsers import ParseError, default_label, load_depth, parse_gff_rnas
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
    p.add_argument("-d", "--depth", nargs="+", default=[],
                   help="Depth sources: samtools depth output, or SAM/BAM alignments "
                        "(stacked if multiple). Can be combined with reads below")
    p.add_argument("-l", "--labels", nargs="+",
                   help="Label(s) for the depth tracks: -d sources first, then each "
                        "-U sample, then each -1/-2 pair, in command-line order")

    reads = p.add_argument_group(
        "map reads with bowtie2 (options mirror bowtie2)",
        "Give FASTQ instead of, or as well as, -d. Each -U is one unpaired sample and "
        "each -1/-2 pair one paired sample; comma-separate files that belong to the "
        "same sample, as bowtie2 does. Needs bowtie2 and bowtie2-build on PATH.")
    reads.add_argument("-x", "--reference", metavar="FASTA",
                       help="Reference genome FASTA; a bowtie2 index beside it is used, "
                            "otherwise one is built under --outdir and reused")
    reads.add_argument("-U", dest="unpaired", action="append", metavar="FASTQ[,FASTQ...]",
                       help="Unpaired reads for one sample (repeat for more samples)")
    reads.add_argument("-1", dest="mate1", action="append", metavar="FASTQ[,FASTQ...]",
                       help="Mate-1 reads for one paired sample (repeat for more samples)")
    reads.add_argument("-2", dest="mate2", action="append", metavar="FASTQ[,FASTQ...]",
                       help="Mate-2 reads, one -2 per -1")
    reads.add_argument("-p", "--threads", type=int, default=1,
                       help="Threads for bowtie2 and bowtie2-build [%(default)s]")
    reads.add_argument("--bowtie2-args", default="",
                       help="Extra options passed to bowtie2 verbatim, e.g. \"--local --very-sensitive\"")
    p.add_argument("--ref",
                   help="Reference name to use from SAM/BAM/depth files when the GFF "
                        "sequence id does not match (single-RNA GFF only)")
    p.add_argument("--rna-order", choices=["length", "gff"], default="length",
                   help="Order segments largest first (the ICTV stacking order) "
                        "or as the GFF lists them [%(default)s]; --rnas sets its own order")
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
    p.add_argument("--bare-x", action="store_true",
                   help="Drop the x-axis label and tick numbers, for a panel that will be "
                        "stacked under another (see bin/stack_figures.py)")
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
    p.add_argument("--border", action=argparse.BooleanOptionalAction, default=False,
                   help="Outline feature glyphs in black (--border), or draw them as "
                        "flat colour as the ICTV figures do (--no-border)")
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


def _read_sets(args) -> list[ReadSet]:
    """Turn -U / -1 / -2 occurrences into one ReadSet per sample."""
    sets: list[ReadSet] = []
    for spec in args.unpaired or []:
        files = spec.split(",")
        sets.append(ReadSet(label=default_read_label(files), unpaired=files))
    m1, m2 = args.mate1 or [], args.mate2 or []
    if len(m1) != len(m2):
        log.error("-1 given %d time(s) but -2 %d time(s); each paired sample needs one of each",
                  len(m1), len(m2))
        sys.exit(1)
    for a, b in zip(m1, m2):
        fa, fb = a.split(","), b.split(",")
        sets.append(ReadSet(label=default_read_label(fa), mate1=fa, mate2=fb))
    if sets and not args.reference:
        log.error("Reads were given (-U/-1/-2) but no reference: add -x FASTA")
        sys.exit(1)
    return sets


def _map_read_sets(args, read_sets: list[ReadSet], rnas) -> list[str]:
    """Run bowtie2 for every sample; return the SAM paths in sample order."""
    missing = missing_bowtie2_tools()
    if missing:
        raise MappingError(f"{', '.join(missing)} not found on PATH; install bowtie2 "
                           "(e.g. conda install -c bioconda bowtie2)")
    for r in read_sets:
        r.validate()

    # the FASTA is the truth about length; warn early if the GFF disagrees
    lengths = fasta_lengths(args.reference)
    for rna in rnas:
        if rna.seqid in lengths and lengths[rna.seqid] != rna.length:
            log.warning("%s is %d bp in the GFF but %d bp in %s",
                        rna.seqid, rna.length, lengths[rna.seqid],
                        os.path.basename(args.reference))
    absent = [r.seqid for r in rnas if r.seqid not in lengths]
    if absent and len(lengths) > 1:
        log.warning("Not in the reference FASTA: %s (headers there: %s)",
                    ", ".join(absent), ", ".join(list(lengths)[:6]))

    index = ensure_bowtie2_index(args.reference, args.outdir, threads=args.threads)
    extra = args.bowtie2_args.split() if args.bowtie2_args else []
    sams = []
    for r in read_sets:
        out_sam = os.path.join(args.outdir, f"{args.name}.{r.label}.sam")
        sams.append(map_reads(r, index, out_sam, threads=args.threads, extra_args=extra))
    return sams


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

    read_sets = _read_sets(args)
    if not args.depth and not read_sets:
        log.error("Give at least one depth source (-d) or reads to map (-U or -1/-2)")
        sys.exit(1)
    n_tracks = len(args.depth) + len(read_sets)
    if args.labels and len(args.labels) != n_tracks:
        log.error("Mismatching number of labels [%d] to depth tracks [%d]",
                  len(args.labels), n_tracks)
        sys.exit(1)

    # --- parse annotations ---
    try:
        rnas = parse_gff_rnas(args.gff)
    except ParseError as exc:
        log.error("%s", exc)
        sys.exit(1)
    if args.rnas:
        by_id = {r.seqid: r for r in rnas}
        missing = [x for x in args.rnas if x not in by_id]
        if missing:
            log.error("--rnas not found in GFF: %s (available: %s)",
                      ", ".join(missing), ", ".join(by_id))
            sys.exit(1)
        rnas = [by_id[x] for x in args.rnas]
    elif args.rna_order == "length":
        # ICTV: segments are stacked largest first. An explicit --rnas is the
        # user's own order and is left alone.
        rnas = sorted(rnas, key=lambda r: (-r.length, r.seqid or ""))
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

    # --- labels: -d sources first, then samples; -l overrides both and names the SAMs ---
    labels = args.labels or ([default_label(f) for f in args.depth]
                             + [r.label for r in read_sets])
    for r, label in zip(read_sets, labels[len(args.depth):]):
        r.label = label

    # --- map reads, if any, so they join the depth sources as SAM files ---
    depth_sources = list(args.depth)
    if read_sets:
        os.makedirs(args.outdir, exist_ok=True)
        try:
            depth_sources += _map_read_sets(args, read_sets, rnas)
        except MappingError as exc:
            log.error("%s", exc)
            sys.exit(1)

    # --- depth tracks ---
    for rna in rnas:
        counts: set[int] = set()
        for label, df in zip(labels, depth_sources):
            try:
                y, n, kind = load_depth(df, rna.length, seqid=rna.seqid,
                                        ref=args.ref, min_mapq=args.min_mapq,
                                        circular=rna.circular)
            except (AlignmentError, ParseError) as exc:
                log.error("%s", exc)
                sys.exit(1)
            where = f"{rna.seqid} in {df}" if multi else df
            if kind == "alignments":
                log.info("Counted %d aligned reads for %s", n, where)
            else:
                log.info("Parsed %d depth entries for %s", n, where)
                counts.add(n)
            rna.add_depth(label, y)
        # a circular genome's depth files may legitimately differ in length
        # (one made against a padded reference, one not)
        if len(counts) > 1 and not rna.circular:
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
