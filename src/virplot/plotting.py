"""Matplotlib plotting for genome annotations and depth.

``Plotter`` holds the state every renderer needs (settings, CLI options) and
the helpers that do not depend on figure geometry: preparing depth tracks,
legends, saving.  ``LinearPlotter`` draws the classic two-panel figure — an
annotation track above a depth track sharing one linear x-axis.  A circular
renderer would subclass ``Plotter`` and override ``render`` with polar axes.
"""

from __future__ import annotations

import argparse
import datetime
import logging
import os

import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Ellipse, Polygon, Rectangle
from matplotlib.textpath import TextPath
import matplotlib.pyplot as plt

from virplot.analysis import smooth_depth
from virplot.layout import ABOVE, BELOW, Placement, place_features, resolve_mode
from virplot.models import RNA, Noncoding
from virplot.settings import DEFAULT_FUNCTION_PALETTE, Settings

log = logging.getLogger(__name__)

# layout constants

ANNOTATION_Y_BASE = 0.5
ANNOTATION_HEIGHT = 0.6
ARROW_HEAD_FRACTION = 0.012     # of genome length; arrowhead length cap
NEAR_GAP_FRACTION = 0.01        # same-colour neighbours closer than this flip apart
ORF_LABEL_FONTSIZE = 7          # ORF name written outside the glyph
LABEL_PAD_FRACTION = 0.02
SMALL_FEATURE_THRESHOLD = 500  # bp; features shorter than this get external labels
# Polyprotein domains: the thin line between two mature proteins; the air a
# domain name needs around it to count as fitting inside its segment; and the
# vertical step between rows of outside labels that would otherwise overlap.
DOMAIN_DIVIDER = dict(color="black", linewidth=0.6, solid_capstyle="butt")
TEXT_FIT_MARGIN = 1.2
OUTSIDE_LABEL_ROW = 0.28        # axis units; ~ one 8 pt line in the annotation panel
# Non-coding landmarks sit on the genome line itself: a region (UTR, IR) is a
# bar this tall centred on the line; a stem-loop is a hairpin icon this tall.
NONCODING_BAR_HEIGHT = 0.22
STEM_LOOP_HEIGHT = 0.42
# Subgenomic RNAs: rows beneath everything else, longest first, each a line
# from its 5' end with an arrowhead there. The gap clears the lowest ORF tier.
SGRNA_ROW_STEP = 0.30           # axis units between consecutive sgRNA rows
SGRNA_TOP_GAP = 0.34            # from the lowest ORF tier down to the first row
SGRNA_LINE = dict(linewidth=1.4, solid_capstyle="butt")
SGRNA_LABEL_FONTSIZE = 7
SGRNA_MARKER_SIZE = 4
CIRC_NONCODING_HEIGHT = 0.07    # radial thickness of an IR arc astride the circle
SMOOTH_WINDOW = 15
# Depth trace, measured from the original README figure: a ~0.3 pt black line
# at alpha ~0.8 over fills at alpha 0.9 (the stacked layers have always used 0.9;
# a single sample now matches one layer instead of a paler wash).
DEPTH_OUTLINE = dict(color="black", linewidth=0.3, alpha=0.8)
DEPTH_FILL_ALPHA = 0.9
Y_HEADROOM = 1.05
PNG_DPI = 400
FIGURE_WIDTH = 12
FIGURE_HEIGHT = 4
HEIGHT_RATIO_ANNOTATION = 1
HEIGHT_RATIO_DEPTH = 1

# circular layout: radii as a fraction of the plotted radius, inside out
CIRC_FIGSIZE = 8
CIRC_R_DEPTH_BASE = 0.30        # depth baseline (zero coverage)
CIRC_R_DEPTH_MAX = 0.66         # depth at the y-limit
CIRC_R_BASELINE = 0.71          # the genome circle; ORF arcs sit just outside it
CIRC_LANE_STEP = 0.095          # each nesting tier steps this far inward
CIRC_MAX_TIERS = 3              # deeper nesting than this collides with the depth band
CIRC_ANN_HEIGHT = 0.085
CIRC_R_OUTSIDE_LABEL = 1.01     # labels for features too narrow to hold text
CIRC_RMAX = 1.22                # leaves room for position ticks; rim is hidden
CIRC_SMALL_FEATURE_FRACTION = 0.055  # of the circle; below this the label goes outside


class Plotter:
    """Base renderer: shared state and geometry-independent helpers.

    When several RNAs are rendered from one run, call ``prepare(rnas)`` first:
    it fixes a common normalisation denominator, y-limit and reference length
    so the separate figures are directly comparable (same y-scale, x-axis
    width proportional to length). For a single RNA it changes nothing.
    """

    def __init__(self, settings: Settings, args: argparse.Namespace):
        self.settings = settings
        self.args = args
        self.shared_denom: float | None = None    # --normalize divisor across RNAs
        self.shared_ymax: float | None = None     # common y-limit across RNAs
        self.max_length: int | None = None        # longest RNA, for width scaling
        self.bottom_seqid: str | None = None      # last segment rendered; keeps the x-axis

    def _bare_x(self, rna: RNA | None) -> bool:
        """True when this panel should omit the x-axis label and tick numbers.

        ``--bare-x`` is for panels that will be stacked by
        bin/stack_figures.py, where one axis row at the bottom serves them
        all. In a run over several segments the last one rendered is that
        bottom row and keeps its axis; a lone panel is bare, since it is only
        asked for when something else carries the axis. The ticks themselves
        stay so the plot area keeps its geometry and the panels still align.
        """
        if not getattr(self.args, "bare_x", False):
            return False
        return (rna is None or self.bottom_seqid is None
                or rna.seqid != self.bottom_seqid)

    def prepare(self, rnas: list[RNA]) -> None:
        """Compute cross-RNA scaling so separate figures share axes conventions."""
        if len(rnas) < 2:
            return
        # the last segment rendered is the bottom row once stacked, so it is
        # the one that keeps the shared x-axis under --bare-x
        self.bottom_seqid = rnas[-1].seqid
        raw_max = max(self._smoothed_total(r).max() for r in rnas) or 1.0
        if self.args.normalize and not self.args.free_y:
            self.shared_denom = raw_max
        if not self.args.free_y:
            self.shared_ymax = 1.0 if self.args.normalize else raw_max
        if not self.args.equal_width:
            self.max_length = max(r.length for r in rnas)

    def render(self, rna: RNA, threshold_results: list[tuple], out_base: str) -> None:
        """Build the figure for ``rna`` and save it as ``<out_base>.<format>``."""
        raise NotImplementedError

    # --- shared helpers -----------------------------------------------------

    def _smoothed_total(self, rna: RNA) -> np.ndarray:
        return self._prepared_tracks(rna, normalize=False)[1]

    def _prepared_tracks(self, rna: RNA, normalize: bool | None = None
                         ) -> tuple[list[np.ndarray], np.ndarray]:
        """Apply --smooth / --normalize to the depth tracks.

        Returns (per-track arrays, combined array). Smoothing is applied per
        track; normalisation divides everything by the combined maximum (of
        this RNA, or of all RNAs after ``prepare``) so the total peaks at 1.
        An RNA with no depth track gets one flat zero track.
        """
        tracks = rna.tracks or [np.zeros(rna.length, dtype=float)]
        if self.args.smooth:
            tracks = [smooth_depth(y, window_size=SMOOTH_WINDOW) for y in tracks]
        total = np.sum(tracks, axis=0) if len(tracks) > 1 else tracks[0]
        if self.args.normalize if normalize is None else normalize:
            denom = self.shared_denom or total.max() or 1.0
            tracks = [y / denom for y in tracks]
            total = total / denom
        return tracks, total

    def _figure_width(self, rna: RNA) -> float:
        if self.max_length:
            return FIGURE_WIDTH * rna.length / self.max_length
        return FIGURE_WIDTH

    def _glyph_edge(self) -> dict:
        """Outline for feature glyphs: none by default (flat colour, as in the
        ICTV figures, ``--no-border``); a thin black edge with ``--border``."""
        if getattr(self.args, "border", False):
            return dict(edgecolor="black", linewidth=1.0)
        return dict(edgecolor="none")

    def _feature_color(self, product: str) -> str:
        return self.settings.feature_color(product)

    @staticmethod
    def _text_on(fill: str, dark_below: float = 0.35) -> str:
        """Black text on a light fill, white on a dark one (relative luminance).

        The default threshold keeps black on every ICTV palette colour (the
        figures write black even on purple) and flips only fills a user has
        mapped really dark; the non-coding grey passes a higher one.
        """
        r, g, b = to_rgb(fill)
        return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < dark_below else "black"

    def _noncoding_color(self, nc: Noncoding) -> str:
        """Grey for every non-coding landmark (ICTV), unless its label is mapped."""
        s = self.settings
        if nc.label in s.color_mapping:
            return s.color_mapping[nc.label]
        return s.function_palette.get("noncoding", DEFAULT_FUNCTION_PALETTE["noncoding"])

    @staticmethod
    def _mid(rna: RNA, nc) -> float:
        spans = rna.feature_spans(nc)
        covered = sum(e - s + 1 for s, e in spans)
        return ((spans[0][0] - 1 + covered / 2) % rna.length) + 1

    @staticmethod
    def _noncoding_labels(rna: RNA, regions: list, loops: list):
        """Decide where each non-coding name is written, once.

        A region that holds a stem-loop (the begomovirus common region with
        its hairpin) is named at the hairpin's tip, as the ICTV figures write
        ``CRA`` over the icon, unless the stem-loop has a name of its own —
        then the region keeps its name on the arc. Returns
        ``({id(loop): text}, [(region, text), ...])``.
        """
        loop_labels = {id(nc): nc.label for nc in loops}
        region_labels = []
        for reg in regions:
            if not reg.label:
                continue
            inside = [lp for lp in loops if not lp.label and any(
                s <= Plotter._mid(rna, lp) <= e for s, e in rna.feature_spans(reg))]
            if inside:
                loop_labels[id(inside[0])] = reg.label
            else:
                region_labels.append((reg, reg.label))
        return loop_labels, region_labels

    @staticmethod
    def _text_width_pt(text: str, fontsize: float) -> float:
        """Width of ``text`` at ``fontsize`` in points, from the font's own metrics.

        ``TextPath`` lays the string out with the default font without needing
        a renderer, so the answer is the same for every backend.
        """
        return TextPath((0, 0), text, size=fontsize).get_extents().width

    def _text_fits(self, text: str, room_pt: float, fontsize: float) -> bool:
        """Whether ``text`` fits in ``room_pt`` points with a little air either side."""
        return self._text_width_pt(text, fontsize) * TEXT_FIT_MARGIN <= room_pt


    def _placements(self, rna: RNA, circular: bool) -> tuple[list[Placement], str]:
        """Lay the features out; returns (placements, rule used).

        In flip mode a glyph also flips when its upstream neighbour is the same
        colour and closer than ``NEAR_GAP_FRACTION`` of the genome: without an
        outline two such boxes would merge into one (ICTV draws BYV's CP
        below CPm for exactly this reason).
        """
        mode = "nest" if circular else resolve_mode(self.settings.overlap_mode, rna.two_strand)
        near = NEAR_GAP_FRACTION * rna.length

        def crowded(prev, feat) -> bool:
            gap = feat.start - prev.end
            return 0 <= gap < near and self._feature_color(prev.product) == self._feature_color(feat.product)

        return place_features(rna.features, rna.feature_spans, mode, flip_if=crowded), mode

    @staticmethod
    def _head_length(span_len: int, seq_len: int) -> float:
        return min(span_len * 0.35, seq_len * ARROW_HEAD_FRACTION)

    def _add_legend(self, ax: plt.Axes, layers: list, labels: list[str]) -> None:
        legend = ax.legend(
            handles=reversed(layers),
            labels=labels,
            loc=self.settings.legend_location,
            fontsize=8,
            frameon=True,
            fancybox=False,
            framealpha=1.0,
            facecolor="white",
            edgecolor="black",
        )
        legend.get_frame().set_linewidth(0.5)

    def _add_title(self, fig: plt.Figure, rna: RNA | None = None):
        """Draw the title if --title was given; returns the artist or None.

        With no title in the YAML, the ICTV form ``name (length nts)`` is used.
        """
        if not self.args.title:
            return None
        title = self.settings.title
        if not title and rna is not None:
            title = f"{rna.name} ({rna.length:,} nts)"
        return fig.text(
            0.5, 0.95, title,
            ha="center", va="bottom", fontsize=14, fontweight="bold",
        )

    def _save(self, fig: plt.Figure, extra_artists: list, out_base: str) -> None:
        """Save the figure as ``<outdir>/<out_base>.<format>``."""
        args = self.args
        os.makedirs(args.outdir, exist_ok=True)

        extra = [a for a in extra_artists if a is not None]
        save_kwargs = dict(bbox_inches="tight", bbox_extra_artists=extra, pad_inches=0.05)

        ext = args.format
        if ext == "png":
            save_kwargs["dpi"] = PNG_DPI
        else:
            # reproducible vector output: no timestamp, and element ids hashed
            # from the file name rather than a random salt, so regenerating an
            # unchanged figure yields a byte-identical file
            save_kwargs["metadata"] = {"Date": None} if ext == "svg" else {"CreationDate": None}
            plt.rcParams["svg.hashsalt"] = out_base

        output_path = os.path.join(args.outdir, f"{out_base}.{ext}")

        if os.path.exists(output_path):
            ts = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
            alt = f"{out_base}-{ts}.{ext}"
            output_path = os.path.join(args.outdir, alt)
            log.warning("File already exists. Saving to file: %s", alt)
        else:
            log.info("Saving to file: %s", os.path.basename(output_path))

        fig.savefig(output_path, format=ext, **save_kwargs)
        log.info("Plot saved to: %s", output_path)


class LinearPlotter(Plotter):
    """Annotation track stacked above a depth track on a shared linear x-axis."""

    def render(self, rna: RNA, threshold_results: list[tuple], out_base: str) -> None:
        fig, (ax_ann, ax_depth) = plt.subplots(
            2, 1,
            figsize=(self._figure_width(rna), FIGURE_HEIGHT),
            sharex=True,
            gridspec_kw={"height_ratios": [HEIGHT_RATIO_ANNOTATION, HEIGHT_RATIO_DEPTH]},
        )
        fig.subplots_adjust(hspace=0)

        extra_artists = self._draw_annotations(ax_ann, rna)

        tracks, total = self._prepared_tracks(rna)
        layers = self._draw_depth(ax_depth, rna, tracks, total)

        if self.args.legend:
            self._add_legend(ax_depth, layers, rna.labels)

        if self.args.shade_breaks:
            self._shade_gaps(ax_depth, threshold_results)

        self._style_depth_axis(ax_depth, total, rna)

        extra_artists.append(self._add_title(fig, rna))

        self._save(fig, extra_artists, out_base)
        plt.close(fig)

    # --- panels -------------------------------------------------------------

    def _draw_annotations(self, ax: plt.Axes, rna: RNA) -> list:
        """Genome line, end marks and feature glyphs.

        Placement follows docs/ictv_drawing_conventions.md: one-strand genomes
        flip a box across the line only when it overlaps its upstream
        neighbour; two-strand genomes keep + above / − below with arrows and
        tier same-side overlaps outward. Returns artists outside the axes so
        the saved bbox includes them.
        """
        args, settings = self.args, self.settings
        seq_len = rna.length
        pad = int(seq_len * LABEL_PAD_FRACTION)
        y0, H = ANNOTATION_Y_BASE, ANNOTATION_HEIGHT
        ends = self._draw_backbone(ax, rna, pad)

        placements, mode = self._placements(rna, circular=False)
        arrows = mode == "tier"
        tiers_up = max([p.tier for p in placements if p.side == ABOVE], default=0) + 1
        tiers_dn = max([p.tier for p in placements if p.side != ABOVE], default=0) + 1

        # points per nucleotide once xlim is (0, seq_len): decides whether a
        # domain name fits inside its segment
        pt_per_nt = ax.get_position().width * ax.figure.get_figwidth() * 72 / seq_len
        # outside labels already placed on each side, per row: [(x0, x1), ...]
        outside_rows: dict[int, list[list[tuple[float, float]]]] = {ABOVE: [], BELOW: []}

        for pl in placements:
            feat = pl.feature
            yb = y0 + pl.tier * H if pl.side == ABOVE else y0 - (pl.tier + 1) * H
            color = self._feature_color(feat.product)
            spans = rna.feature_spans(feat)
            head_span = spans[-1] if feat.forward else spans[0]
            # outside-label height for this side: clear of every tier on it
            ly_out = (y0 + tiers_up * H + H / 2 if pl.side == ABOVE
                      else y0 - tiers_dn * H - H / 2)
            for (s_, e_) in spans:
                if arrows and (s_, e_) == head_span:
                    glyph = Polygon(self._arrow_points(s_, e_, yb, H, feat.forward, seq_len),
                                    closed=True, facecolor=color, **self._glyph_edge())
                else:
                    glyph = Rectangle((s_, yb), e_ - s_, H, facecolor=color, **self._glyph_edge())
                ax.add_patch(glyph)
                if feat.domains:
                    self._draw_domains(ax, feat, s_, e_, yb, H, glyph, ly_out, pt_per_nt,
                                       labels=not args.no_label,
                                       rows=outside_rows[pl.side], up=pl.side == ABOVE)

            self._draw_mechanism(ax, feat, spans, yb, H, pl.side == ABOVE)

            if args.no_label or not feat.show_label:
                continue
            covered = sum(e_ - s_ + 1 for s_, e_ in spans)
            mid = ((spans[0][0] - 1 + covered / 2) % seq_len) + 1 if rna.circular else feat.midpoint
            small = covered < SMALL_FEATURE_THRESHOLD
            if not feat.domains:                    # a polyprotein's domains label its box
                ax.text(mid, ly_out if small else yb + H / 2, feat.product,
                        ha="center", va="center", fontsize=settings.annotation_fontsize,
                        color="black" if small else self._text_on(color))
            if feat.gene and not small:
                gy = yb + H + 0.12 if pl.side == ABOVE else yb - 0.12
                ax.text(mid, gy, feat.gene, ha="center",
                        va="bottom" if pl.side == ABOVE else "top",
                        fontsize=ORF_LABEL_FONTSIZE, color="0.25")

        ends.extend(self._draw_noncoding(ax, rna, placements, labels=not args.no_label))

        # extra rows of outside domain labels push the panel's edge out
        extra_up = max(len(outside_rows[ABOVE]) - 1, 0) * OUTSIDE_LABEL_ROW
        extra_dn = max(len(outside_rows[BELOW]) - 1, 0) * OUTSIDE_LABEL_ROW

        # sgRNA ladder hangs below the lowest ORF tier and its outside labels
        sg_top = y0 - tiers_dn * H - extra_dn - SGRNA_TOP_GAP
        sg_bottom = self._draw_sgrnas(ax, rna, sg_top, labels=not args.no_label)

        ax.set_xlim(0, seq_len)
        ax.set_ylim(min(-1.5, sg_bottom - 0.5, y0 - tiers_dn * H - 0.9 - extra_dn),
                    max(2.0, y0 + tiers_up * H + 0.9 + extra_up))
        ax.axis("off")
        return ends

    def _draw_sgrnas(self, ax: plt.Axes, rna: RNA, top: float,
                     labels: bool = True) -> float:
        """Subgenomic RNA rows beneath the genome, longest first (ICTV §1).

        Each row is a line from the sgRNA's 5' end to its 3' end with an
        arrowhead at the 5' end, labelled just outside that end. Rows are
        annotation only: a plant virus sgRNA is co-linear with the genome, so
        its reads cannot be told from genomic reads and it has no depth track
        of its own — but each 5' end predicts a step in the shared depth
        trace below, which is the point of drawing them on the same x-axis.

        Returns the y of the lowest row drawn (``top`` when there are none),
        so the caller can extend the panel to fit them.
        """
        ladder = rna.sgrna_ladder
        if not ladder:
            return top

        color = self.settings.sgrna_color
        pad = rna.length * LABEL_PAD_FRACTION
        y = top
        for sg in ladder:
            y -= SGRNA_ROW_STEP
            ax.plot([sg.start, sg.end], [y, y], color=color, **SGRNA_LINE)
            ax.plot([sg.start], [y], marker=">", color=color,
                    markersize=SGRNA_MARKER_SIZE, linestyle="none")
            if labels:
                ax.text(sg.start - pad * 0.35, y, sg.label, ha="right", va="center",
                        fontsize=SGRNA_LABEL_FONTSIZE, color=color)
        return y

    def _draw_domains(self, ax: plt.Axes, feat, s_: int, e_: int, yb: float, H: float,
                      glyph, ly_out: float, pt_per_nt: float, labels: bool,
                      rows: list, up: bool) -> None:
        """Split a polyprotein's box into its mature proteins (ICTV §1).

        Each domain is painted in its own colour over the parent glyph and
        clipped to it, so an arrowhead keeps its shape and any unannotated
        stretch keeps the polyprotein's colour; thin vertical lines mark the
        boundaries. A domain name goes inside its segment when it fits, else
        outside the box at the small-feature label height, as Potyviridae
        Fig. 2 writes ``6K1``/``6K2`` beside the box. Outside labels that
        would overlap one already placed on this side step out one row
        (``rows`` keeps the occupied x-ranges per row across features).
        """
        fontsize = self.settings.annotation_fontsize
        for d in feat.domains_within(s_, e_):
            fill = self._feature_color(d.product)
            seg = Rectangle((d.start, yb), d.end - d.start, H, facecolor=fill, **self._glyph_edge())
            ax.add_patch(seg)
            seg.set_clip_path(glyph)
            if not labels:
                continue
            mid = (d.start + d.end) / 2
            if self._text_fits(d.product, (d.end - d.start) * pt_per_nt, fontsize):
                ax.text(mid, yb + H / 2, d.product, ha="center", va="center",
                        fontsize=fontsize, color=self._text_on(fill))
                continue
            half = self._text_width_pt(d.product, fontsize) * TEXT_FIT_MARGIN / pt_per_nt / 2
            row = self._free_row(rows, mid - half, mid + half)
            ly = ly_out + row * OUTSIDE_LABEL_ROW * (1 if up else -1)
            ax.text(mid, ly, d.product, ha="center", va="center",
                    fontsize=fontsize, color="black")
        for x in feat.dividers_within(s_, e_):
            (line,) = ax.plot([x, x], [yb, yb + H], zorder=3, **DOMAIN_DIVIDER)
            line.set_clip_path(glyph)                 # stays inside an arrowhead

    @staticmethod
    def _free_row(rows: list[list[tuple[float, float]]], x0: float, x1: float) -> int:
        """Index of the first row where ``[x0, x1]`` overlaps nothing; records it there."""
        for i, taken in enumerate(rows):
            if all(x1 < a or x0 > b for a, b in taken):
                taken.append((x0, x1))
                return i
        rows.append([(x0, x1)])
        return len(rows) - 1

    def _draw_mechanism(self, ax: plt.Axes, feat, spans, yb: float, H: float, above: bool) -> None:
        """Mark how a continuation ORF is reached, as the ICTV figures do.

        Frameshift: the box has already flipped across the line; write ``+1 FS``
        / ``−1 FS`` at the junction, outside the box on its far side.
        Readthrough: the extension sits on its partner's side; draw a thin bar
        at the read-through stop and write ``RT`` beside it.

        LIMIT: these marks exist in the linear layout only. ``CircularPlotter``
        places frameshift and readthrough ORFs correctly (the nesting rule does
        not care how an ORF is reached) but draws no ``FS``/``RT`` text; see
        the note in its ``_draw_annotations``.
        """
        if feat.mechanism not in ("frameshift", "readthrough"):
            return
        junction = spans[0][0] if feat.forward else spans[-1][1]
        far_y = yb + H + 0.08 if above else yb - 0.08
        va = "bottom" if above else "top"
        if feat.mechanism == "frameshift":
            # the continuation has flipped, so its partner ends at the junction on
            # the other side; write the label there, outside the partner's box,
            # ending at the junction (CTV: "+1FS" under the 3' end of ORF1a)
            sign = ("" if feat.shift is None else
                    f"{'+' if feat.shift > 0 else chr(0x2212)}{abs(feat.shift)} ")
            partner_y = (ANNOTATION_Y_BASE - 0.08) if above else (ANNOTATION_Y_BASE + H + 0.08)
            ax.text(junction, partner_y, f"{sign}FS", ha="right" if feat.forward else "left",
                    va="top" if above else "bottom", fontsize=ORF_LABEL_FONTSIZE, color="black")
        else:
            ax.plot([junction, junction], [yb - 0.06, yb + H + 0.06], color="black",
                    linewidth=1.0, solid_capstyle="butt", zorder=4)
            ax.text(junction, far_y, "RT", ha="center", va=va,
                    fontsize=ORF_LABEL_FONTSIZE, color="black")

    def _draw_noncoding(self, ax: plt.Axes, rna: RNA, placements, labels: bool) -> list:
        """UTRs, intergenic regions and stem-loops, on the genome line (ICTV §1).

        A region is a grey bar astride the line, under the ORF boxes; its name
        goes in small text on whichever side of the line has no box at that
        position (below first). A stem-loop is a hairpin icon rising from the
        line (dropping below it when a box is above), named at its tip, as the
        ambisense-segment and geminivirus figures draw their hairpins.
        """
        y0, H = ANNOTATION_Y_BASE, ANNOTATION_HEIGHT
        extra: list = []

        def box_at(x: float, above: bool) -> bool:
            return any((pl.side == ABOVE) == above
                       and any(s_ <= x <= e_ for s_, e_ in rna.feature_spans(pl.feature))
                       for pl in placements)

        regions = [nc for nc in rna.noncoding if nc.kind != "stem_loop"]
        loops = [nc for nc in rna.noncoding if nc.kind == "stem_loop"]
        loop_labels, region_labels = self._noncoding_labels(rna, regions, loops)
        region_text = {id(nc): text for nc, text in region_labels}

        for nc in regions:
            color = self._noncoding_color(nc)
            h = NONCODING_BAR_HEIGHT
            for s_, e_ in rna.feature_spans(nc):
                ax.add_patch(Rectangle((s_, y0 - h / 2), e_ - s_, h, facecolor=color,
                                       zorder=1.5, **self._glyph_edge()))
            text = region_text.get(id(nc), "")
            if labels and text:
                mid = self._mid(rna, nc)
                below = not box_at(mid, above=False) or box_at(mid, above=True)
                ax.text(mid, y0 - h / 2 - 0.08 if below else y0 + h / 2 + 0.08, text,
                        ha="center", va="top" if below else "bottom",
                        fontsize=ORF_LABEL_FONTSIZE, color="0.25")

        for nc in loops:
            mid = self._mid(rna, nc)
            up = not box_at(mid, above=True) or box_at(mid, above=False)
            sign = 1 if up else -1
            tip = y0 + sign * STEM_LOOP_HEIGHT
            ax.plot([mid, mid], [y0, tip - sign * 0.06], color="black", linewidth=1.0,
                    zorder=4, solid_capstyle="butt")
            ax.plot([mid], [tip], marker="o", markersize=5, markerfacecolor="white",
                    markeredgecolor="black", markeredgewidth=1.0, zorder=4)
            text = loop_labels.get(id(nc), "")
            if labels and text:
                ax.text(mid, tip + sign * 0.12, text, ha="center",
                        va="bottom" if up else "top", fontsize=ORF_LABEL_FONTSIZE,
                        color="black")
        return extra

    def _draw_backbone(self, ax: plt.Axes, rna: RNA, pad: int) -> list:
        """The genome line with its end marks (or continuation marks if circular)."""
        settings = self.settings
        seq_len, y0 = rna.length, ANNOTATION_Y_BASE
        if rna.circular:
            ax.plot([0, seq_len], [y0, y0], color="black", linewidth=1.2)
            for x0, x1 in ((-pad, 0), (seq_len, seq_len + pad)):
                ax.plot([x0, x1], [y0, y0], color="black", linewidth=1.2,
                        linestyle=(0, (2, 2)), clip_on=False)
            return [
                ax.text(-pad * 1.25, y0, "\u21ba", va="center", ha="right", fontsize=11),
                ax.text(seq_len + pad * 1.25, y0, "\u21bb", va="center", ha="left", fontsize=11),
            ]

        ax.plot([-pad, seq_len + pad], [y0, y0], color="black", linewidth=1.2)
        ends = []
        five = settings.end_5_label
        if five.strip().lower() == "vpg":
            # ICTV: a grey oval labelled VPg at the 5' end
            oval = Ellipse((-pad * 0.9, y0), width=pad * 1.4, height=0.8,
                           facecolor="0.75", edgecolor="black", linewidth=0.8, clip_on=False)
            ax.add_patch(oval)
            ends.append(ax.text(-pad * 0.9, y0, "VPg", ha="center", va="center",
                                fontsize=7, clip_on=False))
        else:
            ends.append(ax.text(-pad * 0.4, y0, five, va="center", ha="right",
                                fontsize=10, fontweight="bold"))
        ends.append(ax.text(seq_len + pad * 0.4, y0, settings.end_3_label,
                            va="center", ha="left", fontsize=10, fontweight="bold"))
        return ends

    def _arrow_points(self, s_: int, e_: int, yb: float, H: float,
                      forward: bool, seq_len: int) -> list[tuple[float, float]]:
        """A box with an arrowhead at its reading end."""
        head = self._head_length(e_ - s_, seq_len)
        if forward:
            return [(s_, yb), (e_ - head, yb), (e_, yb + H / 2), (e_ - head, yb + H), (s_, yb + H)]
        return [(e_, yb), (s_ + head, yb), (s_, yb + H / 2), (s_ + head, yb + H), (e_, yb + H)]

    def _draw_depth(self, ax: plt.Axes, rna: RNA,
                    tracks: list[np.ndarray], total: np.ndarray) -> list:
        """Single line+fill, or stacked areas for several tracks.

        Returns the layer artists for legend construction.
        """
        s = self.settings
        x = rna.positions

        if len(tracks) == 1:
            layers = [ax.fill_between(x, total, color=s.depth_line_color,
                                      alpha=DEPTH_FILL_ALPHA, linewidth=0)]
        else:
            k = len(tracks)
            colors = (s.stacked_area_colors[:k]
                      + [s.default_color] * max(0, k - len(s.stacked_area_colors)))
            layers = ax.stackplot(x, *reversed(tracks), colors=colors, alpha=0.9, step="pre")
        # the outline follows the fill's shape: stepped when the fill is
        ax.plot(x, total, drawstyle="steps-pre" if len(tracks) > 1 else "default",
                **DEPTH_OUTLINE)

        ax.set_xlim(x[0], x[-1])
        return layers

    def _shade_gaps(self, ax: plt.Axes, threshold_results: list[tuple]) -> None:
        for _, _, gaps, _ in threshold_results:
            for g in gaps:
                ax.axvspan(g["start_bp"], g["end_bp"],
                           color=self.settings.shade_color, alpha=0.15, lw=0)

    def _style_depth_axis(self, ax: plt.Axes, total: np.ndarray,
                          rna: RNA | None = None) -> None:
        args = self.args
        ax.set_ylabel("Read Depth", fontsize=10)
        if self._bare_x(rna):
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Genome Position (bp)", fontsize=10)

        if args.yscale == "symlog":
            ax.set_yscale("symlog", linthresh=args.linthresh, linscale=1)

        ymax = (self.shared_ymax if self.shared_ymax is not None else total.max()) or 1.0
        ax.set_ylim(0, ymax * Y_HEADROOM)

        if args.grid:
            ax.grid(True, linestyle="--", linewidth=0.3)


class CircularPlotter(Plotter):
    """One polar axes: feature arcs in an outer ring, depth as an inner band.

    Position *p* maps to θ = 2π(p−1)/L, with position 1 at the top and
    increasing clockwise, so the whole molecule closes on itself and a
    feature or a read crossing the origin is drawn continuously.

    Three radii depend on the figure being drawn and are set by
    ``_draw_annotations`` / ``_draw_depth`` for the label and axis helpers
    that run after them; they start at the single-lane defaults.
    """

    def __init__(self, settings: Settings, args: argparse.Namespace):
        super().__init__(settings, args)
        self._ring_top = CIRC_R_BASELINE + 0.02 + CIRC_ANN_HEIGHT
        self._outside_label_r = CIRC_R_OUTSIDE_LABEL
        self._depth_ymax = 1.0

    def render(self, rna: RNA, threshold_results: list[tuple], out_base: str) -> None:
        fig = plt.figure(figsize=(CIRC_FIGSIZE, CIRC_FIGSIZE))
        ax = fig.add_subplot(projection="polar")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_rorigin(0)
        ax.set_ylim(0, CIRC_RMAX)          # pin it: autoscale pads below the
                                           # smallest radius and shows an inner spine
        ax.set_yticks([])
        ax.spines["polar"].set_visible(False)

        if self.args.yscale == "symlog":
            log.warning("--yscale symlog is ignored in the circular layout")

        if rna.sgrnas:
            # Deliberate: the ICTV figures put no transcript rows on circular
            # genomes. Geminivirus transcription is bidirectional from the IR
            # with overlapping transcripts rather than a 3'-coterminal set,
            # and nanovirus components carry one ORF each, so the ladder is a
            # linear-genome convention. Nesting near-complete arcs would also
            # be unreadable. See docs/ictv_drawing_conventions.md §6 K.
            log.warning("%d subgenomic RNA row(s) not drawn: the sgRNA ladder is "
                        "a linear-layout convention (use --layout linear)",
                        len(rna.sgrnas))

        extra_artists = self._draw_annotations(ax, rna)

        tracks, total = self._prepared_tracks(rna)
        layers = self._draw_depth(ax, rna, tracks, total)

        if self.args.shade_breaks:
            self._shade_gaps(ax, rna, threshold_results)

        self._draw_axis(ax, rna, total)

        if self.args.legend:
            self._add_legend(ax, layers, rna.labels)

        extra_artists.append(self._add_title(fig, rna))

        self._save(fig, extra_artists, out_base)
        plt.close(fig)

    # --- geometry -----------------------------------------------------------

    @staticmethod
    def _theta(pos, seq_len: int):
        """1-based genome position -> angle in radians."""
        return 2 * np.pi * (np.asarray(pos) - 1) / seq_len

    def _depth_radius(self, values: np.ndarray, ymax: float) -> np.ndarray:
        """Scale depth into the radial band [base, max]."""
        span = CIRC_R_DEPTH_MAX - CIRC_R_DEPTH_BASE
        scaled = np.clip(values / ymax, 0.0, 1.0) if ymax else np.zeros_like(values)
        return CIRC_R_DEPTH_BASE + scaled * span

    @staticmethod
    def _close(theta: np.ndarray, *arrays: np.ndarray):
        """Repeat the first sample at θ=2π so filled areas close at the origin."""
        theta_c = np.append(theta, 2 * np.pi)
        return (theta_c, *(np.append(a, a[0]) for a in arrays))

    # --- panels -------------------------------------------------------------

    def _draw_annotations(self, ax, rna: RNA) -> list:
        """Genome circle, ORF arcs with arrowheads, origin stem-loop.

        Modes C1/C2 of docs/ictv_drawing_conventions.md: arcs sit just outside
        the circle, run clockwise for virion-sense ORFs and anticlockwise for
        complementary-sense ones (so each half of a geminivirus circle carries
        one strand), and an arc overlapping one already placed nests inward.

        LIMIT: ``Feature.mechanism`` is not drawn here. A frameshift or
        readthrough ORF is nested like any other arc, but the ``+1 FS`` /
        ``RT`` marks that ``LinearPlotter._draw_mechanism`` writes have no
        circular counterpart yet — the junction would need a radial tick and a
        label placed clear of the arc labels. Rare on circular genomes (the
        ICTV circular figures show none), so deferred.
        """
        args = self.args
        extra: list = []
        L = rna.length

        placements, _ = self._placements(rna, circular=True)
        deepest = max((p.tier for p in placements), default=0)
        if deepest >= CIRC_MAX_TIERS:
            log.warning("%d nesting levels needed but only %d fit; deepest arcs are "
                        "drawn on the innermost lane", deepest + 1, CIRC_MAX_TIERS)
        deepest = min(deepest, CIRC_MAX_TIERS - 1)
        # The innermost tier in use is anchored to the genome circle; tier 0 is
        # the outermost, so "AC4 inside AC1" reads as the ICTV figures draw it.
        #   outside   : arcs sit just outside the circle (default; the depth
        #               ring occupies the inside, so nesting cannot go there)
        #   on_circle : the innermost arcs straddle the circle line, as in the
        #               Begomovirus / Mastrevirus figures where the arcs *are*
        #               the genome
        if self.settings.circular_arcs == "on_circle":
            innermost = CIRC_R_BASELINE - CIRC_ANN_HEIGHT / 2
        else:
            innermost = CIRC_R_BASELINE + 0.02
        lane0 = innermost + deepest * CIRC_LANE_STEP
        self._ring_top = lane0 + CIRC_ANN_HEIGHT
        self._outside_label_r = self._ring_top + 0.03

        # the genome circle, under the arcs
        theta = np.linspace(0, 2 * np.pi, 361)
        ax.plot(theta, np.full_like(theta, CIRC_R_BASELINE), color="black",
                linewidth=1.0, zorder=2)

        # non-coding regions: a thick grey arc astride the circle (Begomovirus
        # CRA), under the ORF arcs; stem-loops are drawn after the ring below
        regions = [nc for nc in rna.noncoding if nc.kind != "stem_loop"]
        for nc in regions:
            r_in = CIRC_R_BASELINE - CIRC_NONCODING_HEIGHT / 2
            for (s_, e_) in rna.feature_spans(nc):
                th, r = self._arc_arrow(self._theta(s_, L), self._theta(e_ + 1, L),
                                        r_in, CIRC_NONCODING_HEIGHT, True, head=0.0)
                ax.fill(th, r, facecolor=self._noncoding_color(nc), zorder=2.5,
                        **self._glyph_edge())

        for pl in placements:
            feat = pl.feature
            tier = min(pl.tier, CIRC_MAX_TIERS - 1)
            lane = lane0 - tier * CIRC_LANE_STEP
            color = self._feature_color(feat.product)
            spans = rna.feature_spans(feat)
            head_span = spans[-1] if feat.forward else spans[0]
            for (s_, e_) in spans:
                th0, th1 = self._theta(s_, L), self._theta(e_ + 1, L)
                head = self._theta(1 + self._head_length(e_ - s_ + 1, L), L) if (s_, e_) == head_span else 0.0
                th, r = self._arc_arrow(th0, th1, lane, CIRC_ANN_HEIGHT, feat.forward, head)
                (glyph,) = ax.fill(th, r, facecolor=color, zorder=3, **self._glyph_edge())
                if feat.domains:
                    extra.extend(self._draw_domains(ax, rna, feat, s_, e_, lane, tier, glyph,
                                                    labels=not args.no_label))

            if not args.no_label:
                horizontal = self.settings.circular_labels == "horizontal"
                if feat.domains:                     # the domains have labelled the arc
                    if feat.gene:
                        extra.append(self._label_horizontal(ax, rna, feat.gene, spans, tier)
                                     if horizontal else
                                     self._label_outside(ax, rna, feat.gene, spans,
                                                         fontsize=ORF_LABEL_FONTSIZE))
                elif horizontal:
                    text = f"{feat.gene} ({feat.product})" if feat.gene else feat.product
                    extra.append(self._label_horizontal(ax, rna, text, spans, tier))
                else:
                    extra.append(self._label_arc(ax, rna, feat.product, lane, spans, fill=color))
                    if feat.gene:
                        extra.append(self._label_outside(ax, rna, feat.gene, spans,
                                                         fontsize=ORF_LABEL_FONTSIZE))

        # position 1 at 12 o'clock: a dashed radius through the depth band
        ax.plot([0, 0], [CIRC_R_DEPTH_BASE, CIRC_R_BASELINE], color="black",
                linewidth=0.8, linestyle=(0, (3, 2)), zorder=4)

        # non-coding labels, then stem-loop icons. A GFF stem_loop puts the
        # hairpin where it belongs (the nick site is rarely exactly at 1);
        # with none, the icon marks the origin as before, unlabelled.
        loops = [nc for nc in rna.noncoding if nc.kind == "stem_loop"]
        loop_labels, region_labels = self._noncoding_labels(rna, regions, loops)
        if not args.no_label:
            for nc, text in region_labels:
                spans = rna.feature_spans(nc)
                extra.append(self._label_horizontal(ax, rna, text, spans)
                             if self.settings.circular_labels == "horizontal" else
                             self._label_arc(ax, rna, text, CIRC_R_BASELINE - CIRC_ANN_HEIGHT / 2,
                                             spans, fill=self._noncoding_color(nc), dark_below=0.5))
        for nc in loops or [None]:
            theta_sl = 0.0 if nc is None else float(self._theta(self._mid(rna, nc), L))
            text = loop_labels.get(id(nc), "") if nc is not None and not args.no_label else ""
            extra.append(self._stem_loop_icon(ax, theta_sl, text))
        return [a for a in extra if a is not None]


    def _stem_loop_icon(self, ax, theta: float, label: str = ""):
        """A hairpin standing on the outside of the ring at ``theta``, optionally named."""
        stem_top = self._ring_top + 0.03
        ax.plot([theta, theta], [stem_top, stem_top + 0.04], color="black", linewidth=1.0,
                zorder=4, clip_on=False)
        ax.plot([theta], [stem_top + 0.055], marker="o", markersize=5, markerfacecolor="white",
                markeredgecolor="black", markeredgewidth=1.0, zorder=4, clip_on=False)
        if not label:
            return None
        ha, va = self._level_anchor(theta)
        return ax.annotate(label, xy=(theta, stem_top + 0.085), xytext=(0, 0),
                           textcoords="offset points", ha=ha, va=va,
                           fontsize=self.settings.annotation_fontsize,
                           annotation_clip=False, zorder=5)

    @staticmethod
    def _level_anchor(theta: float) -> tuple[str, str]:
        """Text anchor for level text placed outside the ring at ``theta``."""
        deg = np.degrees(theta) % 360                     # 0 = top, clockwise
        if deg < 15 or deg > 345:
            return "center", "bottom"
        if 165 < deg < 195:
            return "center", "top"
        return ("left", "center") if deg < 180 else ("right", "center")

    @staticmethod
    def _arc_arrow(th0: float, th1: float, lane: float, height: float,
                   forward: bool, head: float) -> tuple[np.ndarray, np.ndarray]:
        """One closed polar path for an arc with an optional arrowhead.

        Body and head used to be two patches, which left a hairline seam and a
        flat head base against a curved wedge. A single densely sampled
        outline (outer arc, head, inner arc back) has neither.
        """
        r_in, r_out = lane, lane + height
        r_lo, r_hi, r_mid = lane - height * 0.18, lane + height * 1.18, lane + height / 2
        head = min(head, th1 - th0)
        b0, b1 = (th0, th1 - head) if forward else (th0 + head, th1)
        n = max(4, int(np.degrees(b1 - b0) * 2))          # ~2 samples per degree
        outer = np.linspace(b0, b1, n)
        inner = outer[::-1]
        if head <= 0:
            th = np.concatenate([outer, inner])
            r = np.concatenate([np.full(n, r_out), np.full(n, r_in)])
        elif forward:                                      # tip at th1 (clockwise end)
            th = np.concatenate([outer, [b1, th1, b1], inner])
            r = np.concatenate([np.full(n, r_out), [r_hi, r_mid, r_lo], np.full(n, r_in)])
        else:                                              # tip at th0 (anticlockwise end)
            th = np.concatenate([outer, inner, [b0, th0, b0]])
            r = np.concatenate([np.full(n, r_out), np.full(n, r_in), [r_lo, r_mid, r_hi]])
        return th, r

    def _label_horizontal(self, ax, rna: RNA, text: str, spans, tier: int = 0):
        """Level text just outside the ring, anchored on the side facing the arc.

        The ICTV circular figures write every label horizontally, outside the
        ring, in the form ``AV1 (CP)``; nested arcs are labelled inside the
        circle there, which the depth ring rules out here, so all go outside.
        A nested arc shares its mid-angle with the arc it nests in, so its
        label is stepped one line away from the circle's equator per tier
        (upwards in the top half, downwards in the bottom half) instead of
        overprinting.
        """
        covered = sum(e - s + 1 for s, e in spans)
        mid = ((spans[0][0] - 1 + covered / 2) % rna.length) + 1
        theta = float(self._theta(mid, rna.length))
        deg = np.degrees(theta) % 360                     # 0 = top, clockwise
        ha, va = self._level_anchor(theta)
        line = self.settings.annotation_fontsize * 1.4          # points
        dy = tier * line * (1 if deg < 90 or deg > 270 else -1)
        return ax.annotate(text, xy=(theta, self._outside_label_r + 0.01),
                           xytext=(0, dy), textcoords="offset points", ha=ha, va=va,
                           fontsize=self.settings.annotation_fontsize,
                           annotation_clip=False, zorder=5)

    def _label_outside(self, ax, rna: RNA, text: str, spans, fontsize: int):
        """Radial text just outside the arc ring, reading outwards."""
        covered = sum(e - s + 1 for s, e in spans)
        mid = ((spans[0][0] - 1 + covered / 2) % rna.length) + 1
        theta = float(self._theta(mid, rna.length))
        deg = np.degrees(theta)
        left = 90 < (deg % 360) < 270
        return ax.text(theta, self._outside_label_r, text,
                       ha="right" if left else "left", va="center",
                       rotation=(-deg + 180) if left else -deg,
                       rotation_mode="anchor", fontsize=fontsize, color="0.25",
                       clip_on=False, zorder=5)

    def _draw_domains(self, ax, rna: RNA, feat, s_: int, e_: int, lane: float, tier: int,
                      glyph, labels: bool) -> list:
        """Split a polyprotein arc into its mature proteins: own colours, thin
        radial dividers, one label per domain by the layout's label rule.
        Mirrors ``LinearPlotter._draw_domains``; see that docstring."""
        L, H = rna.length, CIRC_ANN_HEIGHT
        out = []
        for d in feat.domains_within(s_, e_):
            th, r = self._arc_arrow(self._theta(d.start, L), self._theta(d.end + 1, L),
                                    lane, H, feat.forward, head=0.0)
            (seg,) = ax.fill(th, r, facecolor=self._feature_color(d.product), zorder=3,
                             **self._glyph_edge())
            seg.set_clip_path(glyph)
            if labels:
                out.append(self._label_horizontal(ax, rna, d.product, [(d.start, d.end)], tier)
                           if self.settings.circular_labels == "horizontal" else
                           self._label_arc(ax, rna, d.product, lane, [(d.start, d.end)],
                                           fill=self._feature_color(d.product)))
        for x in feat.dividers_within(s_, e_):
            th = float(self._theta(x + 0.5, L))
            (line,) = ax.plot([th, th], [lane, lane + H], zorder=4, **DOMAIN_DIVIDER)
            line.set_clip_path(glyph)
        return out

    def _label_arc(self, ax, rna: RNA, text: str, lane: float,
                   spans: list[tuple[int, int]], fill: str = "white",
                   dark_below: float = 0.35):
        """Label along the arc (in a colour that reads on ``fill``), or
        radially outside it when the arc is narrow."""
        fontsize = self.settings.annotation_fontsize
        covered = sum(e - s + 1 for s, e in spans)
        # midway along the covered arc, which may run through the origin
        mid = ((spans[0][0] - 1 + covered / 2) % rna.length) + 1
        theta = float(self._theta(mid, rna.length))
        deg = np.degrees(theta)
        fraction = min(covered, rna.length) / rna.length

        if fraction >= CIRC_SMALL_FEATURE_FRACTION:
            rotation = -deg
            if 90 < (deg % 360) < 270:          # keep text the right way up
                rotation += 180
            return ax.text(theta, lane + CIRC_ANN_HEIGHT / 2, text,
                           ha="center", va="center", rotation=rotation,
                           rotation_mode="anchor", fontsize=fontsize, zorder=5,
                           color=self._text_on(fill, dark_below))

        # narrow feature: read outwards from the rim
        left = 90 < (deg % 360) < 270
        return ax.text(theta, self._outside_label_r, text,
                       ha="right" if left else "left", va="center",
                       rotation=(-deg + 180) if left else -deg,
                       rotation_mode="anchor", fontsize=fontsize,
                       clip_on=False, zorder=5)

    def _draw_depth(self, ax, rna: RNA, tracks: list[np.ndarray],
                    total: np.ndarray) -> list:
        """Depth as a filled radial band; several tracks stack cumulatively."""
        s = self.settings
        ymax = (self.shared_ymax if self.shared_ymax is not None else total.max()) or 1.0
        ymax *= Y_HEADROOM
        theta = self._theta(rna.positions, rna.length)
        base = CIRC_R_DEPTH_BASE

        ax.plot(*self._close(theta, np.full_like(total, base, dtype=float)),
                color="black", linewidth=0.5, alpha=0.4, zorder=1)

        if len(tracks) == 1:
            r = self._depth_radius(total, ymax)
            th, rc = self._close(theta, r)
            layers = [ax.fill_between(th, base, rc, color=s.depth_line_color,
                                      alpha=DEPTH_FILL_ALPHA, linewidth=0, zorder=2)]
            ax.plot(th, rc, zorder=3, **DEPTH_OUTLINE)
        else:
            k = len(tracks)
            colors = (s.stacked_area_colors[:k]
                      + [s.default_color] * max(0, k - len(s.stacked_area_colors)))
            layers = []
            lower = np.zeros_like(total, dtype=float)
            for track, color in zip(reversed(tracks), colors):
                upper = lower + track
                th, rl, ru = self._close(theta,
                                         self._depth_radius(lower, ymax),
                                         self._depth_radius(upper, ymax))
                layers.append(ax.fill_between(th, rl, ru, color=color, alpha=0.9, zorder=2))
                lower = upper
            th, rc = self._close(theta, self._depth_radius(total, ymax))
            ax.plot(th, rc, zorder=3, **DEPTH_OUTLINE)

        self._depth_ymax = ymax
        return layers

    def _shade_gaps(self, ax, rna: RNA, threshold_results: list[tuple]) -> None:
        height = CIRC_R_DEPTH_MAX - CIRC_R_DEPTH_BASE
        for _, _, gaps, _ in threshold_results:
            for g in gaps:
                start, end = g["start_bp"], g["end_bp"]
                width = self._theta(end + 1, rna.length) - self._theta(start, rna.length)
                centre = self._theta(start, rna.length) + width / 2
                ax.bar(centre, height, width=width, bottom=CIRC_R_DEPTH_BASE,
                       color=self.settings.shade_color, alpha=0.15, linewidth=0, zorder=1)

    def _draw_axis(self, ax, rna: RNA, total: np.ndarray) -> None:
        """Genome-position ticks round the rim, and a depth scale at the origin."""
        step = _nice_step(rna.length / 8)
        ticks = np.arange(0, rna.length, step)
        ax.set_xticks(self._theta(ticks + 1, rna.length))
        ax.set_xticklabels([f"{int(t):,}" if t else "1" for t in ticks], fontsize=8)
        ax.tick_params(axis="x", pad=2)
        ax.grid(self.args.grid, axis="x", linestyle="--", linewidth=0.3, alpha=0.6)

        # depth scale: floor and ceiling on the origin radius, so the radial
        # extent is readable without a second axis
        ymax = self._depth_ymax
        fmt = (lambda v: f"{v:.2g}") if self.args.normalize else (lambda v: f"{v:.0f}")
        unit = "" if self.args.normalize else "\u00d7"
        # the ceiling label normally stands on its radius, up against the
        # genome circle; when a non-coding arc straddles the origin (a
        # geminivirus common region) it hangs inward instead, clear of the arc
        near = 0.03 * rna.length
        covered = any(s_ <= 1 + near or e_ >= rna.length - near
                      for nc in rna.noncoding if nc.kind != "stem_loop"
                      for s_, e_ in rna.feature_spans(nc))
        for radius, text, va in ((CIRC_R_DEPTH_BASE, f" {fmt(0.0)}", "bottom"),
                                 (CIRC_R_DEPTH_MAX, f" {fmt(ymax)}{unit} depth",
                                  "top" if covered else "bottom")):
            ax.text(0, radius, text, ha="left", va=va,
                    fontsize=7, color="0.35", zorder=5)

        # the hub is empty: name the molecule there rather than over the trace
        ax.text(0, 0, f"{rna.name}\n{rna.length:,} nts", ha="center", va="center",
                fontsize=9, color="0.35", linespacing=1.5, zorder=5)


def _nice_step(target: float) -> int:
    """Round a tick spacing up to 1, 2 or 5 x a power of ten."""
    if target <= 0:
        return 1
    exp = 10 ** int(np.floor(np.log10(target)))
    for mult in (1, 2, 5, 10):
        if mult * exp >= target:
            return max(1, int(mult * exp))
    return max(1, int(10 * exp))
