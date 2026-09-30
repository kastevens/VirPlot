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
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

from virplot.analysis import smooth_depth
from virplot.models import RNA
from virplot.settings import Settings

log = logging.getLogger(__name__)

# layout constants

ANNOTATION_Y_BASE = 0.5
ANNOTATION_HEIGHT = 0.6
LABEL_PAD_FRACTION = 0.02
SMALL_FEATURE_THRESHOLD = 500  # bp; features shorter than this get external labels
SMOOTH_WINDOW = 15
Y_HEADROOM = 1.05
PNG_DPI = 400
FIGURE_WIDTH = 12
FIGURE_HEIGHT = 4
HEIGHT_RATIO_ANNOTATION = 1
HEIGHT_RATIO_DEPTH = 1

# circular layout: radii as a fraction of the plotted radius, inside out
CIRC_FIGSIZE = 8
CIRC_R_DEPTH_BASE = 0.30        # depth baseline (zero coverage)
CIRC_R_DEPTH_MAX = 0.72         # depth at the y-limit
CIRC_R_ANN_LANES = (0.80, 0.895)  # bottom radius of the inner / outer feature lane
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

    def prepare(self, rnas: list[RNA]) -> None:
        """Compute cross-RNA scaling so separate figures share axes conventions."""
        if len(rnas) < 2:
            return
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
        tracks = rna.tracks
        if self.args.smooth:
            tracks = [smooth_depth(y, window_size=SMOOTH_WINDOW) for y in tracks]
        return np.sum(tracks, axis=0) if len(tracks) > 1 else tracks[0]

    def _prepared_tracks(self, rna: RNA) -> tuple[list[np.ndarray], np.ndarray]:
        """Apply --smooth / --normalize to the depth tracks.

        Returns (per-track arrays, combined array). Smoothing is applied per
        track; normalisation divides everything by the combined maximum (of
        this RNA, or of all RNAs after ``prepare``) so the total peaks at 1.
        """
        tracks = rna.tracks
        if self.args.smooth:
            tracks = [smooth_depth(y, window_size=SMOOTH_WINDOW) for y in tracks]
        total = np.sum(tracks, axis=0) if len(tracks) > 1 else tracks[0]
        if self.args.normalize:
            denom = self.shared_denom or total.max() or 1.0
            tracks = [y / denom for y in tracks]
            total = total / denom
        return tracks, total

    def _figure_width(self, rna: RNA) -> float:
        if self.max_length:
            return FIGURE_WIDTH * rna.length / self.max_length
        return FIGURE_WIDTH

    def _feature_color(self, product: str) -> str:
        return self.settings.color_mapping.get(product, self.settings.default_color)

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

    def _add_title(self, fig: plt.Figure):
        """Draw the YAML title if --title was given; returns the artist or None."""
        if not self.args.title:
            return None
        return fig.text(
            0.5, 0.95, self.settings.title,
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

        self._style_depth_axis(ax_depth, total)

        extra_artists.append(self._add_title(fig))

        self._save(fig, extra_artists, out_base)
        plt.close(fig)

    # --- panels -------------------------------------------------------------

    def _draw_annotations(self, ax: plt.Axes, rna: RNA) -> list:
        """Genome line, 5'/3' marks and feature rectangles.

        Returns the end-mark artists so the saved bbox includes them.
        """
        args = self.args
        seq_len = rna.length
        pad = int(seq_len * LABEL_PAD_FRACTION)
        y0 = ANNOTATION_Y_BASE

        if rna.circular:
            # No 5'/3' ends on a circle: show the backbone continuing past both
            # edges instead, so the join at the origin is visible.
            ax.plot([0, seq_len], [y0, y0], color="black", linewidth=1.2)
            for x0, x1 in ((-pad, 0), (seq_len, seq_len + pad)):
                ax.plot([x0, x1], [y0, y0], color="black", linewidth=1.2,
                        linestyle=(0, (2, 2)), clip_on=False)
            ends = [
                ax.text(-pad * 1.25, y0, "\u21ba", va="center", ha="right", fontsize=11),
                ax.text(seq_len + pad * 1.25, y0, "\u21bb", va="center", ha="left", fontsize=11),
            ]
        else:
            ax.plot([-pad, seq_len + pad], [y0, y0], color="black", linewidth=1.2)
            ends = [
                ax.text(-pad * 0.4, y0, "5'", va="center", ha="right",
                        fontsize=10, fontweight="bold"),
                ax.text(seq_len + pad * 0.4, y0, "3'", va="center", ha="left",
                        fontsize=10, fontweight="bold"),
            ]

        for i, feat in enumerate(rna.features):
            upper = i % 2 == 0
            y = ANNOTATION_Y_BASE if upper else ANNOTATION_Y_BASE - ANNOTATION_HEIGHT
            rect = Rectangle(
                (feat.start, y), feat.length, ANNOTATION_HEIGHT,
                facecolor=self._feature_color(feat.product),
                edgecolor="none" if args.no_border else "black",
            )
            ax.add_patch(rect)

            if not args.no_label:
                if feat.length < SMALL_FEATURE_THRESHOLD:
                    ly = (y + ANNOTATION_HEIGHT * 1.5 if upper
                          else y - ANNOTATION_HEIGHT * 0.5)
                else:
                    ly = y + ANNOTATION_HEIGHT / 2
                ax.text(feat.midpoint, ly, feat.product, ha="center", va="center",
                        fontsize=self.settings.annotation_fontsize, color="black")

        ax.set_xlim(0, seq_len)
        ax.set_ylim(-1.5, 2.0)
        ax.axis("off")
        return ends

    def _draw_depth(self, ax: plt.Axes, rna: RNA,
                    tracks: list[np.ndarray], total: np.ndarray) -> list:
        """Single line+fill, or stacked areas for several tracks.

        Returns the layer artists for legend construction.
        """
        s = self.settings
        x = rna.positions

        if len(tracks) == 1:
            layers = ax.plot(x, total, color=s.depth_line_color, linewidth=0.8, alpha=0.9)
            ax.fill_between(x, total, color=s.depth_line_color, alpha=0.3)
        else:
            k = len(tracks)
            colors = (s.stacked_area_colors[:k]
                      + [s.default_color] * max(0, k - len(s.stacked_area_colors)))
            layers = ax.stackplot(x, *reversed(tracks), colors=colors, alpha=0.9, step="pre")
            ax.plot(x, total, color="black", linewidth=0.3, alpha=0.8, label="Combined depth")

        ax.set_xlim(x[0], x[-1])
        return layers

    def _shade_gaps(self, ax: plt.Axes, threshold_results: list[tuple]) -> None:
        for _, _, gaps, _ in threshold_results:
            for g in gaps:
                ax.axvspan(g["start_bp"], g["end_bp"],
                           color=self.settings.shade_color, alpha=0.15, lw=0)

    def _style_depth_axis(self, ax: plt.Axes, total: np.ndarray) -> None:
        args = self.args
        ax.set_ylabel("Read Depth", fontsize=10)
        ax.set_xlabel("Genome Position (bp)", fontsize=10)

        if args.yscale == "symlog":
            ax.set_yscale("symlog", linthresh=args.linthresh, linscale=1)

        ymax = self.shared_ymax if self.shared_ymax is not None else total.max()
        ax.set_ylim(0, ymax * Y_HEADROOM if total.size else 1)

        if args.grid:
            ax.grid(True, linestyle="--", linewidth=0.3)


class CircularPlotter(Plotter):
    """One polar axes: feature arcs in an outer ring, depth as an inner band.

    Position *p* maps to θ = 2π(p−1)/L, with position 1 at the top and
    increasing clockwise, so the whole molecule closes on itself and a
    feature or a read crossing the origin is drawn continuously.
    """

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

        extra_artists = self._draw_annotations(ax, rna)

        tracks, total = self._prepared_tracks(rna)
        layers = self._draw_depth(ax, rna, tracks, total)

        if self.args.shade_breaks:
            self._shade_gaps(ax, rna, threshold_results)

        self._draw_axis(ax, rna, total)

        if self.args.legend:
            self._add_legend(ax, layers, rna.labels)

        extra_artists.append(self._add_title(fig))

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
        """Feature arcs in two lanes, plus the origin marker."""
        args = self.args
        extra: list = []

        for i, feat in enumerate(rna.features):
            lane = CIRC_R_ANN_LANES[i % 2]
            color = self._feature_color(feat.product)
            for start, end in rna.feature_spans(feat):
                width = self._theta(end + 1, rna.length) - self._theta(start, rna.length)
                centre = self._theta(start, rna.length) + width / 2
                ax.bar(centre, CIRC_ANN_HEIGHT, width=width, bottom=lane,
                       facecolor=color, edgecolor="none" if args.no_border else "black",
                       linewidth=0.6, align="center", zorder=3)

            if not args.no_label:
                extra.append(self._label_feature(ax, rna, feat, lane,
                                                 rna.feature_spans(feat)))

        # origin: a radial tick at position 1
        ax.plot([0, 0], [CIRC_R_DEPTH_BASE, CIRC_R_ANN_LANES[1] + CIRC_ANN_HEIGHT],
                color="black", linewidth=0.8, linestyle=(0, (3, 2)), zorder=4)
        return [a for a in extra if a is not None]

    def _label_feature(self, ax, rna: RNA, feat, lane: float,
                       spans: list[tuple[int, int]]):
        """Label along the arc, or radially outside it when the arc is narrow."""
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
            return ax.text(theta, lane + CIRC_ANN_HEIGHT / 2, feat.product,
                           ha="center", va="center", rotation=rotation,
                           rotation_mode="anchor", fontsize=fontsize, zorder=5)

        # narrow feature: read outwards from the rim
        left = 90 < (deg % 360) < 270
        return ax.text(theta, CIRC_R_OUTSIDE_LABEL, feat.product,
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
            ax.fill_between(th, base, rc, color=s.depth_line_color, alpha=0.3, zorder=2)
            layers = ax.plot(th, rc, color=s.depth_line_color, linewidth=0.8,
                             alpha=0.9, zorder=2)
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
            ax.plot(th, rc, color="black", linewidth=0.3, alpha=0.8, zorder=3,
                    label="Combined depth")

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
        ymax = getattr(self, "_depth_ymax", total.max() or 1.0)
        fmt = (lambda v: f"{v:.2g}") if self.args.normalize else (lambda v: f"{v:.0f}")
        unit = "" if self.args.normalize else "\u00d7"
        for radius, text in ((CIRC_R_DEPTH_BASE, f" {fmt(0.0)}"),
                             (CIRC_R_DEPTH_MAX, f" {fmt(ymax)}{unit} depth")):
            ax.text(0, radius, text, ha="left", va="bottom",
                    fontsize=7, color="0.35", zorder=5)

        # the hub is empty: name the molecule there rather than over the trace
        ax.text(0, 0, f"{rna.name}\n{rna.length:,} nt", ha="center", va="center",
                fontsize=9, color="0.35", linespacing=1.5, zorder=5)


def _nice_step(target: float) -> int:
    """Round a tick spacing up to 1, 2 or 5 x a power of ten."""
    if target <= 0:
        return 1
    exp = 10 ** int(np.floor(np.log10(target)))
    for mult in (1, 2, 5, 10):
        if mult * exp >= target:
            return int(mult * exp)
    return int(10 * exp)
