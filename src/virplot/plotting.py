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
