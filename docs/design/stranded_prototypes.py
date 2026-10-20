"""Prototype figure for docs/stranded_depth.md. Not part of VirPlot.

Writes docs/stranded_prototypes.svg (docs/*.svg is tracked; docs/figures/
is emptied by make_figures.sh on every run, so the figure does not live there).

Splits the real GLRaV-13 H8881 reads (examples/glrav13.bam) by strand of
origin and draws them six ways: two linear (mirrored, linear and symlog) and
four circular (A-D). GLRaV-13 is a linear genome; the circular panels use its
depth only to test the geometry with a real strand ratio.

    PYTHONPATH=src python3 docs/design/stranded_prototypes.py
"""

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from virplot.alignments import SKIP_FLAGS, open_alignments
from virplot.plotting import CIRC_R_BASELINE, CIRC_R_DEPTH_BASE, CIRC_R_DEPTH_MAX

HERE = os.path.dirname(os.path.abspath(__file__))
BAM = os.path.join(HERE, "..", "..", "examples", "glrav13.bam")
OUT = os.path.normpath(os.path.join(HERE, "..", "stranded_prototypes.svg"))

PC, MC = "#7fc8d6", "#d68a7f"                  # + strand, - strand
LINE = dict(color="black", linewidth=0.3, alpha=0.8)
LINTHRESH = 10.0


def strand_depth(path: str, length: int) -> tuple[np.ndarray, np.ndarray]:
    """(+ strand, - strand) depth for a reverse-stranded single-end library.

    The read is antisense to the RNA it came from, so a read aligned in
    reverse (flag 0x10) counts towards the + strand.
    """
    diff = {True: np.zeros(length + 1, np.int64), False: np.zeros(length + 1, np.int64)}
    for aln in open_alignments(path).records:
        if aln.flag & SKIP_FLAGS:
            continue
        d = diff[bool(aln.flag & 0x10)]
        pos = aln.pos - 1
        for n, op in aln.cigar:
            if op in "M=X":
                s, e = max(pos, 0), min(pos + n, length)
                if s < e:
                    d[s] += 1
                    d[e] -= 1
            if op in "MDN=X":
                pos += n
    return np.cumsum(diff[True][:-1]), np.cumsum(diff[False][:-1])


def symlog(v):
    v = np.asarray(v, dtype=float)
    return np.where(v < LINTHRESH, v / LINTHRESH, 1 + np.log10(np.maximum(v, LINTHRESH) / LINTHRESH))


def main() -> None:
    length = 17564
    plus, minus = strand_depth(BAM, length)
    k = np.ones(15) / 15
    plus, minus = np.convolve(plus, k, "same"), np.convolve(minus, k, "same")
    # every 16th base after smoothing (the 15-base window already set the detail)
    step = 16
    plus, minus = plus[::step], minus[::step]
    x = np.arange(1, length + 1, step)
    theta = np.append(2 * np.pi * (x - 1) / length, 2 * np.pi)
    n = len(x)
    close = lambda a: np.append(a, a[0])
    base, top = CIRC_R_DEPTH_BASE, CIRC_R_DEPTH_MAX
    mid, half = (base + top) / 2, (top - base) / 2

    fig = plt.figure(figsize=(13, 9.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.45], hspace=0.3, wspace=0.3)

    for col, (title, log) in enumerate((("Linear, mirrored, shared linear scale", False),
                                        ("Linear, mirrored, shared symlog scale", True))):
        ax = fig.add_subplot(gs[0, 2 * col:2 * col + 2])
        ax.fill_between(x, plus, color=PC, lw=0, rasterized=True)
        ax.fill_between(x, -minus, color=MC, lw=0, rasterized=True)
        ax.plot(x, plus, rasterized=True, **LINE)
        ax.plot(x, -minus, rasterized=True, **LINE)
        ax.axhline(0, color="black", lw=0.6)
        if log:
            ax.set_yscale("symlog", linthresh=LINTHRESH)
        ax.set_xlim(1, length)
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_ylabel("depth (+ up, − down)", fontsize=8)
        ax.tick_params(labelsize=7)

    def polar(col, title):
        ax = fig.add_subplot(gs[1, col], projection="polar")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_ylim(0, 0.78)
        ax.axis("off")
        ax.plot(theta, np.full(n + 1, CIRC_R_BASELINE), color="black", lw=1.2)
        ax.set_title(title, fontsize=9, pad=2)
        return ax

    def band(ax, r0, r1, color):
        # polar paths are subdivided per segment; raster bands keep the SVG small
        ax.fill_between(theta, close(r0), close(r1), color=color, lw=0, zorder=2, rasterized=True)
        ax.plot(theta, close(r1), zorder=3, rasterized=True, **LINE)

    flat = np.full(n, mid)
    shared = max(plus.max(), minus.max()) * 1.05

    ax = polar(0, "A. Mirrored, shared linear scale")
    band(ax, flat, mid + plus / shared * half, PC)
    band(ax, flat, mid - minus / shared * half, MC)

    ax = polar(1, f"B. Mirrored, own scale per strand\n(− magnified {plus.max() / minus.max():.1f}×)")
    band(ax, flat, mid + plus / (plus.max() * 1.05) * half, PC)
    band(ax, flat, mid - minus / (minus.max() * 1.05) * half, MC)

    ax = polar(2, "C. Two bands, both outward,\nown scales (+ outer, − inner)")
    gap = 0.02
    hb = (top - base - gap) / 2
    band(ax, np.full(n, base), base + minus / (minus.max() * 1.05) * hb, MC)
    b2 = base + hb + gap
    band(ax, np.full(n, b2), b2 + plus / (plus.max() * 1.05) * hb, PC)

    ax = polar(3, "D. Mirrored, shared symlog scale\n(dotted: 10×, 100×, 1000×)")
    s_top = symlog(shared)
    band(ax, flat, mid + symlog(plus) / s_top * half, PC)
    band(ax, flat, mid - symlog(minus) / s_top * half, MC)
    for d in (10, 100, 1000):
        for sign in (1, -1):
            ax.plot(theta, np.full(n + 1, mid + sign * symlog(d) / s_top * half),
                    color="0.45", lw=0.35, ls=":", zorder=4)

    fig.suptitle("Stranded depth prototypes: GLRaV-13 H8881, real reverse-stranded reads "
                 "(+ blue, − red; circular panels test geometry only)", fontsize=10, y=0.995)
    plt.rcParams["svg.hashsalt"] = "virplot"     # stable ids, so reruns give the same file
    plt.rcParams["svg.fonttype"] = "none"        # text as text, not glyph paths
    fig.savefig(OUT, bbox_inches="tight", dpi=150, metadata={"Date": None})
    print(f"wrote {OUT}  (+ mean {plus.mean():.0f}x, - mean {minus.mean():.0f}x)")


if __name__ == "__main__":
    main()
