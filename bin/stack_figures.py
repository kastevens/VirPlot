#!/usr/bin/env python3
"""Stack VirPlot segment figures into one ICTV-style multipartite figure.

The ICTV convention for a segmented genome is one row per segment, stacked
largest first on a shared scale (Geminiviridae DNA-A/DNA-B, LIYV RNA-1/RNA-2,
TSWV L/M/S, nanovirus components). VirPlot renders one file per RNA rather
than one stacked figure, because it already renders them *composably*: within
a run the figure width is proportional to genome length and matplotlib's
margins are fractional, so the panels share an exact x-scale (measured across
TSWV's 2,916-8,897 bp spread: 0.0000% scale drift). Stacking is therefore a
document operation, not a re-plot, and this script is it.

What it does beyond pasting:

* orders panels largest first, as the convention asks;
* **aligns the plot areas**. Panels can start at slightly different x when
  their y-axis tick labels differ in width (a sample peaking at 9x and one
  peaking at 1500x), because the figures are saved with a tight bounding box.
  Each panel is shifted so every plot area begins at the same x, which hand
  pasting cannot do;
* leaves the panels themselves untouched - each is nested in its own <svg>
  with its own coordinate system, so nothing is rescaled or redrawn.

Render the segments in ONE virplot run so they share a scale and a depth
y-limit, then stack them:

    virplot -g tswv.gff3 -d tswv.sam -y tswv.yml -f svg -o out --name tswv
    python3 bin/stack_figures.py -i out/tswv.*.svg -o out/tswv_stacked.svg

Pass --bare-x to virplot to drop the repeated x-axis from all but the
bottom panel.
"""

import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET

SVG_NS = "http://www.w3.org/2000/svg"
SVG = "{%s}" % SVG_NS
_PT = re.compile(r"([\d.]+)pt")
_RECT = re.compile(r'<rect x="([\d.]+)" y="([\d.]+)" '
                   r'width="([\d.]+)" height="([\d.]+)"')

ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")


class StackError(ValueError):
    """Raised for a file that is not a VirPlot SVG panel."""


class Panel:
    """One rendered segment: its page size and where its plot area sits."""

    def __init__(self, path: str) -> None:
        self.path = path
        try:
            tree = ET.parse(path)
        except ET.ParseError as exc:
            raise StackError(f"{os.path.basename(path)} is not valid XML: {exc}") from None
        self.root = tree.getroot()
        if self.root.tag != SVG + "svg":
            raise StackError(f"{os.path.basename(path)} is not an SVG")

        self.width = self._length("width")
        self.height = self._length("height")
        self.plot_x, self.plot_width = self._plot_area(path)
        self.shift = 0.0                      # set by align()

    def _length(self, attr: str) -> float:
        raw = self.root.get(attr) or ""
        match = _PT.search(raw)
        if not match:
            raise StackError(f"{os.path.basename(self.path)} has no {attr} in points")
        return float(match.group(1))

    def _plot_area(self, path: str) -> tuple[float, float]:
        """x and width of the largest <rect>: the depth axes' background.

        VirPlot's linear figure draws two axes backgrounds of equal width (the
        annotation panel and the depth panel), so the largest rect is the plot
        area whose left edge must line up across panels.
        """
        with open(path) as fp:
            rects = _RECT.findall(fp.read())
        if not rects:
            raise StackError(f"{os.path.basename(path)} has no plot area "
                             "(is it a VirPlot figure?)")
        x, _y, w, h = max(rects, key=lambda r: float(r[2]) * float(r[3]))
        return float(x), float(w)

    @property
    def name(self) -> str:
        return os.path.basename(self.path)


def order_panels(panels: list[Panel], by_length: bool = True) -> list[Panel]:
    """Largest first (the ICTV rule), or as given.

    Plot-area width is proportional to genome length within a run, so the
    widest plot area is the longest segment without re-reading the GFF.
    """
    if not by_length:
        return list(panels)
    return sorted(panels, key=lambda p: (-p.plot_width, p.name))


def align(panels: list[Panel], tolerance: float = 0.01) -> float:
    """Shift panels so every plot area starts at the same x.

    Returns the largest correction applied, in points. A tight bounding box
    crops each figure to its own content, so a panel whose depth axis carries
    wider tick labels starts further right; without this the segments would be
    offset from each other by that difference.
    """
    widest = max(p.plot_x for p in panels)
    worst = 0.0
    for panel in panels:
        panel.shift = widest - panel.plot_x
        worst = max(worst, panel.shift)
    return worst if worst > tolerance else 0.0


def stack(panels: list[Panel], out_path: str, gap: float = 4.0) -> None:
    """Write one SVG with each panel nested in a translated group."""
    width = max(p.width + p.shift for p in panels)
    height = sum(p.height for p in panels) + gap * (len(panels) - 1)

    svg = ET.Element(SVG + "svg", {
        "width": f"{width:.4f}pt", "height": f"{height:.4f}pt",
        "viewBox": f"0 0 {width:.4f} {height:.4f}", "version": "1.1",
    })
    y = 0.0
    for panel in panels:
        group = ET.SubElement(svg, SVG + "g", {
            "transform": f"translate({panel.shift:.4f},{y:.4f})"})
        # nest, rather than merge: each panel keeps its own coordinate system,
        # so nothing is rescaled and element ids cannot collide across panels
        inner = ET.SubElement(group, SVG + "svg", {
            "width": f"{panel.width:.4f}", "height": f"{panel.height:.4f}",
            "viewBox": panel.root.get("viewBox", f"0 0 {panel.width} {panel.height}"),
            "overflow": "visible",
        })
        for child in list(panel.root):
            inner.append(child)
        y += panel.height + gap

    ET.ElementTree(svg).write(out_path, encoding="utf-8", xml_declaration=True)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Stack VirPlot segment figures into one multipartite figure",
        epilog="Render the segments in one virplot run so they share a scale.")
    p.add_argument("-i", "--input", nargs="+", required=True, metavar="SVG",
                   help="Segment figures to stack (one virplot run)")
    p.add_argument("-o", "--output", default="./stacked.svg",
                   help="Output path [%(default)s]")
    p.add_argument("--gap", type=float, default=4.0,
                   help="Vertical gap between panels, in points [%(default)s]")
    p.add_argument("--order", choices=["length", "given"], default="length",
                   help="Row order: largest segment first, or as listed [%(default)s]")
    p.add_argument("--no-align", action="store_true",
                   help="Do not shift panels to align their plot areas")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if len(args.input) < 2:
        print("[stack] Error: give at least two figures to stack", file=sys.stderr)
        sys.exit(1)

    try:
        panels = [Panel(path) for path in args.input]
    except (StackError, OSError) as exc:
        print(f"[stack] Error: {exc}", file=sys.stderr)
        sys.exit(1)

    panels = order_panels(panels, by_length=args.order == "length")
    if not args.no_align:
        worst = align(panels)
        if worst:
            print(f"[stack] Aligned plot areas; largest shift {worst:.2f} pt")

    stack(panels, args.output, gap=args.gap)
    print(f"[stack] {len(panels)} panels -> {args.output}")
    for panel in panels:
        print(f"[stack]   {panel.name} (plot area {panel.plot_width:.1f} pt)")


if __name__ == "__main__":
    main()
