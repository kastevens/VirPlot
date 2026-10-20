"""Composing segment figures into one stacked multipartite figure.

VirPlot renders one file per RNA instead of a stacked figure, and
`bin/stack_figures.py` joins them afterwards. That is only sound because the
separate panels already share an exact x-scale, so the first test here guards
that property: if it ever breaks, composing silently produces a figure whose
rows do not line up, and the design decision recorded in
docs/ictv_drawing_conventions.md §6 G stops being true.
"""

import argparse
import importlib.util
import re
import xml.etree.ElementTree as ET

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from virplot.cli import build_parser
from virplot.models import RNA, Feature
from virplot.plotting import LinearPlotter
from virplot.settings import Settings

_SPEC = importlib.util.spec_from_file_location(
    "stack_figures", __file__.rsplit("/tests/", 1)[0] + "/bin/stack_figures.py")
stack_figures = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(stack_figures)

SVG = "{http://www.w3.org/2000/svg}"
_PT = re.compile(r"([\d.]+)pt")
_RECT = re.compile(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" height="([\d.]+)"')


def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, yscale="linear", linthresh=10.0, bare_x=False,
                format="svg", outdir=".", name="x", layout="linear")
    base.update(over)
    return argparse.Namespace(**base)


def segment(seqid: str, length: int, depth: int = 10) -> RNA:
    rna = RNA(name=seqid, seqid=seqid, length=length,
              features=[Feature(10, max(20, length // 2), "+", "CP", gene="ORF1")])
    rna.add_depth("s1", np.full(length, depth))
    return rna


def render_all(tmp_path, rnas, **over) -> list[str]:
    """Render every RNA through one plotter, as a single virplot run does."""
    args = make_args(outdir=str(tmp_path), **over)
    plotter = LinearPlotter(Settings(), args)
    plotter.prepare(rnas)
    paths = []
    for rna in rnas:
        base = f"fig.{rna.seqid}"
        plotter.render(rna, [], base)
        paths.append(str(tmp_path / f"{base}.svg"))
    return paths


def plot_area(path: str) -> tuple[float, float]:
    """(x, width) of the largest rect — the depth axes' background."""
    rects = _RECT.findall(open(path).read())
    x, _y, w, _h = max(rects, key=lambda r: float(r[2]) * float(r[3]))
    return float(x), float(w)


# --- the property the whole decision rests on ---------------------------------

def test_separately_rendered_segments_share_an_exact_x_scale(tmp_path):
    """Figure width is proportional to length and margins are fractional, so
    points-per-base is identical across panels. This is what makes stacking a
    document operation rather than a re-plot."""
    lengths = {"A": 8897, "B": 4821, "C": 2916}
    rnas = [segment(k, v) for k, v in lengths.items()]
    paths = render_all(tmp_path, rnas)

    scales = []
    for path, (seqid, length) in zip(paths, lengths.items()):
        _x, width = plot_area(path)
        scales.append(width / length)

    assert max(scales) - min(scales) < 1e-6, f"x-scale drifted across panels: {scales}"


def test_segments_rendered_together_start_at_the_same_x(tmp_path):
    """Equal depth magnitudes give equal y-tick label widths, so the plot
    areas already line up; the composer only has to correct the cases where
    they do not."""
    rnas = [segment("A", 6000), segment("B", 3000)]
    xs = [plot_area(p)[0] for p in render_all(tmp_path, rnas)]
    assert max(xs) - min(xs) < 0.01


# --- the composer -------------------------------------------------------------

def test_panels_are_ordered_largest_first(tmp_path):
    rnas = [segment("small", 2000), segment("big", 8000), segment("mid", 5000)]
    panels = [stack_figures.Panel(p) for p in render_all(tmp_path, rnas)]
    ordered = stack_figures.order_panels(panels)
    assert [p.plot_width for p in ordered] == sorted(
        (p.plot_width for p in panels), reverse=True)


def test_order_given_keeps_the_caller_s_order(tmp_path):
    rnas = [segment("small", 2000), segment("big", 8000)]
    panels = [stack_figures.Panel(p) for p in render_all(tmp_path, rnas)]
    kept = stack_figures.order_panels(panels, by_length=False)
    assert [p.name for p in kept] == [p.name for p in panels]


def test_align_shifts_every_panel_onto_the_widest_left_edge(tmp_path):
    """A panel whose depth axis carries wider tick labels starts further
    right once the figure is cropped to a tight bounding box."""
    rnas = [segment("A", 6000, depth=9), segment("B", 6000, depth=150000)]
    panels = [stack_figures.Panel(p) for p in render_all(tmp_path, rnas, free_y=True)]
    assert panels[0].plot_x != panels[1].plot_x          # the case being corrected

    worst = stack_figures.align(panels)
    assert worst > 0
    aligned = [p.plot_x + p.shift for p in panels]
    assert max(aligned) - min(aligned) < 1e-9


def test_no_shift_is_reported_when_panels_already_line_up(tmp_path):
    panels = [stack_figures.Panel(p)
              for p in render_all(tmp_path, [segment("A", 6000), segment("B", 3000)])]
    assert stack_figures.align(panels) == 0.0


def test_stacked_output_nests_each_panel_and_sums_their_heights(tmp_path):
    rnas = [segment("A", 8000), segment("B", 4000), segment("C", 2000)]
    panels = [stack_figures.Panel(p) for p in render_all(tmp_path, rnas)]
    out = str(tmp_path / "stacked.svg")
    stack_figures.stack(stack_figures.order_panels(panels), out, gap=4.0)

    root = ET.parse(out).getroot()
    groups = root.findall(SVG + "g")
    assert len(groups) == 3
    assert all(len(g.findall(SVG + "svg")) == 1 for g in groups)   # panels nested, not merged

    height = float(_PT.search(root.get("height")).group(1))
    assert height == pytest.approx(sum(p.height for p in panels) + 8.0)
    width = float(_PT.search(root.get("width")).group(1))
    assert width == pytest.approx(max(p.width for p in panels))


def test_panels_are_stacked_top_to_bottom_in_order(tmp_path):
    rnas = [segment("A", 8000), segment("B", 4000)]
    panels = stack_figures.order_panels(
        [stack_figures.Panel(p) for p in render_all(tmp_path, rnas)])
    out = str(tmp_path / "stacked.svg")
    stack_figures.stack(panels, out, gap=0.0)

    offsets = [float(re.search(r"translate\([\d.-]+,([\d.-]+)\)", g.get("transform")).group(1))
               for g in ET.parse(out).getroot().findall(SVG + "g")]
    assert offsets == sorted(offsets)
    assert offsets[0] == 0.0
    assert offsets[1] == pytest.approx(panels[0].height)


def test_a_non_virplot_svg_is_rejected(tmp_path):
    bad = tmp_path / "bad.svg"
    bad.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="10pt" height="10pt"></svg>')
    with pytest.raises(stack_figures.StackError, match="no plot area"):
        stack_figures.Panel(str(bad))


def test_a_file_that_is_not_xml_is_rejected(tmp_path):
    bad = tmp_path / "bad.svg"
    bad.write_text("not xml at all")
    with pytest.raises(stack_figures.StackError, match="not valid XML"):
        stack_figures.Panel(str(bad))


# --- the CLI flags that serve composing ---------------------------------------

def test_rna_order_defaults_to_largest_first():
    assert build_parser().parse_args(["-g", "g", "-d", "d", "-y", "y"]).rna_order == "length"


def test_bare_x_drops_the_axis_from_every_panel_but_the_bottom(tmp_path):
    rnas = [segment("A", 8000), segment("B", 4000), segment("C", 2000)]
    args = make_args(outdir=str(tmp_path), bare_x=True)
    plotter = LinearPlotter(Settings(), args)
    plotter.prepare(rnas)
    assert plotter.bottom_seqid == "C"
    assert plotter._bare_x(rnas[0]) and plotter._bare_x(rnas[1])
    assert not plotter._bare_x(rnas[2])              # the bottom row keeps it


def test_a_lone_panel_is_bare_because_something_else_carries_the_axis(tmp_path):
    rna = segment("A", 8000)
    plotter = LinearPlotter(Settings(), make_args(outdir=str(tmp_path), bare_x=True))
    plotter.prepare([rna])
    assert plotter._bare_x(rna)


def test_without_bare_x_every_panel_keeps_its_axis(tmp_path):
    rnas = [segment("A", 8000), segment("B", 4000)]
    plotter = LinearPlotter(Settings(), make_args(outdir=str(tmp_path)))
    plotter.prepare(rnas)
    assert not any(plotter._bare_x(r) for r in rnas)


def test_bare_x_keeps_the_ticks_so_panels_still_align(tmp_path):
    """Only the label and the numbers go; dropping the ticks would change the
    plot area and break the alignment the composer relies on."""
    rnas = [segment("A", 6000), segment("B", 3000)]
    plain = [plot_area(p) for p in render_all(tmp_path / "plain", rnas)]
    rnas = [segment("A", 6000), segment("B", 3000)]
    bare = [plot_area(p) for p in render_all(tmp_path / "bare", rnas, bare_x=True)]
    assert [w for _x, w in plain] == [pytest.approx(w) for _x, w in bare]
