"""--border / --no-border: feature glyphs are flat colour unless asked for."""

import argparse

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

from virplot.cli import build_parser
from virplot.models import RNA, Feature
from virplot.plotting import CircularPlotter, LinearPlotter
from virplot.settings import Settings


def parse(*extra):
    return build_parser().parse_args(["-g", "x.gff3", "-y", "x.yml", *extra])


def test_no_border_is_the_default():
    assert parse().border is False


def test_no_border_flag_still_accepted():
    assert parse("--no-border").border is False


def test_border_flag_turns_outlines_on():
    assert parse("--border").border is True


def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, yscale="linear", linthresh=10.0,
                format="png", outdir=".", name="x", layout="linear")
    base.update(over)
    return argparse.Namespace(**base)


def rna(circular=False):
    feats = [Feature(1, 1000, "+", "RdRp"), Feature(2000, 3000, "-", "CP")]
    return RNA(name="r", seqid="r", length=5000, circular=circular, features=feats)


def glyph_edges(plotter_cls, circular, **over):
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar" if circular else None)
    plotter_cls(Settings(), make_args(**over))._draw_annotations(ax, rna(circular))
    edges = [p.get_edgecolor() for p in ax.patches]
    plt.close(fig)
    return edges


@pytest.mark.parametrize("cls,circular", [(LinearPlotter, False), (CircularPlotter, True)])
def test_glyphs_have_no_outline_by_default(cls, circular):
    edges = glyph_edges(cls, circular)
    assert edges and all(e[3] == 0 for e in edges)          # fully transparent edge


@pytest.mark.parametrize("cls,circular", [(LinearPlotter, False), (CircularPlotter, True)])
def test_border_outlines_glyphs_in_black(cls, circular):
    edges = glyph_edges(cls, circular, border=True)
    assert edges and all(tuple(e) == to_rgba("black") for e in edges)


def test_args_without_border_attribute_mean_no_border():
    # Namespaces built by hand (as other tests do) predate the option.
    args = make_args()
    assert not hasattr(args, "border")
    assert LinearPlotter(Settings(), args)._glyph_edge() == {"edgecolor": "none"}
