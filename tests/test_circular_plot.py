"""Tests for CircularPlotter and the --layout selection."""

import argparse
import os

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import virplot.cli as cli_mod
from virplot.cli import main
from virplot.models import RNA, Feature
from virplot.plotting import (
    CIRC_R_DEPTH_BASE,
    CIRC_R_DEPTH_MAX,
    CIRC_RMAX,
    CircularPlotter,
    LinearPlotter,
)
from virplot.settings import Settings


def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, yscale="linear", linthresh=10.0,
                format="png", outdir=".", name="x", layout="circular")
    base.update(over)
    return argparse.Namespace(**base)


def plotter(**over):
    return CircularPlotter(Settings(), make_args(**over))


# --- feature_spans (origin crossing) --------------------------------------

def test_feature_spans_simple():
    rna = RNA(name="r", length=100, circular=True)
    assert rna.feature_spans(Feature(10, 40, "+", "p")) == [(10, 40)]


def test_feature_spans_end_past_length_splits():
    rna = RNA(name="r", length=100, circular=True)
    assert rna.feature_spans(Feature(90, 110, "+", "p")) == [(90, 100), (1, 10)]


def test_feature_spans_start_after_end_splits():
    rna = RNA(name="r", length=100, circular=True)
    assert rna.feature_spans(Feature(90, 10, "+", "p")) == [(90, 100), (1, 10)]


def test_feature_spans_linear_rna_is_never_split():
    rna = RNA(name="r", length=100, circular=False)
    assert rna.feature_spans(Feature(90, 110, "+", "p")) == [(90, 110)]


def test_feature_spans_feature_wrapping_right_round_covers_once():
    rna = RNA(name="r", length=100, circular=True)
    assert rna.feature_spans(Feature(10, 130, "+", "p")) == [(1, 100)]


# --- geometry --------------------------------------------------------------

def test_theta_maps_first_base_to_zero_and_wraps_once():
    assert CircularPlotter._theta(1, 100) == 0.0
    assert CircularPlotter._theta(101, 100) == pytest.approx(2 * np.pi)
    assert CircularPlotter._theta(51, 100) == pytest.approx(np.pi)


def test_depth_radius_spans_the_band_and_clips():
    p = plotter()
    r = p._depth_radius(np.array([0.0, 30.0, 60.0, 90.0]), 60.0)
    assert r[0] == pytest.approx(CIRC_R_DEPTH_BASE)
    assert r[2] == pytest.approx(CIRC_R_DEPTH_MAX)
    assert r[1] == pytest.approx((CIRC_R_DEPTH_BASE + CIRC_R_DEPTH_MAX) / 2)
    assert r[3] == pytest.approx(CIRC_R_DEPTH_MAX)      # beyond ymax is clipped


def test_close_repeats_first_sample_at_two_pi():
    theta = np.array([0.0, np.pi])
    th, vals = CircularPlotter._close(theta, np.array([5.0, 9.0]))
    assert th[-1] == pytest.approx(2 * np.pi)
    assert vals[-1] == 5.0 and len(th) == len(vals) == 3


# --- drawing ---------------------------------------------------------------

def circular_rna(length=300, features=(), depth=None):
    rna = RNA(name="c", seqid="c", length=length, circular=True,
              features=[Feature(*f) for f in features])   # (start, end, strand, product[, gene])
    rna.add_depth("s", np.full(length, 10) if depth is None else depth)
    return rna


def test_origin_crossing_feature_draws_two_arcs():
    rna = circular_rna(features=[(280, 20, "+", "wrapped")])
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar")
    p = plotter()
    p._draw_annotations(ax, rna)
    # one closed path per arc, so a feature through the origin contributes two
    from matplotlib.patches import Polygon
    assert sum(isinstance(a, Polygon) for a in ax.patches) == 2
    plt.close(fig)


def test_radial_limits_are_pinned(tmp_path):
    """Autoscale would pad below the smallest radius and expose an inner spine."""
    rna = circular_rna()
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar")
    p = plotter()
    ax.set_rorigin(0)
    ax.set_ylim(0, CIRC_RMAX)
    p._draw_depth(ax, rna, rna.tracks, rna.total_depth())
    assert ax.get_ylim() == (0, CIRC_RMAX)
    plt.close(fig)


def test_depth_ymax_uses_shared_limit_when_prepared():
    a = circular_rna(depth=np.full(300, 40))
    b = circular_rna(depth=np.full(300, 10))
    p = plotter()
    p.prepare([a, b])
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar")
    p._draw_depth(ax, b, b.tracks, b.total_depth())
    assert p._depth_ymax == pytest.approx(40 * 1.05)   # a's peak, not b's
    plt.close(fig)


# --- CLI layout selection --------------------------------------------------

def examples_dir():
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def run(tmp_path, gff, depth, yaml, *extra, name="out"):
    ex = examples_dir()
    main(["-g", os.path.join(ex, gff), "-d", os.path.join(ex, depth),
          "-y", os.path.join(ex, yaml), "-o", str(tmp_path),
          "--name", name, "-f", "svg", *extra])


def record_layouts(fn):
    """Run fn with both plotters stubbed; return the layouts actually used."""
    used = []
    originals = (CircularPlotter.render, LinearPlotter.render)
    CircularPlotter.render = lambda self, *a, **k: used.append("circular")
    LinearPlotter.render = lambda self, *a, **k: used.append("linear")
    try:
        fn()
    finally:
        CircularPlotter.render, LinearPlotter.render = originals
    return used


def test_layout_circular_renders_a_file(tmp_path):
    run(tmp_path, "grbv.gff3", "grbv.sam", "grbv.yml", "--layout", "circular")
    assert os.listdir(tmp_path) == ["out.svg"]


def test_layout_defaults_to_linear_even_for_circular_genome(tmp_path):
    used = record_layouts(lambda: run(tmp_path, "grbv.gff3", "grbv.sam", "grbv.yml"))
    assert used == ["linear"]


def test_layout_auto_follows_topology(tmp_path):
    used = record_layouts(
        lambda: run(tmp_path, "grbv.gff3", "grbv.sam", "grbv.yml", "--layout", "auto"))
    assert used == ["circular"]

    used = record_layouts(
        lambda: run(tmp_path, "sample.gff3", "sample.dep", "spec.yml", "--layout", "auto"))
    assert used == ["linear"]


def test_layout_circular_works_on_a_linear_genome(tmp_path):
    run(tmp_path, "sample.gff3", "sample.dep", "spec.yml", "--layout", "circular")
    assert os.listdir(tmp_path) == ["out.svg"]


# --- ICTV-style options: horizontal labels, arcs on the circle -------------

def test_settings_circular_style_keys(tmp_path):
    from virplot.settings import load_settings
    y = tmp_path / "s.yml"
    y.write_text("circular_labels: HORIZONTAL\ncircular_arcs: on_circle\n")
    s = load_settings(str(y))
    assert (s.circular_labels, s.circular_arcs) == ("horizontal", "on_circle")
    y.write_text("circular_labels: sideways\ncircular_arcs: inside\n")
    s = load_settings(str(y))
    assert (s.circular_labels, s.circular_arcs) == ("tangential", "outside")   # defaults on bad values


def test_horizontal_labels_are_level_and_name_the_orf():
    rna = circular_rna(features=[(10, 120, "+", "CP", "V1"), (150, 250, "-", "Rep", "C1")])
    fig = plt.figure(); ax = fig.add_subplot(projection="polar")
    p = CircularPlotter(Settings(circular_labels="horizontal"), make_args())
    p._draw_annotations(ax, rna)
    labels = [t for t in ax.texts if "(" in t.get_text()]
    assert sorted(t.get_text() for t in labels) == ["C1 (Rep)", "V1 (CP)"]
    assert all(t.get_rotation() == 0 for t in labels)
    plt.close(fig)


def test_nested_arc_label_is_offset_not_overprinted():
    # C3 nests inside C1 and shares its mid-angle
    rna = circular_rna(features=[(2250, 3044, "-", "Rep", "C1"), (2408, 2890, "-", "REn", "C3")], length=3206)
    fig = plt.figure(); ax = fig.add_subplot(projection="polar")
    p = CircularPlotter(Settings(circular_labels="horizontal"), make_args())
    p._draw_annotations(ax, rna)
    offsets = {t.get_text(): t.xyann for t in ax.texts if "(" in t.get_text()}
    assert offsets["C1 (Rep)"][1] == 0 and offsets["C3 (REn)"][1] != 0
    plt.close(fig)


def test_on_circle_puts_innermost_lane_astride_the_baseline():
    from virplot.plotting import CIRC_R_BASELINE, CIRC_ANN_HEIGHT
    rna = circular_rna(features=[(10, 120, "+", "CP")])
    for mode, expect_top in (("outside", CIRC_R_BASELINE + 0.02 + CIRC_ANN_HEIGHT),
                             ("on_circle", CIRC_R_BASELINE + CIRC_ANN_HEIGHT / 2)):
        fig = plt.figure(); ax = fig.add_subplot(projection="polar")
        p = CircularPlotter(Settings(circular_arcs=mode), make_args())
        p._draw_annotations(ax, rna)
        assert p._ring_top == pytest.approx(expect_top), mode
        plt.close(fig)
