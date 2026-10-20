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


def test_feature_spans_wrap_spellings_are_split_even_on_a_linear_rna():
    # --topology linear on a circular genome, or a GFF that forgot Is_circular:
    # the pieces are drawn where the coordinates say rather than off the axis
    rna = RNA(name="r", length=100, circular=False)
    assert rna.feature_spans(Feature(90, 110, "+", "p")) == [(90, 100), (1, 10)]
    assert rna.feature_spans(Feature(90, 10, "+", "p")) == [(90, 100), (1, 10)]
    assert rna.feature_spans(Feature(10, 50, "+", "p")) == [(10, 50)]


def test_feature_spans_reduces_padded_start_modulo_length():
    rna = RNA(name="r", length=1000, circular=True)
    assert rna.feature_spans(Feature(1050, 1200, "+", "p")) == [(50, 200)]


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


def test_layout_defaults_to_the_shape_of_the_molecule(tmp_path):
    """A GFF that says Is_circular=true draws as a circle without being asked;
    a linear one stays a track."""
    used = record_layouts(lambda: run(tmp_path, "grbv.gff3", "grbv.sam", "grbv.yml"))
    assert used == ["circular"]

    used = record_layouts(lambda: run(tmp_path, "sample.gff3", "sample.dep", "spec.yml"))
    assert used == ["linear"]


def test_layout_linear_forces_a_track_for_a_circular_genome(tmp_path):
    """The escape hatch. Editing Is_circular out of the GFF would also turn off
    depth wrapping and origin-crossing features, so presentation needs its own
    control — see docs/ictv_drawing_conventions.md §6 H."""
    used = record_layouts(
        lambda: run(tmp_path, "grbv.gff3", "grbv.sam", "grbv.yml", "--layout", "linear"))
    assert used == ["linear"]


def test_spec_yml_can_set_the_layout_and_the_cli_overrides_it(tmp_path):
    ex = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")
    yml = tmp_path / "pinned.yml"
    yml.write_text(open(os.path.join(ex, "grbv.yml")).read() + "\nlayout: linear\n")

    used = record_layouts(lambda: run(tmp_path, "grbv.gff3", "grbv.sam", str(yml)))
    assert used == ["linear"]                       # the YAML pinned it

    used = record_layouts(
        lambda: run(tmp_path, "grbv.gff3", "grbv.sam", str(yml), "--layout", "circular"))
    assert used == ["circular"]                     # the CLI wins


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


def test_circular_label_helpers_work_before_annotations_are_drawn():
    # per-figure radii start at sane defaults rather than AttributeError
    import argparse
    from virplot.plotting import CircularPlotter
    from virplot.settings import Settings
    args = argparse.Namespace(smooth=False, normalize=False, free_y=False, equal_width=False,
                              legend=False, shade_breaks=False, title=False, grid=False,
                              no_label=False, yscale="linear", linthresh=10.0,
                              format="png", outdir=".", name="x", layout="circular")
    fig = plt.figure(); ax = fig.add_subplot(projection="polar")
    rna = RNA(name="r", length=1000, circular=True)
    CircularPlotter(Settings(), args)._label_horizontal(ax, rna, "x", [(1, 100)])
    plt.close(fig)


def test_tiny_genome_tick_step_is_at_least_one():
    from virplot.plotting import _nice_step
    assert _nice_step(4 / 8) == 1 and _nice_step(0) == 1


def test_circular_warns_that_frameshift_marks_are_linear_only(tmp_path, caplog):
    """GRBV has no frameshift; BYV does. Under the new default a circular
    genome reaches this layout unasked, so the dropped mark must be said."""
    import logging
    with caplog.at_level(logging.WARNING):
        run(tmp_path, "byv.gff3", "byv.sam", "byv.yml", "--layout", "circular")
    assert "not marked" in caplog.text and "--layout linear" in caplog.text


def test_circular_is_quiet_when_nothing_is_dropped(tmp_path, caplog):
    import logging
    with caplog.at_level(logging.WARNING):
        run(tmp_path, "grbv.gff3", "grbv.sam", "grbv.yml", "--layout", "circular")
    assert "not marked" not in caplog.text
