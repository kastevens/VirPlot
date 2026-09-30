"""Tests for the ICTV drawing conventions: layout rules, palette, settings, glyphs."""

import argparse
import os
import textwrap

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon, Rectangle

from virplot.layout import ABOVE, BELOW, place_features, resolve_mode, spans_overlap
from virplot.models import RNA, Feature
from virplot.parsers import parse_gff_rnas
from virplot.plotting import LinearPlotter
from virplot.settings import Settings, classify_function, lighter, load_settings


def F(start, end, strand="+", product="x", gene=None):
    return Feature(start, end, strand, product, gene)


def rna_with(*feats, length=10000, circular=False):
    return RNA(name="r", seqid="r", length=length, circular=circular, features=list(feats))


def layout(rna, mode):
    return {p.feature.product: (p.side, p.tier) for p in
            place_features(rna.features, rna.feature_spans, mode)}


# --- overlap and mode resolution -----------------------------------------

def test_spans_overlap_inclusive_and_multi_span():
    assert spans_overlap([(1, 10)], [(10, 20)])
    assert not spans_overlap([(1, 9)], [(10, 20)])
    assert spans_overlap([(90, 100), (1, 5)], [(3, 4)])         # wrapped feature


def test_resolve_mode_auto_follows_strands():
    assert resolve_mode("auto", two_strand=False) == "flip"
    assert resolve_mode("auto", two_strand=True) == "tier"
    assert resolve_mode("tier", two_strand=False) == "tier"


def test_two_strand_property():
    assert not rna_with(F(1, 10), F(20, 30)).two_strand
    assert rna_with(F(1, 10), F(20, 30, "-")).two_strand
    assert not rna_with(F(1, 10, "-"), F(20, 30, "-")).two_strand   # all minus: one strand


# --- L1: flip on overlap ----------------------------------------------------

def test_flip_first_orf_above_and_non_overlapping_stay():
    r = rna_with(F(1, 100, product="a"), F(200, 300, product="b"), F(400, 500, product="c"))
    assert layout(r, "flip") == {"a": (ABOVE, 0), "b": (ABOVE, 0), "c": (ABOVE, 0)}


def test_flip_alternates_through_an_overlapping_run():
    r = rna_with(F(1, 100, product="a"), F(90, 200, product="b"),
                 F(190, 300, product="c"), F(400, 500, product="d"))
    assert layout(r, "flip") == {"a": (ABOVE, 0), "b": (BELOW, 0),
                                 "c": (ABOVE, 0), "d": (ABOVE, 0)}


def test_flip_long_upstream_orf_does_not_cause_a_collision():
    # d overlaps c (below) so wants above, but a still occupies above there
    r = rna_with(F(1, 3000, product="a"), F(100, 500, product="b"),
                 F(600, 700, product="c"), F(650, 800, product="d"))
    out = layout(r, "flip")
    assert out["a"] == (ABOVE, 0) and out["b"] == (BELOW, 0) and out["c"] == (BELOW, 0)
    assert out["d"] == (ABOVE, 1)                # the flipped side, one tier out


def test_flip_sorts_by_start_first():
    r = rna_with(F(500, 600, product="late"), F(1, 100, product="early"))
    placed = place_features(r.features, r.feature_spans, "flip")
    assert [p.feature.product for p in placed] == ["early", "late"]


# --- L2: tier by strand -----------------------------------------------------

def test_tier_side_is_fixed_by_strand_and_overlap_tiers():
    r = rna_with(F(1, 100, product="p1"), F(50, 150, product="p2"),
                 F(120, 200, product="p3"), F(300, 400, "-", "m1"), F(350, 450, "-", "m2"))
    out = layout(r, "tier")
    assert out["p1"] == (ABOVE, 0) and out["p2"] == (ABOVE, 1) and out["p3"] == (ABOVE, 0)
    assert out["m1"] == (BELOW, 0) and out["m2"] == (BELOW, 1)


def test_tier_never_flips_across_the_line():
    r = rna_with(F(1, 100, product="p1"), F(50, 150, product="p2"), F(60, 160, product="p3"))
    out = layout(r, "tier")
    assert {s for s, _ in out.values()} == {ABOVE}
    assert sorted(t for _, t in out.values()) == [0, 1, 2]


# --- C1/C2: nest inward -----------------------------------------------------

def test_nest_puts_overlapping_arc_one_tier_in():
    r = rna_with(F(292, 807, product="V2"), F(710, 1384, product="V1"),
                 F(1365, 1736, product="V3"), F(2250, 3044, "-", "C1"),
                 F(2408, 2890, "-", "C3"), length=3206, circular=True)
    out = layout(r, "nest")
    # largest first: V1 and C1 take the outer lane, the ORFs overlapping them nest inside
    assert out["V1"][1] == 0 and out["V2"][1] == 1 and out["V3"][1] == 1
    assert out["C1"][1] == 0 and out["C3"][1] == 1            # AC4-inside-AC1 pattern


def test_nest_respects_origin_crossing_spans():
    r = rna_with(F(950, 50, product="wrap"), F(20, 40, product="inside"), length=1000, circular=True)
    out = layout(r, "nest")
    assert out["wrap"][1] == 0 and out["inside"][1] == 1


# --- palette ---------------------------------------------------------------

def test_classify_function_families():
    assert classify_function("RdRp") == "replicase"
    assert classify_function("replication-associated protein") == "replicase"
    assert classify_function("coat protein") == "capsid"
    assert classify_function("CPm") == "capsid"
    assert classify_function("movement protein") == "movement"
    assert classify_function("HSP70h") == "hsp70"
    assert classify_function("RNA silencing suppressor") == "suppressor"
    assert classify_function("p22") is None            # a bare pN name carries no function
    assert classify_function("V1 protein") is None


def test_feature_color_precedence_and_putative_tint():
    s = Settings(color_mapping={"RdRp": "#123456"})
    assert s.feature_color("RdRp") == "#123456"              # explicit wins
    assert s.feature_color("CP") == s.function_palette["capsid"]
    assert s.feature_color("putative coat protein") == lighter(s.function_palette["capsid"])
    assert s.feature_color("mystery") == s.default_color


def test_lighter_moves_towards_white():
    assert lighter("#000000", 0.5) == "#808080"
    assert lighter("#ffffff") == "#ffffff"


# --- settings keys ---------------------------------------------------------

def test_load_settings_new_keys(tmp_path):
    y = tmp_path / "s.yml"
    y.write_text(textwrap.dedent("""\
        overlap_mode: TIER
        end_5_label: "5' m7G"
        end_3_label: "A(n)"
        function_palette:
          capsid: '#000000'
    """))
    s = load_settings(str(y))
    assert s.overlap_mode == "tier"
    assert (s.end_5_label, s.end_3_label) == ("5' m7G", "A(n)")
    assert s.function_palette["capsid"] == "#000000"
    assert s.function_palette["replicase"]                    # defaults kept


def test_load_settings_bad_overlap_mode_falls_back(tmp_path):
    y = tmp_path / "s.yml"
    y.write_text("overlap_mode: sideways\n")
    assert load_settings(str(y)).overlap_mode == "auto"


# --- GFF gene attribute ----------------------------------------------------

def test_gene_attribute_parsed(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "A\t.\tregion\t1\t500\t.\t+\t.\tID=A\n"
        "A\t.\tCDS\t1\t300\t.\t+\t0\tID=c1;gene=ORF1a;product=RdRp\n"
        "A\t.\tCDS\t320\t400\t.\t-\t0\tID=c2;product=CP\n"
    )
    (rna,) = parse_gff_rnas(str(gff))
    assert rna.features[0].gene == "ORF1a" and rna.features[1].gene is None
    assert rna.two_strand


# --- linear glyphs -----------------------------------------------------------

def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=True, grid=False,
                no_label=False, no_border=False, yscale="linear", linthresh=10.0,
                format="png", outdir=".", name="x", layout="linear")
    base.update(over)
    return argparse.Namespace(**base)


def draw_linear(rna, settings=None):
    fig, ax = plt.subplots()
    p = LinearPlotter(settings or Settings(), make_args())
    p._draw_annotations(ax, rna)
    return fig, ax, p


def test_two_strand_genome_uses_arrows_one_strand_uses_boxes():
    fig, ax, _ = draw_linear(rna_with(F(1, 1000), F(2000, 3000, "-")))
    assert sum(isinstance(a, Polygon) for a in ax.patches) == 2
    plt.close(fig)
    fig, ax, _ = draw_linear(rna_with(F(1, 1000), F(2000, 3000)))
    assert sum(isinstance(a, Rectangle) for a in ax.patches) == 2
    assert not any(isinstance(a, Polygon) for a in ax.patches)
    plt.close(fig)


def test_arrow_points_head_at_reading_end():
    p = LinearPlotter(Settings(), make_args())
    fwd = p._arrow_points(100, 1000, 0.0, 1.0, True, 10000)
    rev = p._arrow_points(100, 1000, 0.0, 1.0, False, 10000)
    assert max(x for x, _ in fwd) == 1000 and fwd[2] == (1000, 0.5)   # tip at end
    assert min(x for x, _ in rev) == 100 and rev[2] == (100, 0.5)     # tip at start


def test_vpg_end_draws_grey_oval():
    fig, ax, _ = draw_linear(rna_with(F(1, 1000)), Settings(end_5_label="VPg"))
    assert any(isinstance(a, Ellipse) for a in ax.patches)
    plt.close(fig)


def test_default_title_is_name_and_length():
    fig, ax, p = draw_linear(rna_with(F(1, 1000), length=15468))
    artist = p._add_title(fig, rna_with(F(1, 1000), length=15468))
    assert artist.get_text() == "r (15,468 nts)"
    plt.close(fig)


def test_ylim_grows_with_tiers():
    fig, ax, _ = draw_linear(rna_with(F(1, 1000), F(50, 1200), F(60, 1300)))   # 3 tiers above
    assert ax.get_ylim()[1] > 2.0
    plt.close(fig)


# --- against the ICTV figure: BYV (Closteroviridae Fig. 2) -----------------

def test_byv_fixture_layout_versus_ictv_figure():
    """Real BYV coordinates under the L1 rule, compared with the published panel.

    The figure draws: ORF1a above, RdRp below, p6 below, Hsp70h above, p64
    below, CPm above, CP below, p20 above, p21 below. The flip rule reproduces
    every overlap-driven flip. It diverges exactly where the figure's choice is
    not overlap-driven: RdRp (a +1 frameshift that abuts ORF1a at 7997/7999
    without overlapping) and CP (70 nt clear of CPm, flipped for legibility).
    Those two need expression semantics / a judgement the data does not carry.
    """
    ex = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")
    (rna,) = parse_gff_rnas(os.path.join(ex, "byv.gff3"))
    assert not rna.two_strand
    out = layout(rna, resolve_mode("auto", rna.two_strand))
    side = {k: v[0] for k, v in out.items()}
    assert all(t == 0 for _, t in out.values())          # no tiering needed

    # overlap-driven flips agree with the figure (relative to each neighbour)
    assert side["Hsp70h"] == -side["p6"]                 # 9608 shared base
    assert side["p64"] == -side["Hsp70h"]
    assert side["CPm"] == -side["p64"]
    assert side["p20"] == -side["CP"]
    assert side["p21"] == -side["p20"]
    # non-overlap: rule keeps the side; figure agrees for RdRp/p6...
    assert side["p6"] == side["RdRp"]
    # ...and disagrees for the two documented cases
    assert side["RdRp"] == side["L-Pro/Mtr/Hel"]         # figure: flipped (+1 FS)
    assert side["CP"] == side["CPm"]                     # figure: flipped (no overlap)
