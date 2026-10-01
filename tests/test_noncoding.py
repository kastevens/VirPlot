"""Non-coding landmarks: UTR / intergenic-region / stem-loop rows, how they are
parsed from RefSeq and hand-written GFF3, and how both layouts draw them."""

import argparse
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

from virplot.models import RNA, Feature, Noncoding
from virplot.parsers import parse_gff_rnas
from virplot.plotting import CircularPlotter, LinearPlotter, Plotter
from virplot.settings import DEFAULT_FUNCTION_PALETTE, Settings, classify_function


def gff(tmp_path, body, length=3000, circular=True):
    p = tmp_path / "t.gff3"
    circ = "true" if circular else "false"
    p.write_text("##gff-version 3\n" + f"A\t.\tregion\t1\t{length}\t.\t+\t.\tID=A;Is_circular={circ}\n"
                 + textwrap.dedent(body))
    (rna,) = parse_gff_rnas(str(p))
    return rna


# --- parsing -------------------------------------------------------------------

def test_dedicated_types_are_always_read_with_default_labels(tmp_path):
    rna = gff(tmp_path, """\
        A	RefSeq	five_prime_UTR	1	100	.	+	.	ID=u5
        A	RefSeq	three_prime_UTR	2900	3000	.	+	.	ID=u3
        A	RefSeq	stem_loop	120	150	.	+	.	ID=sl;Note=conserved stem-loop structure
        A	RefSeq	origin_of_replication	2950	2960	.	+	.	ID=ori
    """, circular=False)
    by = {(n.kind, n.start): n for n in rna.noncoding}
    assert by[("region", 1)].label == "5′ UTR"
    assert by[("region", 2900)].label == "3′ UTR"
    assert by[("stem_loop", 120)].label == ""            # a Note describes, it does not name
    assert by[("region", 2950)].label == "ori"
    assert rna.features == []


def test_catch_all_rows_are_read_only_when_their_words_say_noncoding(tmp_path):
    rna = gff(tmp_path, """\
        A	RefSeq	CDS	327	1070	.	+	0	ID=cp;product=CP
        A	RefSeq	sequence_feature	2501	326	.	+	.	ID=a;Note=common region;gbkey=misc_feature
        A	RefSeq	sequence_feature	600	640	.	+	.	ID=b;Note=RdRp motif;gbkey=misc_feature
        A	RefSeq	regulatory_region	1500	1506	.	-	.	ID=c;regulatory_class=polyA_signal_sequence
        A	.	misc_feature	2000	2100	.	+	.	ID=d;Name=IR
        A	.	repeat_region	2200	2210	.	+	.	ID=e;Note=iteron
    """)
    assert [(n.start, n.label) for n in rna.noncoding] == [(2501, "common region"), (2000, "IR")]
    assert len(rna.features) == 1


def test_label_precedence_name_product_gene_then_note(tmp_path):
    rna = gff(tmp_path, """\
        A	.	misc_feature	10	20	.	+	.	Name=CRA;Note=common region
        A	.	misc_feature	30	40	.	+	.	product=LIR;Note=intergenic
        A	.	misc_feature	50	60	.	+	.	gene=SIR;Note=intergenic
        A	.	misc_feature	70	80	.	+	.	Note=intergenic region
    """)
    assert [n.label for n in rna.noncoding] == ["CRA", "LIR", "SIR", "intergenic region"]


def test_short_region_described_as_hairpin_becomes_a_stem_loop(tmp_path):
    rna = gff(tmp_path, """\
        A	.	misc_feature	100	130	.	+	.	Name=IR;Note=stem-loop with the nonanucleotide
        A	.	misc_feature	100	400	.	+	.	Name=LIR;Note=long intergenic region containing a stem-loop
    """)
    assert [n.kind for n in rna.noncoding] == ["stem_loop", "region"]


def test_origin_crossing_region_gets_two_spans(tmp_path):
    rna = gff(tmp_path, """\
        A	.	misc_feature	2801	3291	.	+	.	Name=IR
    """)
    (ir,) = rna.noncoding
    assert rna.feature_spans(ir) == [(2801, 3000), (1, 291)]


def test_classifier_knows_the_geminivirus_region_names():
    for name in ("IR", "LIR", "SIR", "CRA", "CRB", "common region", "5' UTR", "ori"):
        assert classify_function(name) == "noncoding", name
    assert classify_function("CP") == "capsid"


# --- drawing -----------------------------------------------------------------

def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, yscale="linear", linthresh=10.0,
                format="png", outdir=".", name="x", layout="linear")
    base.update(over)
    return argparse.Namespace(**base)


def begomo(circular=True):
    return RNA(name="r", seqid="r", length=2600, circular=circular,
               features=[Feature(327, 1070, "+", "CP", gene="AV1"),
                         Feature(1543, 2500, "-", "Rep", gene="AC1")],
               noncoding=[Noncoding(2501, 326, "+", "CRA", "region"),
                          Noncoding(110, 142, "+", "", "stem_loop")])


def texts(ax):
    return [t.get_text() for t in ax.texts]


GREY = DEFAULT_FUNCTION_PALETTE["noncoding"]


def test_linear_draws_grey_bars_on_the_line_and_a_hairpin():
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, begomo())
    bars = [p for p in ax.patches if to_hex(p.get_facecolor()) == GREY]
    assert len(bars) == 2                                     # origin-crossing: two pieces
    assert all(abs(p.get_y() + p.get_height() / 2 - 0.5) < 1e-9 for p in bars)   # astride y0
    markers = [l for l in ax.lines if l.get_marker() == "o"]
    assert len(markers) == 1                                  # the hairpin head
    plt.close(fig)


def test_region_holding_an_unnamed_stem_loop_is_named_once_at_the_hairpin():
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, begomo())
    cra = [t for t in ax.texts if t.get_text() == "CRA"]
    assert len(cra) == 1
    assert abs(cra[0].get_position()[0] - 126.5) < 1e-6         # over the hairpin, not the bar
    plt.close(fig)


def test_named_stem_loop_keeps_its_own_name_and_the_region_keeps_its_own():
    r = begomo()
    r.noncoding[1] = Noncoding(110, 142, "+", "nick", "stem_loop")
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, r)
    assert {"CRA", "nick"} <= set(texts(ax))
    plt.close(fig)


def test_linear_hairpin_drops_below_when_a_box_is_above():
    r = RNA(name="r", seqid="r", length=5000,
            features=[Feature(100, 2000, "+", "P")],
            noncoding=[Noncoding(1000, 1030, "+", "hp", "stem_loop")])
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, r)
    (head,) = [l for l in ax.lines if l.get_marker() == "o"]
    assert head.get_ydata()[0] < 0.5
    plt.close(fig)


def test_noncoding_takes_no_part_in_layout():
    r = begomo(circular=False)
    plotter = LinearPlotter(Settings(), make_args())
    placements, _ = plotter._placements(r, circular=False)
    assert len(placements) == 2


def test_noncoding_colour_is_grey_unless_mapped():
    p = Plotter(Settings(color_mapping={"CRA": "#123456"}), make_args())
    assert p._noncoding_color(Noncoding(1, 2, "+", "CRA")) == "#123456"
    assert p._noncoding_color(Noncoding(1, 2, "+", "IR")) == GREY


def test_text_on_dark_fill_is_white_only_for_really_dark_fills():
    assert Plotter._text_on("#f5e663") == "black"
    assert Plotter._text_on("#7b5ea7") == "black"             # ICTV purple keeps black text
    assert Plotter._text_on("#102030") == "white"
    assert Plotter._text_on(GREY, dark_below=0.5) == "white"


def test_circular_draws_ir_arc_astride_the_circle_and_moves_the_hairpin():
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="polar")
    CircularPlotter(Settings(), make_args(layout="circular"))._draw_annotations(ax, begomo())
    grey = [p for p in ax.patches if to_hex(p.get_facecolor()) == GREY]
    assert len(grey) == 2
    heads = [l for l in ax.lines if l.get_marker() == "o"]
    assert len(heads) == 1
    theta = heads[0].get_xdata()[0]
    assert abs(theta - 2 * 3.141592653589793 * 125.5 / 2600) < 1e-6   # at the stem-loop, not at 1
    assert "CRA" in texts(ax)
    plt.close(fig)


def test_circular_without_stem_loop_keeps_the_origin_icon():
    r = begomo(); r.noncoding = r.noncoding[:1]
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="polar")
    CircularPlotter(Settings(), make_args(layout="circular"))._draw_annotations(ax, r)
    (head,) = [l for l in ax.lines if l.get_marker() == "o"]
    assert head.get_xdata()[0] == 0.0
    plt.close(fig)


def test_circular_horizontal_mode_labels_region_level():
    r = begomo(); r.noncoding = r.noncoding[:1]
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="polar")
    s = Settings(circular_labels="horizontal")
    CircularPlotter(s, make_args(layout="circular"))._draw_annotations(ax, r)
    (cra,) = [t for t in ax.texts if t.get_text() == "CRA"]
    assert cra.get_rotation() == 0
    plt.close(fig)


def test_grbv_example_carries_its_intergenic_region():
    (rna,) = parse_gff_rnas("examples/grbv.gff3")
    (ir,) = rna.noncoding
    assert ir.label == "IR" and rna.feature_spans(ir) == [(3045, 3206), (1, 291)]
