"""Frameshift and readthrough: parsing from RefSeq attributes and Note=, layout, marks."""

import argparse
import textwrap

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from virplot.layout import ABOVE, BELOW, place_features
from virplot.models import RNA, Feature
from virplot.parsers import parse_gff_rnas
from virplot.plotting import LinearPlotter
from virplot.settings import Settings


def gff(tmp_path, body, length=10000):
    p = tmp_path / "t.gff3"
    p.write_text("##gff-version 3\n" + f"A\t.\tregion\t1\t{length}\t.\t+\t.\tID=A\n"
                 + textwrap.dedent(body))
    (rna,) = parse_gff_rnas(str(p))
    return rna


# --- RefSeq join with ribosomal slippage --------------------------------------

def test_refseq_join_plus_one_frameshift(tmp_path):
    rna = gff(tmp_path, """\
        A	RefSeq	CDS	108	7997	.	+	0	ID=cds-X;exception=ribosomal slippage;product=fusion protein
        A	RefSeq	CDS	7999	9393	.	+	0	ID=cds-X;exception=ribosomal slippage;product=fusion protein
    """)
    a, b = rna.features
    assert (a.mechanism, b.mechanism) == (None, "frameshift")
    assert b.shift == 1                          # skips nt 7998
    assert a.show_label and not b.show_label     # same product: label it once


def test_refseq_join_minus_one_frameshift(tmp_path):
    # coronavirus-style: second segment re-reads the last base of the first
    rna = gff(tmp_path, """\
        A	RefSeq	CDS	266	13468	.	+	0	ID=orf1ab;exception=ribosomal slippage;product=ORF1ab
        A	RefSeq	CDS	13468	21555	.	+	0	ID=orf1ab;exception=ribosomal slippage;product=ORF1ab
    """, length=30000)
    assert rna.features[1].mechanism == "frameshift" and rna.features[1].shift == -1


def test_refseq_join_on_minus_strand_orders_by_translation(tmp_path):
    rna = gff(tmp_path, """\
        A	RefSeq	CDS	1000	5000	.	-	0	ID=j;exception=ribosomal slippage;product=pp
        A	RefSeq	CDS	5002	9000	.	-	0	ID=j;exception=ribosomal slippage;product=pp
    """)
    first, cont = rna.features
    assert first.start == 5002 and cont.start == 1000       # 5' segment first on -
    assert cont.mechanism == "frameshift" and cont.shift == 1  # nt 5001 skipped


def test_join_without_slippage_note_stays_separate_boxes(tmp_path):
    rna = gff(tmp_path, """\
        A	.	CDS	100	500	.	+	0	ID=spliced;product=Rep
        A	.	CDS	800	1200	.	+	0	ID=spliced;product=Rep
    """)
    assert [f.mechanism for f in rna.features] == [None, None]


# --- Note= convention ---------------------------------------------------------

def test_note_frameshift_with_sign(tmp_path):
    rna = gff(tmp_path, """\
        A	.	CDS	100	4000	.	+	0	ID=a;gene=ORF1a;product=L-Pro/Mtr/Hel
        A	.	CDS	4002	6000	.	+	0	ID=b;gene=ORF1b;product=RdRp;Note=+1 frameshift from ORF1a
    """)
    assert rna.features[1].mechanism == "frameshift" and rna.features[1].shift == 1
    assert rna.features[1].show_label                        # its own product


def test_note_frameshift_without_sign_infers_from_neighbour(tmp_path):
    rna = gff(tmp_path, """\
        A	.	CDS	100	4000	.	+	0	ID=a;product=P1a
        A	.	CDS	4000	6000	.	+	0	ID=b;product=RdRp;Note=ribosomal frameshift
    """)
    assert rna.features[1].shift == -1


# --- readthrough ---------------------------------------------------------------

def test_transl_except_readthrough_is_trimmed_to_the_extension(tmp_path):
    # TMV-style: 126K at 69..3419, 183K at 69..4919 reading through the 126K stop
    rna = gff(tmp_path, """\
        A	RefSeq	CDS	69	3419	.	+	0	ID=p126;product=126 kDa replicase
        A	RefSeq	CDS	69	4919	.	+	0	ID=p183;product=183 kDa RdRp;transl_except=(pos:3417..3419%2Caa:OTHER)
    """)
    short, rt = rna.features
    assert short.mechanism is None and (short.start, short.end) == (69, 3419)
    assert rt.mechanism == "readthrough" and (rt.start, rt.end) == (3420, 4919)


def test_note_readthrough_kept_as_given_when_no_partner_shares_its_start(tmp_path):
    rna = gff(tmp_path, """\
        A	.	CDS	100	1000	.	+	0	ID=orf3;product=CP
        A	.	CDS	1001	1800	.	+	0	ID=orf5;product=RTD;Note=readthrough of ORF3 stop
    """)
    rt = rna.features[1]
    assert rt.mechanism == "readthrough" and (rt.start, rt.end) == (1001, 1800)


# --- layout -----------------------------------------------------------------

def spans_of(rna):
    return rna.feature_spans


def test_frameshift_flips_without_overlap():
    r = RNA(name="r", length=10000, features=[
        Feature(100, 4000, "+", "P1a"), Feature(4002, 6000, "+", "RdRp", mechanism="frameshift", shift=1)])
    out = {p.feature.product: p.side for p in place_features(r.features, r.feature_spans, "flip")}
    assert out["RdRp"] == -out["P1a"]


def test_readthrough_stays_on_partner_side_even_when_same_colour_and_close():
    r = RNA(name="r", length=10000, features=[
        Feature(100, 1000, "+", "CP"), Feature(1001, 1800, "+", "RTD", mechanism="readthrough")])
    always_flip = lambda prev, feat: True                      # e.g. same colour, adjacent
    out = {p.feature.product: p.side for p in
           place_features(r.features, r.feature_spans, "flip", flip_if=always_flip)}
    assert out["RTD"] == out["CP"]


# --- marks -----------------------------------------------------------------

def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, yscale="linear", linthresh=10.0,
                format="png", outdir=".", name="x", layout="linear")
    base.update(over)
    return argparse.Namespace(**base)


def texts(ax):
    return [t.get_text() for t in ax.texts]


def test_frameshift_label_drawn():
    r = RNA(name="r", length=10000, features=[
        Feature(100, 4000, "+", "P1a"), Feature(4002, 6000, "+", "RdRp", mechanism="frameshift", shift=1)])
    fig, ax = plt.subplots()
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, r)
    assert "+1 FS" in texts(ax)
    plt.close(fig)


def test_minus_one_uses_a_real_minus_sign():
    r = RNA(name="r", length=10000, features=[
        Feature(100, 4000, "+", "P1a"), Feature(4000, 6000, "+", "RdRp", mechanism="frameshift", shift=-1)])
    fig, ax = plt.subplots()
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, r)
    assert "−1 FS" in texts(ax)
    plt.close(fig)


def test_readthrough_bar_and_label_drawn():
    r = RNA(name="r", length=10000, features=[
        Feature(100, 1000, "+", "CP"), Feature(1001, 1800, "+", "RTD", mechanism="readthrough")])
    fig, ax = plt.subplots()
    before = len(ax.lines)
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, r)
    assert "RT" in texts(ax)
    bars = [l for l in ax.lines[before:] if len(l.get_xdata()) == 2 and l.get_xdata()[0] == l.get_xdata()[1] == 1001]
    assert len(bars) == 1
    plt.close(fig)


def test_join_segment_with_repeated_product_is_labelled_once():
    r = RNA(name="r", length=10000, features=[
        Feature(100, 4000, "+", "fusion"), Feature(4002, 6000, "+", "fusion", mechanism="frameshift",
                                                    shift=1, show_label=False)])
    fig, ax = plt.subplots()
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, r)
    assert texts(ax).count("fusion") == 1
    plt.close(fig)
