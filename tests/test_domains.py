"""Polyprotein domains: parsing mature-protein rows, clipping, and the dividers
and labels drawn for them in both layouts."""

import argparse
import logging
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

from virplot.models import RNA, Domain, Feature
from virplot.parsers import parse_gff_rnas
from virplot.plotting import DOMAIN_DIVIDER, CircularPlotter, LinearPlotter
from virplot.settings import Settings


def gff(tmp_path, body, length=10000, head="A\t.\tregion\t1\t{L}\t.\t+\t.\tID=A\n"):
    p = tmp_path / "t.gff3"
    p.write_text("##gff-version 3\n" + head.format(L=length) + textwrap.dedent(body))
    (rna,) = parse_gff_rnas(str(p))
    return rna


# --- parsing -------------------------------------------------------------------

def test_refseq_rows_attach_by_parent_in_genome_order(tmp_path):
    rna = gff(tmp_path, """\
        A	RefSeq	CDS	185	9376	.	+	0	ID=cds-X;product=polyprotein
        A	RefSeq	mature_protein_region_of_CDS	1037	2404	.	+	.	ID=id-X:285..740;Parent=cds-X;product=HC-Pro
        A	RefSeq	mature_protein_region_of_CDS	185	1036	.	+	.	ID=id-X:1..284;Parent=cds-X;product=P1
    """)
    (pp,) = rna.features
    assert [d.product for d in pp.domains] == ["P1", "HC-Pro"]
    assert pp.domains[0] == Domain(185, 1036, "P1")


def test_mat_peptide_spelling_and_containment_without_parent(tmp_path):
    rna = gff(tmp_path, """\
        A	.	CDS	100	5000	.	+	0	ID=pp;product=polyprotein
        A	.	CDS	6000	7000	.	+	0	ID=cp;product=CP
        A	.	mat_peptide	100	2000	.	+	.	product=Pro
        A	.	mat_peptide	2001	4997	.	+	.	product=RdRp
    """)
    pp, cp = rna.features
    assert [d.product for d in pp.domains] == ["Pro", "RdRp"]
    assert cp.domains == ()


def test_containment_prefers_the_smallest_cds_on_the_same_strand(tmp_path):
    rna = gff(tmp_path, """\
        A	.	CDS	100	9000	.	+	0	ID=big;product=polyprotein
        A	.	CDS	3000	4000	.	+	0	ID=small;product=inner
        A	.	CDS	3000	4000	.	-	0	ID=minus;product=other strand
        A	.	mat_peptide	3100	3900	.	+	.	product=dom
    """)
    by = {f.product: f for f in rna.features}
    assert [d.product for d in by["inner"].domains] == ["dom"]
    assert by["polyprotein"].domains == () and by["other strand"].domains == ()


def test_domain_matching_no_cds_is_warned_and_dropped(tmp_path):
    log = logging.getLogger("virplot.parsers")
    records = []
    handler = logging.Handler(); handler.emit = records.append
    log.addHandler(handler)
    try:
        rna = gff(tmp_path, """\
            A	.	CDS	100	1000	.	+	0	ID=a;product=P
            A	.	mat_peptide	5000	6000	.	+	.	Parent=nope;product=lost
        """)
    finally:
        log.removeHandler(handler)
    assert rna.features[0].domains == ()
    assert any("matches no CDS" in r.getMessage() for r in records)


def test_join_segments_take_the_domains_they_overlap_most(tmp_path):
    # coronavirus-style pp1ab: nsp10 in ORF1a, nsp13 in ORF1b, nsp12 straddling
    rna = gff(tmp_path, """\
        A	RefSeq	CDS	266	13468	.	+	0	ID=j;exception=ribosomal slippage;product=pp1ab
        A	RefSeq	CDS	13468	21555	.	+	0	ID=j;exception=ribosomal slippage;product=pp1ab
        A	RefSeq	mature_protein_region_of_CDS	12686	13024	.	+	.	Parent=j;product=nsp10
        A	RefSeq	mature_protein_region_of_CDS	13025	16236	.	+	.	Parent=j;product=nsp12
        A	RefSeq	mature_protein_region_of_CDS	16237	18039	.	+	.	Parent=j;product=nsp13
    """, length=30000)
    a, b = rna.features
    assert [d.product for d in a.domains] == ["nsp10"]
    assert [d.product for d in b.domains] == ["nsp12", "nsp13"]   # nsp12 mostly in ORF1b


# --- clipping ----------------------------------------------------------------

def poly():
    return Feature(185, 9376, "+", "polyprotein", domains=(
        Domain(185, 1036, "P1"), Domain(1037, 2404, "HC-Pro"), Domain(2405, 9373, "rest")))


def test_domains_within_clips_and_snaps_the_stop_codon():
    f = poly()
    assert f.domains_within(185, 9376)[-1] == Domain(2405, 9376, "rest")   # 9373 -> 9376
    assert f.domains_within(2000, 3000) == [Domain(2000, 2404, "HC-Pro"), Domain(2405, 3000, "rest")]
    assert f.domains_within(9380, 9400) == []


def test_dividers_are_between_bases_and_never_at_the_span_ends():
    f = poly()
    assert f.dividers_within(185, 9376) == [1036.5, 2404.5]
    assert f.dividers_within(2000, 3000) == [2404.5]
    # a gap before the first domain earns a line where the domain begins
    g = Feature(1, 5000, "+", "pp", domains=(Domain(1000, 4997, "only"),))
    assert g.dividers_within(1, 5000) == [999.5]


def test_feature_without_domains_has_no_dividers():
    assert Feature(1, 100, "+", "x").dividers_within(1, 100) == []


# --- drawing -----------------------------------------------------------------

def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, yscale="linear", linthresh=10.0,
                format="png", outdir=".", name="x", layout="linear")
    base.update(over)
    return argparse.Namespace(**base)


def rna_with_poly(circular=False):
    return RNA(name="r", seqid="r", length=10000, circular=circular, features=[poly()])


def texts(ax):
    return [t.get_text() for t in ax.texts]


def divider_lines(ax):
    return [l for l in ax.lines if l.get_linewidth() == DOMAIN_DIVIDER["linewidth"]
            and to_hex(l.get_color()) == "#000000" and len(l.get_xdata()) == 2
            and l.get_xdata()[0] == l.get_xdata()[1]]


def test_linear_draws_one_divider_per_boundary_and_labels_domains_not_the_polyprotein():
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, rna_with_poly())
    assert len(divider_lines(ax)) == 2
    t = texts(ax)
    assert {"P1", "HC-Pro", "rest"} <= set(t) and "polyprotein" not in t
    plt.close(fig)


def test_linear_domain_segments_take_their_own_colour():
    s = Settings(color_mapping={"P1": "#112233", "polyprotein": "#445566"})
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(s, make_args())._draw_annotations(ax, rna_with_poly())
    faces = {to_hex(p.get_facecolor()) for p in ax.patches}
    assert {"#112233", "#445566"} <= faces
    plt.close(fig)


def test_linear_narrow_domain_label_goes_outside_and_colliding_ones_stagger():
    # three adjacent 60-nt domains: none fits inside; their outside labels
    # would overprint, so they land on successive rows
    f = Feature(1, 10000, "+", "pp", domains=(
        Domain(5000, 5059, "AAAAAA"), Domain(5060, 5119, "BBBBBB"), Domain(5120, 5179, "CCCCCC")))
    r = RNA(name="r", seqid="r", length=10000, features=[f])
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, r)
    ys = {t.get_text(): t.get_position()[1] for t in ax.texts if t.get_text() in ("AAAAAA", "BBBBBB", "CCCCCC")}
    box_top = 0.5 + 0.6
    assert all(y > box_top for y in ys.values())            # all outside, above
    assert len(set(round(y, 3) for y in ys.values())) == 3  # three distinct rows
    plt.close(fig)


def test_no_label_suppresses_domain_text_but_keeps_dividers():
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args(no_label=True))._draw_annotations(ax, rna_with_poly())
    assert texts(ax) == [] or not ({"P1", "HC-Pro"} & set(texts(ax)))
    assert len(divider_lines(ax)) == 2
    plt.close(fig)


def test_circular_draws_radial_dividers_and_domain_labels():
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="polar")
    CircularPlotter(Settings(), make_args(layout="circular"))._draw_annotations(ax, rna_with_poly(circular=True))
    assert len(divider_lines(ax)) == 2
    t = texts(ax)
    assert {"P1", "HC-Pro", "rest"} <= set(t) and "polyprotein" not in t
    plt.close(fig)


def test_circular_horizontal_mode_labels_each_domain_level():
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="polar")
    s = Settings(circular_labels="horizontal")
    CircularPlotter(s, make_args(layout="circular"))._draw_annotations(ax, rna_with_poly(circular=True))
    labels = [a for a in ax.texts if a.get_text() in ("P1", "HC-Pro", "rest")]
    assert len(labels) == 3 and all(a.get_rotation() == 0 for a in labels)
    plt.close(fig)


def test_frameshift_label_carries_any_shift_size():
    r = RNA(name="r", length=10000, features=[
        Feature(100, 4000, "+", "P3"), Feature(2900, 3100, "+", "PIPO", mechanism="frameshift", shift=2)])
    fig, ax = plt.subplots()
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, r)
    assert "+2 FS" in texts(ax)
    plt.close(fig)
