"""Subgenomic RNA rows: the GFF3 curation convention that marks them, the
ladder's order, and how the linear panel draws them (the circular layout does
not — see docs/ictv_drawing_conventions.md §6 K)."""

import argparse
import logging
import textwrap

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

from virplot.models import RNA, Feature, SubgenomicRNA
from virplot.parsers import parse_gff_rnas
from virplot.plotting import CircularPlotter, LinearPlotter
from virplot.settings import Settings


def gff(tmp_path, body, length=3000, circular=False):
    p = tmp_path / "t.gff3"
    circ = "true" if circular else "false"
    p.write_text("##gff-version 3\n"
                 + f"A\t.\tregion\t1\t{length}\t.\t+\t.\tID=A;Is_circular={circ}\n"
                 + textwrap.dedent(body))
    (rna,) = parse_gff_rnas(str(p))
    return rna


def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, yscale="linear", linthresh=10.0,
                format="png", outdir=".", name="x", layout="linear")
    base.update(over)
    return argparse.Namespace(**base)


# --- the curation convention ---------------------------------------------------

def test_a_transcript_row_is_an_sgrna_only_when_it_says_so(tmp_path):
    """RefSeq writes real mRNA rows (spliced mastrevirus transcripts) that are
    not sgRNAs; the marker is required, never inferred from the type alone."""
    rna = gff(tmp_path, """\
        A	RefSeq	mRNA	100	3000	.	+	.	ID=m1;gene=V2
        A	VirPlot	mRNA	2000	3000	.	+	.	ID=m2;gene=CP;Note=sgRNA
    """)
    assert [s.label for s in rna.sgrnas] == ["CP"]


def test_the_marker_is_accepted_in_any_of_its_spellings(tmp_path):
    rna = gff(tmp_path, """\
        A	.	mRNA	500	3000	.	+	.	ID=a;gene=one;Note=sgRNA
        A	.	transcript	600	3000	.	+	.	ID=b;gene=two;Note=sg RNA for ORF2
        A	.	ncRNA	700	3000	.	+	.	ID=c;gene=three;Note=subgenomic RNA
        A	.	misc_RNA	800	3000	.	+	.	ID=d;gene=four;Note=sub-genomic
    """)
    assert [s.label for s in rna.sgrna_ladder] == ["one", "two", "three", "four"]


def test_an_sgrna_row_is_not_swallowed_by_the_noncoding_catch_all(tmp_path):
    """sequence_feature / misc_feature feed the non-coding reader too; a row
    marked as an sgRNA must become a ladder row, not a bar on the line."""
    rna = gff(tmp_path, """\
        A	.	sequence_feature	900	3000	.	+	.	ID=s;gene=CP;Note=sgRNA
    """)
    assert [s.label for s in rna.sgrnas] == ["CP"]
    assert rna.noncoding == []


def test_labels_strip_an_sgrna_prefix_and_fall_back_to_position(tmp_path):
    rna = gff(tmp_path, """\
        A	.	mRNA	900	3000	.	+	.	ID=a;gene=sgRNA CP;Note=sgRNA
        A	.	mRNA	1000	3000	.	+	.	ID=b;Name=subgenomic RNA-p20;Note=sgRNA
        A	.	mRNA	1100	3000	.	+	.	ID=c;Note=sgRNA
    """)
    assert [s.label for s in rna.sgrna_ladder] == ["CP", "p20", "sgRNA3"]


def test_an_sgrna_row_before_its_region_line_is_still_kept(tmp_path):
    """GFF3 has no ordering rule, so rows are collected first and attached
    afterwards — the same guarantee CDS and non-coding rows have."""
    p = tmp_path / "unordered.gff3"
    p.write_text("##gff-version 3\n"
                 "A\t.\tmRNA\t900\t3000\t.\t+\t.\tID=a;gene=CP;Note=sgRNA\n"
                 "A\t.\tregion\t1\t3000\t.\t+\t.\tID=A;Is_circular=false\n")
    (rna,) = parse_gff_rnas(str(p))
    assert [s.label for s in rna.sgrnas] == ["CP"]
    assert rna.sgrnas[0].end == 3000          # clamped against the length it learned later


def test_percent_encoded_attributes_are_decoded(tmp_path):
    """NCBI writes %2C / %3B in product= and Note=; the label must not keep them."""
    rna = gff(tmp_path, """\
        A	.	mRNA	900	3000	.	+	.	ID=a;gene=CP%2C major;Note=sgRNA
    """)
    assert [s.label for s in rna.sgrnas] == ["CP, major"]


def test_the_fallback_numbers_rows_five_to_three(tmp_path):
    """Written out of order, the sgRNA{n} fallback still counts along the genome."""
    rna = gff(tmp_path, """\
        A	.	mRNA	2000	3000	.	+	.	ID=late;Note=sgRNA
        A	.	mRNA	500	3000	.	+	.	ID=early;Note=sgRNA
    """)
    assert [(s.start, s.label) for s in rna.sgrna_ladder] == [(500, "sgRNA1"),
                                                             (2000, "sgRNA2")]


# --- coordinates ---------------------------------------------------------------

def test_the_three_prime_end_defaults_to_the_genome_end(tmp_path):
    """The common case is a 3'-coterminal set, often curated as a bare 5'
    coordinate; an end at or past the genome end means 'the genome's 3' end'."""
    rna = gff(tmp_path, """\
        A	.	mRNA	900	900	.	+	.	ID=a;gene=bare;Note=sgRNA
        A	.	mRNA	950	9999	.	+	.	ID=b;gene=past;Note=sgRNA
    """)
    assert {s.label: s.end for s in rna.sgrnas} == {"bare": 3000, "past": 3000}


def test_a_five_prime_proximal_sgrna_keeps_its_own_end(tmp_path):
    """Closteroviruses make 5'-terminal sgRNAs as well, so an earlier end that
    the curator gave is honoured rather than stretched to the 3' end."""
    rna = gff(tmp_path, """\
        A	.	mRNA	100	800	.	+	.	ID=a;gene=five_prime;Note=sgRNA
    """)
    (sg,) = rna.sgrnas
    assert (sg.start, sg.end) == (100, 800)


def test_an_sgrna_starting_past_the_genome_is_dropped(tmp_path, caplog):
    with caplog.at_level(logging.WARNING):
        rna = gff(tmp_path, """\
            A	.	mRNA	5000	5200	.	+	.	ID=a;gene=bad;Note=sgRNA
        """)
    assert rna.sgrnas == []
    assert "past the 3000 bp genome" in caplog.text


def test_the_ladder_is_ordered_longest_first():
    rna = RNA(name="r", seqid="r", length=1000, sgrnas=[
        SubgenomicRNA(800, 1000, "short"),
        SubgenomicRNA(200, 1000, "long"),
        SubgenomicRNA(500, 1000, "mid"),
    ])
    assert [s.label for s in rna.sgrna_ladder] == ["long", "mid", "short"]


# --- drawing -------------------------------------------------------------------

def closterovirus(n=3):
    """A small 3'-coterminal ladder under two ORFs."""
    return RNA(name="r", seqid="r", length=1000,
               features=[Feature(10, 500, "+", "RdRp", gene="ORF1"),
                         Feature(520, 900, "+", "CP", gene="ORF2")],
               sgrnas=[SubgenomicRNA(400 + i * 150, 1000, f"sg{i}") for i in range(n)])


def texts(ax):
    return [t.get_text() for t in ax.texts]


def test_linear_draws_one_row_per_sgrna_descending_from_the_orfs():
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, closterovirus())

    color = to_hex(Settings().sgrna_color)
    rows = [ln for ln in ax.lines
            if to_hex(ln.get_color()) == color and ln.get_linestyle() == "-"
            and len(ln.get_xdata()) == 2]
    assert len(rows) == 3
    ys = [ln.get_ydata()[0] for ln in rows]
    assert ys == sorted(ys, reverse=True)                 # each row below the last
    # every row ends at the 3' end and starts at its own 5' end
    assert [tuple(ln.get_xdata()) for ln in rows] == [(400, 1000), (550, 1000), (700, 1000)]
    plt.close(fig)


def test_the_ladder_hangs_below_the_lowest_orf_and_the_panel_grows_to_fit():
    fig, ax = plt.subplots(figsize=(12, 4))
    plotter = LinearPlotter(Settings(), make_args())
    plotter._draw_annotations(ax, closterovirus(n=6))
    bottom = ax.get_ylim()[0]

    fig2, ax2 = plt.subplots(figsize=(12, 4))
    bare = closterovirus(n=0)
    plotter._draw_annotations(ax2, bare)
    assert bottom < ax2.get_ylim()[0]                     # six rows pushed it down
    plt.close(fig)
    plt.close(fig2)


def test_each_row_is_labelled_and_no_label_silences_them():
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args())._draw_annotations(ax, closterovirus())
    assert {"sg0", "sg1", "sg2"} <= set(texts(ax))

    fig2, ax2 = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(), make_args(no_label=True))._draw_annotations(
        ax2, closterovirus())
    assert not {"sg0", "sg1", "sg2"} & set(texts(ax2))
    plt.close(fig)
    plt.close(fig2)


def test_a_genome_without_sgrnas_draws_no_rows_and_keeps_its_panel():
    fig, ax = plt.subplots(figsize=(12, 4))
    plotter = LinearPlotter(Settings(), make_args())
    top = plotter._draw_sgrnas(ax, closterovirus(n=0), top=-1.0)
    assert top == -1.0
    assert len(ax.lines) == 0
    plt.close(fig)


def test_the_row_colour_comes_from_the_yaml():
    fig, ax = plt.subplots(figsize=(12, 4))
    LinearPlotter(Settings(sgrna_color="#123456"),
                  make_args())._draw_annotations(ax, closterovirus(n=1))
    assert any(to_hex(ln.get_color()) == "#123456" for ln in ax.lines)
    plt.close(fig)


def test_circular_skips_the_ladder_and_says_why(tmp_path, caplog):
    """ICTV draws no transcript rows on circular genomes: geminivirus
    transcription is bidirectional from the IR, not a 3'-coterminal set."""
    rna = closterovirus()
    rna.circular = True
    rna.add_depth("s1", np.full(rna.length, 10))
    args = make_args(layout="circular", format="png", outdir=str(tmp_path))
    plotter = CircularPlotter(Settings(), args)
    plotter.prepare([rna])
    with caplog.at_level(logging.WARNING):
        plotter.render(rna, [], "c")
    assert "use --layout linear" in caplog.text


# --- the shipped fixture -------------------------------------------------------

def test_byv_example_carries_the_closterovirus_ladder():
    (rna,) = parse_gff_rnas("examples/byv.gff3")
    ladder = rna.sgrna_ladder
    assert [s.label for s in ladder] == ["p6", "Hsp70h", "p64", "CPm", "CP", "p20", "p21"]
    # the two 5' termini mapped in Vitushkina et al. 2002
    assert (ladder[0].start, ladder[1].start) == (9402, 9467)
    assert all(s.end == rna.length for s in ladder)       # 3'-coterminal
    # the ladder is annotation only: it adds no depth track
    assert rna.depth == []
