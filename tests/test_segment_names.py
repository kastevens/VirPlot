"""Segment names: a region's Name= titles the figure and labels the circular hub.

The sequence id stays the key for depth matching and output file names; only
what the reader sees changes (RNA-1 rather than NC_003617.1).
"""

import argparse
import logging
import os
import textwrap

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from virplot.cli import main
from virplot.models import RNA
from virplot.parsers import parse_gff_rnas
from virplot.plotting import CircularPlotter, LinearPlotter
from virplot.settings import Settings

NAMED_GFF = textwrap.dedent("""\
    ##gff-version 3
    NC_003617.1\tRefSeq\tregion\t1\t8118\t.\t+\t.\tID=NC_003617.1;Name=RNA-1;Is_circular=false
    NC_003618.1\tRefSeq\tregion\t1\t7193\t.\t+\t.\tID=NC_003618.1;Is_circular=false
    SEG3\t.\tregion\t1\t900\t.\t+\t.\tID=SEG3;Name=RNA%203
""")


def examples_dir():
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=True, grid=False,
                no_label=False, yscale="linear", linthresh=10.0,
                format="svg", outdir=".", name="x", layout="auto")
    base.update(over)
    return argparse.Namespace(**base)


# --- parsing -----------------------------------------------------------------

@pytest.fixture
def rnas(tmp_path):
    gff = tmp_path / "named.gff3"
    gff.write_text(NAMED_GFF)
    return parse_gff_rnas(str(gff))


def test_region_name_is_the_display_name(rnas):
    assert rnas[0].name == "RNA-1"
    assert rnas[0].seqid == "NC_003617.1"        # the join key is untouched


def test_region_without_name_falls_back_to_the_seqid(rnas):
    assert rnas[1].name == rnas[1].seqid == "NC_003618.1"


def test_region_name_is_percent_decoded(rnas):
    assert rnas[2].name == "RNA 3"


def test_shipped_segmented_examples_carry_their_segment_names():
    ex = examples_dir()
    names = {os.path.basename(g): [r.name for r in parse_gff_rnas(os.path.join(ex, g))]
             for g in ("liyv.gff3", "tswv.gff3", "tgmv.gff3")}
    assert names == {"liyv.gff3": ["RNA-1", "RNA-2"],
                     "tswv.gff3": ["L", "M", "S"],
                     "tgmv.gff3": ["DNA-A", "DNA-B"]}


# --- drawing -----------------------------------------------------------------

def test_linear_title_is_the_segment_name_and_length():
    rna = RNA(name="RNA-1", seqid="NC_003617.1", length=8118)
    fig = plt.figure()
    title = LinearPlotter(Settings(), args())._add_title(fig, rna)
    assert title.get_text() == "RNA-1 (8,118 nts)"
    plt.close(fig)


def test_yaml_title_still_wins_over_the_segment_name():
    rna = RNA(name="RNA-1", seqid="NC_003617.1", length=8118)
    fig = plt.figure()
    title = LinearPlotter(Settings(title="LIYV"), args())._add_title(fig, rna)
    assert title.get_text() == "LIYV"
    plt.close(fig)


def test_circular_hub_names_the_component():
    rna = RNA(name="DNA-A", seqid="NC_001507.1", length=2588, circular=True)
    rna.add_depth("s", np.full(2588, 5))
    p = CircularPlotter(Settings(), args(layout="circular"))
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar")
    p._draw_depth(ax, rna, rna.tracks, rna.total_depth())
    p._draw_axis(ax, rna, rna.total_depth())
    assert "DNA-A\n2,588 nts" in [t.get_text() for t in ax.texts]
    plt.close(fig)


# --- CLI: --rnas by name; files still named by seqid ----------------------------

def run(tmp_path, gff, depth, yaml, *extra, name="out"):
    ex = examples_dir()
    main(["-g", os.path.join(ex, gff), "-d", os.path.join(ex, depth),
          "-y", os.path.join(ex, yaml), "-o", str(tmp_path),
          "--name", name, "-f", "svg", *extra])


def test_rnas_accepts_a_segment_name_and_files_keep_the_seqid(tmp_path):
    run(tmp_path, "tswv.gff3", "tswv.sam", "tswv.yml", "--rnas", "S", "NC_002052.1",
        name="tswv")
    assert sorted(os.listdir(tmp_path)) == ["tswv.NC_002051.1.svg", "tswv.NC_002052.1.svg"]


def test_unknown_rna_lists_seqids_with_their_names(tmp_path, caplog):
    with caplog.at_level(logging.ERROR), pytest.raises(SystemExit):
        run(tmp_path, "liyv.gff3", "liyv.bam", "liyv.yml", "--rnas", "RNA-3")
    assert "NC_003617.1 (RNA-1)" in caplog.text and "NC_003618.1 (RNA-2)" in caplog.text


def test_a_name_shared_by_two_regions_is_not_a_selector(tmp_path, caplog):
    gff = tmp_path / "dup.gff3"
    gff.write_text("##gff-version 3\n"
                   "A\t.\tregion\t1\t100\t.\t+\t.\tID=A;Name=DNA-U\n"
                   "B\t.\tregion\t1\t100\t.\t+\t.\tID=B;Name=DNA-U\n")
    ex = examples_dir()
    with caplog.at_level(logging.ERROR), pytest.raises(SystemExit):
        main(["-g", str(gff), "-d", os.path.join(ex, "sample.sam"),
              "-y", os.path.join(ex, "spec.yml"), "-o", str(tmp_path), "--rnas", "DNA-U"])
    assert "not found" in caplog.text
