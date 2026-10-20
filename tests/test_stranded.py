"""Stranded depth: counting by strand, the mirrored drawings, and the CLI.

docs/stranded_depth.md has the design. The two properties everything rests
on: the strands sum to the unstranded depth, and on a circle equal depth on
either side of the zero ring covers equal area.
"""

import argparse
import logging
import os
import re

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from virplot.alignments import (
    OrientedDepth,
    depth_from_alignments,
    oriented_depth_from_alignments,
    read1_reverse,
)
from virplot.cli import main
from virplot.models import RNA
from virplot.parsers import ParseError, load_depth_tracks
from virplot.plotting import (
    CIRC_R_DEPTH_BASE,
    CIRC_R_DEPTH_MAX,
    CIRC_R_STRAND_ZERO,
    Y_HEADROOM,
    CircularPlotter,
    LinearPlotter,
)
from virplot.settings import Settings

HEADER = "@HD\tVN:1.6\n@SQ\tSN:v1\tLN:20\n"
SE_FWD, SE_REV = 0, 16
R1_FWD, R1_REV = 0x1 | 0x40, 0x1 | 0x40 | 0x10
R2_FWD, R2_REV = 0x1 | 0x80, 0x1 | 0x80 | 0x10


def sam_line(flag, pos, cigar="4M", name="r"):
    qlen = sum(int(n) for n, op in re.findall(r"(\d+)([MIS=X])", cigar)) or 1
    return f"{name}\t{flag}\tv1\t{pos}\t60\t{cigar}\t*\t0\t0\t{'A' * qlen}\t{'I' * qlen}\n"


def write_sam(tmp_path, records, name="s.sam"):
    p = tmp_path / name
    p.write_text(HEADER + "".join(records))
    return str(p)


def examples_dir():
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, yscale="linear", linthresh=10.0,
                format="svg", outdir=".", name="x", layout="auto", bare_x=False)
    base.update(over)
    return argparse.Namespace(**base)


def stranded_rna(plus, minus, length=None, name="v", circular=False):
    plus, minus = np.asarray(plus, float), np.asarray(minus, float)
    rna = RNA(name=name, seqid=name, length=length or len(plus), circular=circular)
    rna.add_depth("s", plus, "+")
    rna.add_depth("s", minus, "-")
    return rna


# --- which way a read points -------------------------------------------------

@pytest.mark.parametrize("flag, reverse", [
    (SE_FWD, False), (SE_REV, True),
    (R1_FWD, False), (R1_REV, True),
    (R2_FWD, True), (R2_REV, False),      # a read 2 speaks for its read 1, flipped
])
def test_read1_reverse_flips_read_2(flag, reverse):
    assert read1_reverse(flag) is reverse


def test_strands_follow_the_library_type():
    od = OrientedDepth(np.array([1, 0]), np.array([0, 5]), 1, 1)
    plus, minus = od.strands("reverse")           # dUTP: reverse reads are + RNA
    assert plus.tolist() == [0, 5] and minus.tolist() == [1, 0]
    plus, minus = od.strands("forward")
    assert plus.tolist() == [1, 0] and minus.tolist() == [0, 5]
    with pytest.raises(ValueError):
        od.strands("no")


def test_both_mates_of_a_fragment_land_on_one_strand(tmp_path):
    p = write_sam(tmp_path, [sam_line(R1_FWD, 1, name="f"), sam_line(R2_REV, 9, name="f")])
    od = oriented_depth_from_alignments(p, 20, seqid="v1")
    assert (od.n_forward, od.n_reverse) == (2, 0)
    plus, minus = od.strands("forward")
    assert plus[:4].tolist() == [1] * 4 and plus[8:12].tolist() == [1] * 4
    assert not minus.any()


def test_strands_sum_to_the_unstranded_depth(tmp_path):
    records = [sam_line(f, pos, cigar) for f, pos, cigar in [
        (SE_FWD, 1, "5M"), (SE_REV, 3, "2M1D3M"), (R1_REV, 6, "4M"), (R2_FWD, 6, "2M2N2M"),
        (R2_REV, 15, "6M"), (SE_REV | 0x400, 2, "4M"),   # a duplicate, skipped as before
    ]]
    p = write_sam(tmp_path, records)
    od = oriented_depth_from_alignments(p, 20, seqid="v1")
    total, counted = depth_from_alignments(p, 20, seqid="v1")
    assert (od.forward + od.reverse).tolist() == total.tolist()
    assert od.counted == counted == 5


def test_glrav13_is_a_reverse_stranded_library():
    od = oriented_depth_from_alignments(os.path.join(examples_dir(), "glrav13.bam"), 17564)
    assert (od.n_reverse, od.n_forward) == (199150, 13135)
    plus, minus = od.strands("reverse")
    assert round(plus.mean()) == 1105 and round(minus.mean()) == 68
    total, _ = depth_from_alignments(os.path.join(examples_dir(), "glrav13.bam"), 17564)
    assert ((plus + minus) == total).all()


# --- loading -----------------------------------------------------------------

def test_load_depth_tracks_splits_alignments_by_strand(tmp_path):
    p = write_sam(tmp_path, [sam_line(SE_REV, 1), sam_line(SE_FWD, 1)])
    split = load_depth_tracks(p, 20, stranded="reverse", seqid="v1")
    assert [s for s, _ in split.tracks] == ["+", "-"]
    assert split.reverse_fraction == 0.5
    whole = load_depth_tracks(p, 20, seqid="v1")
    assert [s for s, _ in whole.tracks] == [None]
    assert whole.tracks[0][1][:4].tolist() == [2] * 4


def test_a_depth_file_cannot_be_split(tmp_path):
    dep = tmp_path / "d.dep"
    dep.write_text("".join(f"v1\t{i}\t3\n" for i in range(1, 21)))
    with pytest.raises(ParseError, match="no strand"):
        load_depth_tracks(str(dep), 20, stranded="reverse", seqid="v1")


# --- model -------------------------------------------------------------------

def test_rna_keeps_strands_apart_and_sums_them():
    rna = stranded_rna([3, 4], [1, 0])
    assert rna.stranded and rna.sample_labels == ["s"]
    assert [t.values.tolist() for t in rna.strand_tracks("-")] == [[1, 0]]
    assert rna.total_depth().tolist() == [4, 4]
    with pytest.raises(ValueError):
        rna.add_depth("s", np.zeros(2), "x")


def test_unstranded_rna_is_not_stranded():
    rna = RNA(name="v", length=2)
    rna.add_depth("s", np.array([1, 2]))
    assert not rna.stranded


# --- linear: mirrored on one linear scale --------------------------------------

def test_linear_axis_is_one_scale_with_each_side_as_tall_as_its_strand():
    rna = stranded_rna([0, 400, 1000, 200], [0, 50, 100, 0])
    p = LinearPlotter(Settings(), args())
    fig, ax = plt.subplots()
    pt, ptot, mt, mtot = p._prepared_strands(rna)
    p._draw_stranded_depth(ax, rna, pt, ptot, mt, mtot)
    labels = p._style_stranded_axis(ax, ptot, mtot, rna)
    assert ax.get_ylim() == pytest.approx((-100 * Y_HEADROOM, 1000 * Y_HEADROOM))
    assert [t.get_text() for t in labels] == ["+ strand", "− strand"]
    fig.canvas.draw()
    shown = [t.get_text() for t in ax.get_yticklabels() if t.get_text()]
    assert shown and not any(s.startswith(("-", "−")) for s in shown)
    plt.close(fig)


def test_segments_share_each_strand_limit():
    a = stranded_rna([10, 900], [5, 40])
    b = stranded_rna([300, 20], [120, 0])
    p = LinearPlotter(Settings(), args())
    p.prepare([a, b])
    assert p.shared_strand_ymax == (900, 120)
    assert p.shared_ymax is None                  # the unstranded limit stays unset


def test_normalize_divides_both_strands_by_one_number():
    rna = stranded_rna([0, 800], [200, 0])
    p = LinearPlotter(Settings(), args(normalize=True))
    _, ptot, _, mtot = p._prepared_strands(rna)
    assert ptot.max() == 1.0 and mtot.max() == pytest.approx(0.25)


# --- circular: mirrored, area-true ---------------------------------------------

def test_zero_ring_splits_the_band_into_equal_areas():
    outer = CIRC_R_DEPTH_MAX ** 2 - CIRC_R_STRAND_ZERO ** 2
    inner = CIRC_R_STRAND_ZERO ** 2 - CIRC_R_DEPTH_BASE ** 2
    assert outer == pytest.approx(inner)


@pytest.mark.parametrize("depth", [0, 1, 37, 250, 999, 1000])
def test_equal_depth_covers_equal_area_on_both_sides(depth):
    r_out = CircularPlotter._strand_radius(np.array([depth]), 1000, 1)[0]
    r_in = CircularPlotter._strand_radius(np.array([depth]), 1000, -1)[0]
    z2 = CIRC_R_STRAND_ZERO ** 2
    assert r_out ** 2 - z2 == pytest.approx(z2 - r_in ** 2)
    assert r_in <= CIRC_R_STRAND_ZERO <= r_out


def test_scale_maximum_reaches_both_edges_of_the_band():
    v = np.array([1000.0])
    assert CircularPlotter._strand_radius(v, 1000, 1)[0] == pytest.approx(CIRC_R_DEPTH_MAX)
    assert CircularPlotter._strand_radius(v, 1000, -1)[0] == pytest.approx(CIRC_R_DEPTH_BASE)


def test_circular_scale_reads_zero_plus_and_minus():
    rna = stranded_rna(np.full(300, 40), np.full(300, 10), circular=True)
    p = CircularPlotter(Settings(), args(layout="circular"))
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar")
    pt, ptot, mt, mtot = p._prepared_strands(rna)
    layers, names = p._draw_stranded_depth(ax, rna, pt, ptot, mt, mtot)
    p._draw_axis(ax, rna, ptot + mtot)
    texts = [t.get_text() for t in ax.texts]
    assert " 0" in texts and " + 42× depth" in texts and " − 42×" in texts
    assert names == ["+ strand", "− strand"] and len(layers) == 2
    plt.close(fig)


# --- CLI ---------------------------------------------------------------------

def run(tmp_path, *extra, gff="glrav13.gff3", depth="glrav13.bam", yaml="glrav13.yml"):
    ex = examples_dir()
    main(["-g", os.path.join(ex, gff), "-d", os.path.join(ex, depth),
          "-y", os.path.join(ex, yaml), "-o", str(tmp_path), "--name", "out",
          "-f", "svg", *extra])


def test_cli_stranded_run_reports_both_strands(tmp_path, caplog):
    with caplog.at_level(logging.INFO):
        run(tmp_path, "--stranded", "reverse", "--smooth")
    assert os.listdir(tmp_path) == ["out.svg"]
    assert "93.8% reverse" in caplog.text
    assert "+ strand mean 1105x, - strand mean 68x" in caplog.text
    assert "looks stranded" not in caplog.text


def test_cli_says_when_an_unstranded_run_looks_stranded(tmp_path, caplog):
    with caplog.at_level(logging.WARNING):
        run(tmp_path, "--smooth")
    assert "looks stranded" in caplog.text and "--stranded reverse" in caplog.text


def test_cli_says_when_a_stranded_run_splits_evenly(tmp_path, caplog):
    with caplog.at_level(logging.WARNING):
        run(tmp_path, "--stranded", "reverse", gff="byv.gff3", depth="byv.sam", yaml="byv.yml")
    assert "splits nearly evenly" in caplog.text


def test_cli_refuses_to_split_a_depth_file(tmp_path, caplog):
    with caplog.at_level(logging.ERROR), pytest.raises(SystemExit):
        run(tmp_path, "--stranded", "reverse", gff="sample.gff3", depth="sample.dep",
            yaml="spec.yml")
    assert "records no strand" in caplog.text


def test_cli_stranded_circular_and_segmented_runs_render(tmp_path):
    run(tmp_path, "--stranded", "forward", "--layout", "circular", "--legend",
        gff="tgmv.gff3", depth="tgmv.sam", yaml="tgmv.yml")
    assert sorted(os.listdir(tmp_path)) == ["out.NC_001507.1.svg", "out.NC_001508.1.svg"]
