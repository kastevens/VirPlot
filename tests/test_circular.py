"""Tests for circular topology: detection, depth wrapping, and the GRBV example."""

import os
import textwrap

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from virplot.alignments import depth_from_alignments
from virplot.cli import main
from virplot.parsers import parse_depth, parse_gff_rnas

HEADER = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:20\n"


def sam(tmp_path, records, name="c.sam"):
    p = tmp_path / name
    p.write_text(HEADER + "".join(
        f"r{i}\t{flag}\tc1\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*\n"
        for i, (flag, pos, cigar) in enumerate(records)))
    return str(p)


# --- Is_circular detection -------------------------------------------------

def test_is_circular_true_sets_flag(tmp_path):
    gff = tmp_path / "c.gff3"
    gff.write_text(
        "A\t.\tregion\t1\t100\t.\t+\t.\tID=A;Is_circular=true\n"
        "B\t.\tregion\t1\t50\t.\t+\t.\tID=B\n"
        "C\t.\tregion\t1\t50\t.\t+\t.\tID=C;Is_circular=false\n"
    )
    rnas = {r.seqid: r for r in parse_gff_rnas(str(gff))}
    assert rnas["A"].circular is True
    assert rnas["B"].circular is False
    assert rnas["C"].circular is False


def test_grbv_example_is_circular():
    ex = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")
    (rna,) = parse_gff_rnas(os.path.join(ex, "grbv.gff3"))
    assert rna.circular is True
    assert rna.length == 3206
    assert rna.seqid == "NC_022002.1"
    # curated per the ICTV convention: ORF name in gene=, function in product=.
    # V3 and C3 have no known function, so the curator put the ORF name in both
    # fields; the outside label is dropped rather than writing V3 twice.
    assert [f.gene for f in rna.features] == ["V2", "V1", None, "C2", "C1", None]
    assert [f.product for f in rna.features] == ["putative MP", "CP", "V3", "Rep", "RepA", "C3"]
    assert [f.strand for f in rna.features] == ["+", "+", "+", "-", "-", "-"]


# --- wrapping in alignment depth ------------------------------------------

def test_read_crossing_origin_wraps(tmp_path):
    p = sam(tmp_path, [(0, 18, "6M")])            # 18,19,20 then 1,2,3
    lin, _ = depth_from_alignments(p, 20, seqid="c1")
    cir, _ = depth_from_alignments(p, 20, seqid="c1", circular=True)
    assert lin.sum() == 3 and lin[17:].tolist() == [1, 1, 1] and lin[0] == 0
    assert cir.sum() == 6
    assert cir[17:].tolist() == [1, 1, 1] and cir[:3].tolist() == [1, 1, 1]


def test_position_past_end_wraps(tmp_path):
    p = sam(tmp_path, [(0, 23, "4M")])            # padded-reference POS: 23 -> 3
    y, n = depth_from_alignments(p, 20, seqid="c1", circular=True)
    assert n == 1
    assert y[2:6].tolist() == [1, 1, 1, 1] and y.sum() == 4


def test_read_longer_than_genome_covers_once(tmp_path):
    p = sam(tmp_path, [(0, 5, "50M")])
    y, _ = depth_from_alignments(p, 20, seqid="c1", circular=True)
    assert y.tolist() == [1] * 20


def test_wrapping_conserves_total_covered_bases(tmp_path):
    rng = np.random.default_rng(3)
    recs = [(0, int(rng.integers(1, 21)), "7M") for _ in range(50)]
    p = sam(tmp_path, recs)
    y, n = depth_from_alignments(p, 20, seqid="c1", circular=True)
    assert n == 50
    assert y.sum() == 50 * 7                      # nothing lost off the end


def test_circular_never_loses_more_than_linear(tmp_path):
    recs = [(0, p0, "6M") for p0 in range(1, 21)]
    p = sam(tmp_path, recs)
    lin, _ = depth_from_alignments(p, 20, seqid="c1")
    cir, _ = depth_from_alignments(p, 20, seqid="c1", circular=True)
    assert cir.sum() > lin.sum()
    assert np.all(cir >= lin)
    assert cir.tolist() == [6] * 20               # uniform tiling of the circle


def test_deletion_still_not_counted_when_circular(tmp_path):
    p = sam(tmp_path, [(0, 17, "2M2D2M")])        # 17,18 skip 19,20 cover 1,2
    y, _ = depth_from_alignments(p, 20, seqid="c1", circular=True)
    assert y[16:18].tolist() == [1, 1]
    assert y[18:20].tolist() == [0, 0]
    assert y[:2].tolist() == [1, 1]


# --- wrapping in depth text files -----------------------------------------

def test_depth_file_positions_past_end_wrap(tmp_path):
    p = tmp_path / "d.dep"
    p.write_text("".join(f"c1\t{i}\t{i}\n" for i in range(1, 25)))
    lin, _ = parse_depth(str(p), 20)
    cir, _ = parse_depth(str(p), 20, circular=True)
    assert lin[:4].tolist() == [1, 2, 3, 4]
    assert cir[:4].tolist() == [1 + 21, 2 + 22, 3 + 23, 4 + 24]


# --- CLI -------------------------------------------------------------------

def examples_dir():
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def test_cli_grbv_runs_and_reports_full_coverage(tmp_path):
    ex = examples_dir()
    main(["-g", os.path.join(ex, "grbv.gff3"), "-d", os.path.join(ex, "grbv.sam"),
          "-y", os.path.join(ex, "grbv.yml"), "-o", str(tmp_path),
          "--name", "grbv", "-f", "svg"])
    assert os.listdir(tmp_path) == ["grbv.svg"]


def test_topology_override_changes_origin_depth(tmp_path):
    ex = examples_dir()
    gff = os.path.join(ex, "grbv.gff3")
    sam_path = os.path.join(ex, "grbv.sam")
    circ, _ = depth_from_alignments(sam_path, 3206, seqid="NC_022002.1", circular=True)
    lin, _ = depth_from_alignments(sam_path, 3206, seqid="NC_022002.1", circular=False)
    # the first read-length of the genome is where linear treatment loses reads
    assert lin[:150].mean() < 0.5 * circ[:150].mean()
    # away from the origin the two agree exactly
    assert np.array_equal(lin[300:3000], circ[300:3000])
