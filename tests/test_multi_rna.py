"""Tests for multi-RNA support: GFF regions, per-seqid depth, CLI outputs, shared scaling."""

import argparse
import os
import textwrap

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from virplot.alignments import AlignmentError
from virplot.cli import main
from virplot.models import RNA
from virplot.parsers import parse_depth, parse_gff, parse_gff_rnas
from virplot.plotting import FIGURE_WIDTH, LinearPlotter
from virplot.settings import Settings

MULTI_GFF = textwrap.dedent("""\
    ##gff-version 3
    A\t.\tregion\t1\t1000\t.\t+\t.\tID=A
    A\t.\tCDS\t100\t400\t.\t+\t0\tID=a1;product=RdRp
    B\t.\tregion\t1\t500\t.\t+\t.\tID=B
    B\t.\tCDS\t50\t450\t.\t+\t0\tID=b1;product=CP
    A\t.\tCDS\t500\t900\t.\t+\t0\tID=a2;product=MP
    C\t.\tCDS\t1\t10\t.\t+\t0\tID=c1;product=orphan
""")


# --- GFF ------------------------------------------------------------------

def test_parse_gff_rnas_groups_by_seqid_in_region_order(tmp_path):
    gff = tmp_path / "m.gff3"
    gff.write_text(MULTI_GFF)
    rnas = parse_gff_rnas(str(gff))
    assert [r.seqid for r in rnas] == ["A", "B"]
    assert rnas[0].length == 1000 and [f.product for f in rnas[0].features] == ["RdRp", "MP"]
    assert rnas[1].length == 500 and [f.product for f in rnas[1].features] == ["CP"]
    assert rnas[0].name == "A"


def test_parse_gff_returns_first_rna_only(tmp_path):
    gff = tmp_path / "m.gff3"
    gff.write_text(MULTI_GFF)
    seq_len, feats = parse_gff(str(gff))
    assert seq_len == 1000 and len(feats) == 2


# --- depth files with several sequences -----------------------------------

def write_dep(tmp_path, rows, name="d.dep"):
    p = tmp_path / name
    p.write_text("".join(f"{s}\t{pos}\t{cov}\n" for s, pos, cov in rows))
    return str(p)


def test_parse_depth_filters_by_seqid(tmp_path):
    p = write_dep(tmp_path, [("A", 1, 5), ("A", 2, 6), ("B", 1, 9), ("B", 2, 8)])
    yA, nA = parse_depth(p, 3, seqid="A")
    yB, nB = parse_depth(p, 2, seqid="B")
    assert yA.tolist() == [5, 6, 0] and nA == 2
    assert yB.tolist() == [9, 8] and nB == 2


def test_parse_depth_without_seqid_uses_everything(tmp_path):
    p = write_dep(tmp_path, [("X", 1, 5), ("X", 2, 6)])
    y, n = parse_depth(p, 2)
    assert y.tolist() == [5, 6] and n == 2


def test_parse_depth_single_sequence_fallback_when_seqid_differs(tmp_path):
    p = write_dep(tmp_path, [("X", 1, 5), ("X", 2, 6)])
    y, n = parse_depth(p, 2, seqid="A")            # only one sequence -> used with warning
    assert y.tolist() == [5, 6]


def test_parse_depth_ambiguous_seqid_errors(tmp_path):
    p = write_dep(tmp_path, [("X", 1, 5), ("Y", 1, 6)])
    with pytest.raises(AlignmentError):
        parse_depth(p, 2, seqid="A")


def test_parse_depth_ref_override(tmp_path):
    p = write_dep(tmp_path, [("X", 1, 5), ("Y", 1, 6)])
    y, _ = parse_depth(p, 1, seqid="A", ref="Y")
    assert y.tolist() == [6]


# --- plotter cross-RNA scaling --------------------------------------------

def make_args(**over):
    base = dict(smooth=False, normalize=False, free_y=False, equal_width=False,
                legend=False, shade_breaks=False, title=False, grid=False,
                no_label=False, no_border=False, yscale="linear", linthresh=10.0,
                format="png", outdir=".", name="x")
    base.update(over)
    return argparse.Namespace(**base)


def two_rnas():
    a = RNA(name="A", seqid="A", length=100)
    a.add_depth("s", np.full(100, 40))
    b = RNA(name="B", seqid="B", length=25)
    b.add_depth("s", np.full(25, 10))
    return [a, b]


def test_prepare_shares_ymax_and_width():
    p = LinearPlotter(Settings(), make_args())
    p.prepare(two_rnas())
    assert p.shared_ymax == 40
    assert p.max_length == 100
    assert p._figure_width(two_rnas()[1]) == pytest.approx(FIGURE_WIDTH / 4)


def test_prepare_single_rna_is_noop():
    p = LinearPlotter(Settings(), make_args())
    p.prepare(two_rnas()[:1])
    assert p.shared_ymax is None and p.max_length is None and p.shared_denom is None


def test_prepare_free_y_and_equal_width_opt_out():
    p = LinearPlotter(Settings(), make_args(free_y=True, equal_width=True))
    p.prepare(two_rnas())
    assert p.shared_ymax is None and p.max_length is None


def test_prepare_normalize_uses_global_denominator():
    p = LinearPlotter(Settings(), make_args(normalize=True))
    rnas = two_rnas()
    p.prepare(rnas)
    assert p.shared_denom == 40 and p.shared_ymax == 1.0
    _, total_b = p._prepared_tracks(rnas[1])
    assert total_b.max() == pytest.approx(0.25)     # 10/40, not 10/10


# --- CLI end to end ---------------------------------------------------------

def examples_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(here), "examples")


def test_cli_multi_rna_writes_one_file_per_rna(tmp_path):
    ex = examples_dir()
    main(["-g", os.path.join(ex, "sample_multi.gff3"),
          "-d", os.path.join(ex, "sample.sam"),
          "-y", os.path.join(ex, "spec.yml"),
          "-o", str(tmp_path), "--name", "segs", "-f", "svg", "--report"])
    names = sorted(os.listdir(tmp_path))
    assert "segs.SyntheticVirus1.svg" in names and "segs.SyntheticVirus2.svg" in names
    assert "segs.SyntheticVirus2.intervals_ge1.csv" in names


def test_cli_rnas_selects_subset(tmp_path):
    ex = examples_dir()
    main(["-g", os.path.join(ex, "sample_multi.gff3"),
          "-d", os.path.join(ex, "sample.sam"),
          "-y", os.path.join(ex, "spec.yml"),
          "-o", str(tmp_path), "--name", "one", "-f", "svg", "--rnas", "SyntheticVirus2"])
    assert os.listdir(tmp_path) == ["one.svg"]      # single RNA -> plain name


def test_cli_single_rna_keeps_plain_filename(tmp_path):
    ex = examples_dir()
    main(["-g", os.path.join(ex, "sample.gff3"),
          "-d", os.path.join(ex, "sample.dep"),
          "-y", os.path.join(ex, "spec.yml"),
          "-o", str(tmp_path), "--name", "plain", "-f", "svg"])
    assert os.listdir(tmp_path) == ["plain.svg"]


def test_cli_ref_with_multiple_rnas_is_rejected(tmp_path):
    ex = examples_dir()
    with pytest.raises(SystemExit):
        main(["-g", os.path.join(ex, "sample_multi.gff3"),
              "-d", os.path.join(ex, "sample.sam"),
              "-y", os.path.join(ex, "spec.yml"),
              "-o", str(tmp_path), "--ref", "SyntheticVirus1"])
