"""Tests for virplot.parsers."""

import textwrap
import numpy as np
import pytest

from virplot.parsers import ParseError, parse_gff, parse_depth


# parse_gff

def test_parse_gff_basic(tmp_path):
    gff = tmp_path / "test.gff3"
    gff.write_text(textwrap.dedent("""\
        ##gff-version 3
        seq1\t.\tregion\t1\t1000\t.\t+\t.\tID=seq1
        seq1\t.\tCDS\t100\t400\t.\t+\t0\tID=cds1;product=RdRp
        seq1\t.\tCDS\t500\t900\t.\t+\t0\tID=cds2;product=CP
    """))
    seq_len, features = parse_gff(str(gff))
    assert seq_len == 1000
    assert len(features) == 2
    assert features[0].product == "RdRp"
    assert features[0].start == 100
    assert features[0].end == 400
    assert features[1].product == "CP"


def test_parse_gff_unknown_product(tmp_path):
    gff = tmp_path / "test.gff3"
    gff.write_text(textwrap.dedent("""\
        seq1\t.\tregion\t1\t500\t.\t+\t.\tID=seq1
        seq1\t.\tCDS\t10\t200\t.\t+\t0\tID=cds1
    """))
    _, features = parse_gff(str(gff))
    assert features[0].product == "unknown"


def test_parse_gff_no_region(tmp_path):
    gff = tmp_path / "test.gff3"
    gff.write_text("seq1\t.\tCDS\t10\t200\t.\t+\t0\tID=cds1;product=X\n")
    with pytest.raises(ParseError):
        parse_gff(str(gff))


def test_parse_gff_skips_non_cds(tmp_path):
    gff = tmp_path / "test.gff3"
    gff.write_text(textwrap.dedent("""\
        seq1\t.\tregion\t1\t1000\t.\t+\t.\tID=seq1
        seq1\t.\tgene\t100\t400\t.\t+\t.\tID=gene1
        seq1\t.\tCDS\t100\t400\t.\t+\t0\tID=cds1;product=RdRp
    """))
    _, features = parse_gff(str(gff))
    assert len(features) == 1


# parse_depth

def test_parse_depth_basic(tmp_path):
    dep = tmp_path / "test.dep"
    dep.write_text("seq1\t1\t10\nseq1\t2\t20\nseq1\t3\t5\n")
    y, n = parse_depth(str(dep), 5)
    assert n == 3
    assert y.shape == (5,)
    np.testing.assert_array_equal(y, [10, 20, 5, 0, 0])


def test_parse_depth_out_of_range_ignored(tmp_path):
    dep = tmp_path / "test.dep"
    dep.write_text("seq1\t1\t10\nseq1\t99\t50\n")
    y, n = parse_depth(str(dep), 3)
    assert n == 2
    np.testing.assert_array_equal(y, [10, 0, 0])


def test_parse_depth_malformed_columns(tmp_path):
    dep = tmp_path / "bad.dep"
    dep.write_text("seq1\t1\n")  # only 2 columns
    with pytest.raises(ParseError):
        parse_depth(str(dep), 10)


def test_parse_depth_non_integer(tmp_path):
    dep = tmp_path / "bad.dep"
    dep.write_text("seq1\tabc\t10\n")
    with pytest.raises(ParseError):
        parse_depth(str(dep), 10)


def test_parse_depth_empty(tmp_path):
    dep = tmp_path / "empty.dep"
    dep.write_text("")
    y, n = parse_depth(str(dep), 5)
    assert n == 0
    np.testing.assert_array_equal(y, [0, 0, 0, 0, 0])


def test_parse_depth_boundary_positions(tmp_path):
    """Position at seq_len is valid; position 0 is out of range (1-based)."""
    dep = tmp_path / "boundary.dep"
    dep.write_text("seq1\t0\t99\nseq1\t5\t42\n")
    y, n = parse_depth(str(dep), 5)
    assert n == 2
    # pos 0 is out of 1-based range → ignored; pos 5 == seq_len → valid
    np.testing.assert_array_equal(y, [0, 0, 0, 0, 42])


def test_parse_gff_skips_malformed_columns(tmp_path):
    """Lines with fewer than 9 tab-separated columns are silently skipped."""
    gff = tmp_path / "bad.gff3"
    gff.write_text(
        "seq1\t.\tregion\t1\t500\t.\t+\t.\tID=seq1\n"
        "this line only has three\tcolumns\there\n"
        "seq1\t.\tCDS\t10\t200\t.\t+\t0\tID=cds1;product=RdRp\n"
    )
    seq_len, features = parse_gff(str(gff))
    assert seq_len == 500
    assert len(features) == 1


# --- row order and attribute escapes --------------------------------------------

def test_parse_gff_rows_before_their_region_line_are_kept(tmp_path):
    # GFF3 does not require the region line first; NCBI usually writes it first
    # but hand-edited and concatenated files often do not
    p = tmp_path / "t.gff3"
    p.write_text(textwrap.dedent("""\
        ##gff-version 3
        s1	.	CDS	10	100	.	+	0	ID=c1;product=CP
        s1	.	stem_loop	1	30	.	+	.	ID=sl
        s1	.	region	1	1000	.	+	.	ID=s1
    """))
    from virplot.parsers import parse_gff_rnas
    (rna,) = parse_gff_rnas(str(p))
    assert [f.product for f in rna.features] == ["CP"]
    assert len(rna.noncoding) == 1


def test_parse_gff_percent_decodes_attribute_values(tmp_path):
    p = tmp_path / "t.gff3"
    p.write_text(textwrap.dedent("""\
        ##gff-version 3
        s1	.	region	1	1000	.	+	.	ID=s1
        s1	.	CDS	10	100	.	+	0	ID=c1;product=polyprotein%2C putative;Note=a%3Bb
        s1	.	mat_peptide	10	60	.	+	.	Parent=c1;product=P1%2FP2
    """))
    from virplot.parsers import parse_gff_rnas
    (rna,) = parse_gff_rnas(str(p))
    (f,) = rna.features
    assert f.product == "polyprotein, putative"
    assert f.domains[0].product == "P1/P2"


def test_wrap_feature_on_linear_region_warns(tmp_path):
    import logging
    from virplot.parsers import parse_gff_rnas
    p = tmp_path / "t.gff3"
    p.write_text("##gff-version 3\ns1\t.\tregion\t1\t1000\t.\t+\t.\tID=s1\n"
                 "s1\t.\tCDS\t900\t100\t.\t+\t0\tID=c;product=Rep\n")
    records = []
    h = logging.Handler(); h.emit = records.append
    logging.getLogger("virplot.parsers").addHandler(h)
    try:
        (rna,) = parse_gff_rnas(str(p))
    finally:
        logging.getLogger("virplot.parsers").removeHandler(h)
    assert any("cross the origin" in r.getMessage() for r in records)
    assert rna.feature_spans(rna.features[0]) == [(900, 1000), (1, 100)]
