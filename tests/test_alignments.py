"""Tests for virplot.alignments (SAM/BAM -> depth)."""

import gzip
import re
import struct
import textwrap

import numpy as np
import pytest

from virplot.alignments import (
    AlignmentError,
    depth_from_alignments,
    is_alignment_file,
    open_alignments,
)
from virplot.parsers import gff_seqid, load_depth

HEADER = "@HD\tVN:1.6\n@SQ\tSN:v1\tLN:20\n@SQ\tSN:v2\tLN:10\n"


def sam_line(flag, ref, pos, mapq, cigar, name="r"):
    qlen = sum(int(n) for n, op in re.findall(r"(\d+)([MIS=X])", cigar)) or 1
    return f"{name}\t{flag}\t{ref}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{'A' * qlen}\t{'I' * qlen}\n"


def write_sam(tmp_path, records, header=HEADER, name="t.sam"):
    p = tmp_path / name
    p.write_text(header + "".join(records))
    return str(p)


def brute_force(records, seq_len, ref="v1", min_mapq=0):
    """Independent, position-by-position depth for cross-checking."""
    y = np.zeros(seq_len, int)
    for line in records:
        f = line.split("\t")
        if f[2] != ref or int(f[1]) & 0x704 or int(f[4]) < min_mapq:
            continue
        pos = int(f[3]) - 1
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            n = int(n)
            if op in "M=X":
                for k in range(pos, pos + n):
                    if 0 <= k < seq_len:
                        y[k] += 1
            if op in "MDN=X":
                pos += n
    return y


# --- basic counting -------------------------------------------------------

def test_single_read_simple_match(tmp_path):
    p = write_sam(tmp_path, [sam_line(0, "v1", 5, 60, "4M")])
    y, n = depth_from_alignments(p, 20, seqid="v1")
    assert n == 1
    assert y.tolist() == [0, 0, 0, 0, 1, 1, 1, 1] + [0] * 12


def test_overlapping_reads_accumulate(tmp_path):
    p = write_sam(tmp_path, [sam_line(0, "v1", 1, 60, "6M"), sam_line(16, "v1", 4, 60, "6M")])
    y, _ = depth_from_alignments(p, 20, seqid="v1")
    assert y[:10].tolist() == [1, 1, 1, 2, 2, 2, 1, 1, 1, 0]


def test_cigar_ops_soft_clip_insertion_deletion_skip(tmp_path):
    recs = [
        sam_line(0, "v1", 3, 60, "2S3M"),        # clip does not consume ref: 3..5
        sam_line(0, "v1", 8, 60, "2M2I2M"),      # insertion: 8..11
        sam_line(0, "v1", 12, 60, "2M3D2M"),     # deletion gap not counted: 12,13,17,18
        sam_line(0, "v1", 1, 60, "1M5N1M"),      # ref skip: 1 and 7
    ]
    p = write_sam(tmp_path, recs)
    y, n = depth_from_alignments(p, 20, seqid="v1")
    assert n == 4
    assert np.array_equal(y, brute_force(recs, 20))
    assert y[13:16].tolist() == [0, 0, 0]        # deleted bases 14..16 have no depth
    assert y[16] == 1                            # 17 is covered again


def test_read_past_end_is_clipped(tmp_path):
    p = write_sam(tmp_path, [sam_line(0, "v1", 18, 60, "10M")])
    y, _ = depth_from_alignments(p, 20, seqid="v1")
    assert y[17:].tolist() == [1, 1, 1]
    assert len(y) == 20


# --- filters ----------------------------------------------------------------

def test_default_flag_filters_match_samtools_depth(tmp_path):
    recs = [
        sam_line(0, "v1", 1, 60, "5M"),
        sam_line(4, "*", 0, 0, "*"),             # unmapped
        sam_line(256, "v1", 1, 60, "5M"),        # secondary
        sam_line(512, "v1", 1, 60, "5M"),        # QC fail
        sam_line(1024, "v1", 1, 60, "5M"),       # duplicate
        sam_line(2048, "v1", 1, 60, "5M"),       # supplementary IS counted (as samtools)
    ]
    p = write_sam(tmp_path, recs)
    y, n = depth_from_alignments(p, 20, seqid="v1")
    assert n == 2
    assert y[:5].tolist() == [2] * 5


def test_min_mapq(tmp_path):
    recs = [sam_line(0, "v1", 1, 60, "5M"), sam_line(0, "v1", 1, 5, "5M")]
    p = write_sam(tmp_path, recs)
    assert depth_from_alignments(p, 20, seqid="v1")[0][0] == 2
    assert depth_from_alignments(p, 20, seqid="v1", min_mapq=20)[0][0] == 1


def test_only_target_reference_counted(tmp_path):
    recs = [sam_line(0, "v1", 1, 60, "5M"), sam_line(0, "v2", 1, 60, "5M")]
    p = write_sam(tmp_path, recs)
    y, n = depth_from_alignments(p, 20, seqid="v1")
    assert n == 1


# --- reference resolution -------------------------------------------------

def test_ref_override(tmp_path):
    p = write_sam(tmp_path, [sam_line(0, "v2", 1, 60, "5M")])
    y, n = depth_from_alignments(p, 10, seqid="v1", ref="v2")
    assert n == 1 and y[0] == 1


def test_unknown_ref_override_errors(tmp_path):
    p = write_sam(tmp_path, [sam_line(0, "v1", 1, 60, "5M")])
    with pytest.raises(AlignmentError):
        depth_from_alignments(p, 20, seqid="v1", ref="nope")


def test_single_reference_used_when_seqid_absent(tmp_path):
    header = "@SQ\tSN:onlyone\tLN:20\n"
    p = write_sam(tmp_path, [sam_line(0, "onlyone", 1, 60, "5M")], header=header)
    y, n = depth_from_alignments(p, 20, seqid="something_else")
    assert n == 1


def test_ambiguous_reference_errors(tmp_path):
    p = write_sam(tmp_path, [sam_line(0, "v1", 1, 60, "5M")])
    with pytest.raises(AlignmentError):
        depth_from_alignments(p, 20, seqid="not_in_header")


def test_malformed_record_errors(tmp_path):
    p = tmp_path / "bad.sam"
    p.write_text(HEADER + "r\t0\tv1\n")
    with pytest.raises(AlignmentError):
        depth_from_alignments(str(p), 20, seqid="v1")


# --- BAM ------------------------------------------------------------------

def bam_bytes(refs, records):
    """Minimal BAM writer: refs = [(name, length)], records = (flag, ref_idx, pos1, mapq, cigar)."""
    ops = "MIDNSHP=X"
    out = bytearray(b"BAM\x01")
    text = b"@HD\tVN:1.6\n"
    out += struct.pack("<i", len(text)) + text
    out += struct.pack("<i", len(refs))
    for name, length in refs:
        nm = name.encode() + b"\x00"
        out += struct.pack("<i", len(nm)) + nm + struct.pack("<i", length)
    for flag, ref_idx, pos1, mapq, cigar in records:
        cig = [(int(n), ops.index(op)) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
        qlen = sum(n for n, o in cig if ops[o] in "MIS=X")
        name = b"r\x00"
        body = struct.pack("<iiBBHHHiiii", ref_idx, pos1 - 1, len(name), mapq, 0,
                           len(cig), flag, qlen, -1, -1, 0)
        body += name
        body += b"".join(struct.pack("<I", (n << 4) | o) for n, o in cig)
        body += b"\x00" * ((qlen + 1) // 2) + b"\xff" * qlen
        out += struct.pack("<i", len(body)) + body
    return gzip.compress(bytes(out))


def test_bam_matches_sam(tmp_path):
    sam_recs = [
        sam_line(0, "v1", 3, 60, "2S3M"),
        sam_line(16, "v1", 8, 60, "2M2I2M"),
        sam_line(0, "v1", 12, 60, "2M3D2M"),
        sam_line(256, "v1", 1, 60, "5M"),
        sam_line(0, "v2", 1, 60, "5M"),
    ]
    bam_recs = [(0, 0, 3, 60, "2S3M"), (16, 0, 8, 60, "2M2I2M"),
                (0, 0, 12, 60, "2M3D2M"), (256, 0, 1, 60, "5M"), (0, 1, 1, 60, "5M")]
    sam = write_sam(tmp_path, sam_recs)
    bam = tmp_path / "t.bam"
    bam.write_bytes(bam_bytes([("v1", 20), ("v2", 10)], bam_recs))

    af = open_alignments(str(bam))
    assert af.references == {"v1": 20, "v2": 10}

    ys, ns = depth_from_alignments(sam, 20, seqid="v1")
    yb, nb = depth_from_alignments(str(bam), 20, seqid="v1")
    assert ns == nb == 3
    assert np.array_equal(ys, yb)


def test_gzip_but_not_bam_errors(tmp_path):
    p = tmp_path / "x.bam"
    p.write_bytes(gzip.compress(b"not a bam"))
    with pytest.raises(AlignmentError):
        open_alignments(str(p))


# --- detection and dispatch ----------------------------------------------

def test_is_alignment_file_by_extension_and_content(tmp_path):
    sam = write_sam(tmp_path, [], name="a.sam")
    dep = tmp_path / "a.dep"
    dep.write_text("v1\t1\t3\n")
    unk_sam = write_sam(tmp_path, [], name="noext")
    unk_dep = tmp_path / "noext2"
    unk_dep.write_text("v1\t1\t3\n")
    assert is_alignment_file(sam)
    assert not is_alignment_file(str(dep))
    assert is_alignment_file(unk_sam)
    assert not is_alignment_file(str(unk_dep))


def test_load_depth_dispatch(tmp_path):
    sam = write_sam(tmp_path, [sam_line(0, "v1", 1, 60, "5M")])
    dep = tmp_path / "a.dep"
    dep.write_text("".join(f"v1\t{i}\t{7}\n" for i in range(1, 21)))
    y1, n1, k1 = load_depth(sam, 20, seqid="v1")
    y2, n2, k2 = load_depth(str(dep), 20)
    assert k1 == "alignments" and n1 == 1 and y1[0] == 1
    assert k2 == "depth" and n2 == 20 and y2[0] == 7


def test_gff_seqid(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(textwrap.dedent("""\
        ##gff-version 3
        seqX\t.\tregion\t1\t100\t.\t+\t.\tID=seqX
        seqX\t.\tCDS\t1\t50\t.\t+\t0\tID=c;product=P
    """))
    assert gff_seqid(str(gff)) == "seqX"


# --- randomised cross-check against brute force ---------------------------

def test_random_reads_match_brute_force(tmp_path):
    rng = np.random.default_rng(7)
    ops = ["10M", "3S7M", "4M2I4M", "3M2D5M", "2M4N4M", "5M1X4M", "6=3X"]
    recs = []
    for i in range(300):
        flag = int(rng.choice([0, 16, 256, 1024], p=[0.45, 0.45, 0.05, 0.05]))
        ref = "v1" if rng.random() < 0.9 else "v2"
        pos = int(rng.integers(1, 20))
        mapq = int(rng.integers(0, 61))
        recs.append(sam_line(flag, ref, pos, mapq, str(rng.choice(ops)), name=f"r{i}"))
    p = write_sam(tmp_path, recs)
    for mq in (0, 30):
        y, _ = depth_from_alignments(p, 20, seqid="v1", min_mapq=mq)
        assert np.array_equal(y, brute_force(recs, 20, min_mapq=mq))
