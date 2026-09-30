"""Tests for read mapping with bowtie2 (analysis.py) and the -U/-1/-2/-x/-p CLI."""

import os
import shutil
import stat
import subprocess
import textwrap

import matplotlib
import pytest

matplotlib.use("Agg")

from virplot import analysis
from virplot.analysis import (
    MappingError,
    ReadSet,
    default_read_label,
    ensure_bowtie2_index,
    fasta_lengths,
    map_reads,
)
from virplot.cli import main

EX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


# --- labels and read sets --------------------------------------------------

def test_default_read_label_strips_suffixes():
    assert default_read_label(["/a/vine7_R1.fastq.gz"]) == "vine7"
    assert default_read_label(["s.fq"]) == "s"
    assert default_read_label(["lane1_1.fastq", "lane2_1.fastq"]) == "lane1"
    assert default_read_label(["reads.fastq.gz"]) == "reads"


def test_readset_validation(tmp_path):
    f = tmp_path / "r.fq"; f.write_text("@r\nA\n+\nI\n")
    ReadSet("ok", unpaired=[str(f)]).validate()
    with pytest.raises(MappingError):
        ReadSet("both", unpaired=[str(f)], mate1=[str(f)], mate2=[str(f)]).validate()
    with pytest.raises(MappingError):
        ReadSet("lopsided", mate1=[str(f), str(f)], mate2=[str(f)]).validate()
    with pytest.raises(MappingError):
        ReadSet("missing", unpaired=[str(tmp_path / "nope.fq")]).validate()
    with pytest.raises(MappingError):
        ReadSet("empty").validate()


def test_fasta_lengths(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">seqA some description\nACGT\nAC\n>seqB\nGGG\n")
    assert fasta_lengths(str(fa)) == {"seqA": 6, "seqB": 3}


# --- index handling ---------------------------------------------------------

class FakeRun:
    """Stands in for subprocess.run: records commands, fakes bowtie2-build output."""

    def __init__(self, returncode=0, stderr=""):
        self.calls = []; self.returncode = returncode; self.stderr = stderr

    def __call__(self, cmd, check=False, capture_output=True, text=True):
        self.calls.append(cmd)
        if cmd[0] == "bowtie2-build" and self.returncode == 0:
            base = cmd[-1]
            for suffix in ("1", "2", "3", "4", "rev.1", "rev.2"):
                open(f"{base}.{suffix}.bt2", "w").close()
        return subprocess.CompletedProcess(cmd, self.returncode, stdout="", stderr=self.stderr)


def with_tools_on_path(fn):
    """Run fn with shutil.which pretending bowtie2 tools exist."""
    real = shutil.which
    shutil.which = lambda name, *a, **k: f"/fake/{name}" if name.startswith("bowtie2") else real(name, *a, **k)
    try:
        return fn()
    finally:
        shutil.which = real


def test_index_beside_fasta_is_reused(tmp_path):
    fa = tmp_path / "ref.fa"; fa.write_text(">c1\nACGT\n")
    for suffix in ("1", "2", "3", "4", "rev.1", "rev.2"):
        (tmp_path / f"ref.fa.{suffix}.bt2").write_text("")
    run = FakeRun()
    base = with_tools_on_path(lambda: ensure_bowtie2_index(str(fa), str(tmp_path / "out"), runner=run))
    assert base == str(fa) and run.calls == []


def test_index_is_built_once_then_cached(tmp_path):
    fa = tmp_path / "ref.fa"; fa.write_text(">c1\nACGT\n")
    out = tmp_path / "out"
    run = FakeRun()
    base = with_tools_on_path(lambda: ensure_bowtie2_index(str(fa), str(out), threads=3, runner=run))
    assert base == str(out / "ref.bt2" / "ref")
    assert run.calls[0][:3] == ["bowtie2-build", "--threads", "3"]
    base2 = with_tools_on_path(lambda: ensure_bowtie2_index(str(fa), str(out), runner=run))
    assert base2 == base and len(run.calls) == 1            # cached


def test_missing_reference_errors(tmp_path):
    with pytest.raises(MappingError):
        ensure_bowtie2_index(str(tmp_path / "none.fa"), str(tmp_path), runner=FakeRun())


# --- bowtie2 command --------------------------------------------------------

def test_map_reads_unpaired_command(tmp_path):
    f = tmp_path / "s.fq"; f.write_text("@r\nA\n+\nI\n")
    run = FakeRun(stderr="10 reads; of these:\n  100.00% overall alignment rate\n")
    out = with_tools_on_path(lambda: map_reads(ReadSet("s", unpaired=[str(f)]), "idx",
                                              str(tmp_path / "s.sam"), threads=4,
                                              extra_args=["--local"], runner=run))
    assert out == str(tmp_path / "s.sam")
    assert run.calls[0] == ["bowtie2", "-p", "4", "-x", "idx", "-U", str(f), "--local",
                            "-S", str(tmp_path / "s.sam")]


def test_map_reads_paired_command_joins_files_with_commas(tmp_path):
    a1 = tmp_path / "a_1.fq"; a2 = tmp_path / "a_2.fq"; b1 = tmp_path / "b_1.fq"; b2 = tmp_path / "b_2.fq"
    for f in (a1, a2, b1, b2): f.write_text("@r\nA\n+\nI\n")
    run = FakeRun()
    with_tools_on_path(lambda: map_reads(ReadSet("a", mate1=[str(a1), str(b1)], mate2=[str(a2), str(b2)]),
                                        "idx", str(tmp_path / "a.sam"), runner=run))
    cmd = run.calls[0]
    assert cmd[cmd.index("-1") + 1] == f"{a1},{b1}" and cmd[cmd.index("-2") + 1] == f"{a2},{b2}"
    assert "-U" not in cmd


def test_bowtie2_failure_is_reported_with_stderr_tail(tmp_path):
    f = tmp_path / "s.fq"; f.write_text("@r\nA\n+\nI\n")
    run = FakeRun(returncode=1, stderr="Error: Could not open read file\n(ERR): bowtie2-align exited with value 1\n")
    with pytest.raises(MappingError) as exc:
        with_tools_on_path(lambda: map_reads(ReadSet("s", unpaired=[str(f)]), "idx",
                                            str(tmp_path / "s.sam"), runner=run))
    assert "exited with value 1" in str(exc.value)


def test_missing_tool_is_a_clear_error(tmp_path):
    f = tmp_path / "s.fq"; f.write_text("@r\nA\n+\nI\n")
    real = shutil.which
    shutil.which = lambda name, *a, **k: None
    try:
        with pytest.raises(MappingError) as exc:
            map_reads(ReadSet("s", unpaired=[str(f)]), "idx", str(tmp_path / "s.sam"), runner=FakeRun())
    finally:
        shutil.which = real
    assert "not found on PATH" in str(exc.value)


# --- CLI end to end with stub executables ------------------------------------

def install_stub_bowtie2(tmp_path):
    """Fake bowtie2/bowtie2-build on PATH: build touches index files, align copies grbv.sam."""
    bindir = tmp_path / "bin"; bindir.mkdir()
    build = bindir / "bowtie2-build"
    build.write_text(textwrap.dedent(f"""\
        #!/bin/sh
        base=$(eval echo \\${{$#}})
        for s in 1 2 3 4 rev.1 rev.2; do : > "$base.$s.bt2"; done
    """))
    align = bindir / "bowtie2"
    align.write_text(textwrap.dedent(f"""\
        #!/bin/sh
        # find -S argument
        while [ $# -gt 0 ]; do
          if [ "$1" = "-S" ]; then out="$2"; fi
          shift
        done
        cp "{os.path.join(EX, 'grbv.sam')}" "$out"
        echo "641 reads; of these:" >&2
        echo "  100.00% overall alignment rate" >&2
    """))
    for f in (build, align):
        f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return str(bindir)


def run_cli_with_stubs(tmp_path, argv):
    bindir = install_stub_bowtie2(tmp_path)
    old = os.environ.get("PATH", "")
    os.environ["PATH"] = bindir + os.pathsep + old
    try:
        main(argv)
    finally:
        os.environ["PATH"] = old


def test_cli_maps_unpaired_reads_and_plots(tmp_path):
    fa = tmp_path / "grbv.fa"; fa.write_text(">NC_022002.1\n" + "A" * 3206 + "\n")
    fq = tmp_path / "vine7.fastq"; fq.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "out"
    run_cli_with_stubs(tmp_path, ["-g", os.path.join(EX, "grbv.gff3"), "-y", os.path.join(EX, "grbv.yml"),
                                  "-x", str(fa), "-U", str(fq), "-p", "2",
                                  "-o", str(out), "--name", "t", "-f", "svg"])
    names = sorted(os.listdir(out))
    assert "t.svg" in names and "t.vine7.sam" in names
    assert os.path.isdir(out / "grbv.bt2")                     # index cached under outdir


def test_cli_paired_and_depth_file_together_with_labels(tmp_path):
    fa = tmp_path / "grbv.fa"; fa.write_text(">NC_022002.1\n" + "A" * 3206 + "\n")
    r1 = tmp_path / "s_R1.fq"; r2 = tmp_path / "s_R2.fq"
    for f in (r1, r2): f.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "out"
    run_cli_with_stubs(tmp_path, ["-g", os.path.join(EX, "grbv.gff3"), "-y", os.path.join(EX, "grbv.yml"),
                                  "-d", os.path.join(EX, "grbv.sam"), "-x", str(fa),
                                  "-1", str(r1), "-2", str(r2), "-l", "sim", "paired",
                                  "-o", str(out), "--name", "t", "-f", "svg", "--legend"])
    assert "t.paired.sam" in os.listdir(out)


def test_cli_reads_without_reference_exit(tmp_path):
    fq = tmp_path / "s.fq"; fq.write_text("@r\nA\n+\nI\n")
    with pytest.raises(SystemExit):
        main(["-g", os.path.join(EX, "grbv.gff3"), "-y", os.path.join(EX, "grbv.yml"),
              "-U", str(fq), "-o", str(tmp_path)])


def test_cli_no_depth_and_no_reads_exit(tmp_path):
    with pytest.raises(SystemExit):
        main(["-g", os.path.join(EX, "grbv.gff3"), "-y", os.path.join(EX, "grbv.yml"), "-o", str(tmp_path)])


def test_cli_unbalanced_mates_exit(tmp_path):
    fq = tmp_path / "s.fq"; fq.write_text("@r\nA\n+\nI\n")
    fa = tmp_path / "r.fa"; fa.write_text(">x\nA\n")
    with pytest.raises(SystemExit):
        main(["-g", os.path.join(EX, "grbv.gff3"), "-y", os.path.join(EX, "grbv.yml"),
              "-x", str(fa), "-1", str(fq), "-o", str(tmp_path)])
