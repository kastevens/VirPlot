#!/usr/bin/env python3
"""
align.py — small alignment pipeline that produces a read-depth plot.

Steps:
  1. bowtie2-build   : index the (unindexed) virus FASTA
  2. bowtie2         : align unpaired reads to the index
  3. samtools sort   : sort the alignments into a BAM file
  4. samtools index  : index the sorted BAM
  5. samtools fastq  : FASTQ of the aligned reads (<reads-prefix>_aligned.fastq)
  6. samtools depth  : per-position depth
  7. matplotlib      : depth plot (PNG)

Usage:
  python align.py -U reads.fastq -V virus.fasta [-o outprefix] [-t threads]
"""

import argparse
import os
import shutil
import subprocess
import sys


def run(cmd, **kwargs):
    """Run a command, echoing it first; exit on failure."""
    print("+ " + " ".join(cmd), file=sys.stderr)
    result = subprocess.run(cmd, **kwargs)
    if result.returncode != 0:
        sys.exit(f"ERROR: command failed with exit code {result.returncode}: {cmd[0]}")
    return result


def check_tools(tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        sys.exit("ERROR: required tool(s) not found on PATH: " + ", ".join(missing))


def reads_prefix(fastq_path):
    """Strip .gz and .fastq/.fq from a reads file path, keeping its directory."""
    prefix = fastq_path
    if prefix.endswith(".gz"):
        prefix = prefix[:-3]
    for ext in (".fastq", ".fq"):
        if prefix.endswith(ext):
            prefix = prefix[: -len(ext)]
            break
    return prefix


def main():
    parser = argparse.ArgumentParser(
        description="Align unpaired reads to a viral genome and plot read depth."
    )
    parser.add_argument("-U", "--unpaired", required=True,
                        help="FASTQ file of unpaired reads")
    parser.add_argument("-V", "--virus", required=True,
                        help="FASTA file of the virus genome (will be indexed)")
    parser.add_argument("-o", "--outprefix", default="align",
                        help="prefix for output files (default: align)")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="number of threads for bowtie2/samtools (default: 4)")
    args = parser.parse_args()

    for path in (args.unpaired, args.virus):
        if not os.path.isfile(path):
            sys.exit(f"ERROR: file not found: {path}")

    check_tools(["bowtie2-build", "bowtie2", "samtools"])

    prefix = args.outprefix
    index_base = prefix + "_index"
    sam = prefix + ".sam"
    bam = prefix + ".sorted.bam"
    depth_txt = prefix + ".depth.txt"
    plot_png = prefix + ".depth.png"
    aligned_fastq = reads_prefix(args.unpaired) + "_aligned.fastq"
    threads = str(args.threads)

    # 1. Build the bowtie2 index
    run(["bowtie2-build", "--threads", threads, args.virus, index_base])

    # 2. Align unpaired reads
    run(["bowtie2", "-p", threads, "-x", index_base, "-U", args.unpaired, "-S", sam])

    # 3. Sort alignments into BAM
    run(["samtools", "sort", "-@", threads, "-o", bam, sam])
    os.remove(sam)

    # 4. Index the BAM
    run(["samtools", "index", bam])

    # 5. Write a FASTQ of the reads that aligned (-F 4 excludes unmapped reads)
    run(["samtools", "fastq", "-@", threads, "-F", "4", "-0", aligned_fastq, bam],
        stdout=subprocess.DEVNULL)

    # 6. Compute per-position depth (-a reports positions with zero depth too)
    with open(depth_txt, "w") as fh:
        run(["samtools", "depth", "-a", bam], stdout=fh)

    # 7. Plot
    plot_depth(depth_txt, plot_png, title=os.path.basename(args.virus))
    print(f"Done. Sorted BAM: {bam}\nAligned reads: {aligned_fastq}\n"
          f"Depth table: {depth_txt}\nDepth plot: {plot_png}")


def plot_depth(depth_txt, plot_png, title=""):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("ERROR: matplotlib is required for plotting (pip install matplotlib)")

    positions, depths = [], []
    with open(depth_txt) as fh:
        for line in fh:
            _chrom, pos, depth = line.rstrip("\n").split("\t")[:3]
            positions.append(int(pos))
            depths.append(int(depth))

    if not positions:
        sys.exit("ERROR: no depth data — did any reads align?")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.fill_between(positions, depths, step="mid", alpha=0.4)
    ax.plot(positions, depths, linewidth=0.8)
    ax.set_xlabel("Genome position (bp)")
    ax.set_ylabel("Read depth")
    ax.set_title(f"Read depth across {title}")
    ax.set_xlim(positions[0], positions[-1])
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(plot_png, dpi=150)


if __name__ == "__main__":
    main()
