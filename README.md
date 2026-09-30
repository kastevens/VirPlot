# VirPlot

This tool generates **SVG, PDF, or PNG plots** that combine viral genome feature annotations from a GFF3 file and sequencing depth from a `samtools depth` file or directly from SAM/BAM alignments.

[![RNA-seq read depth across the Beet yellows virus genome](https://i.imgur.com/bDxrx1I.png)](https://i.imgur.com/bDxrx1I.png)

Developed and maintained by Haoran (Henry) Li for [Foundation Plant Services](https://fps.ucdavis.edu/index.cfm) at [UC Davis](https://www.ucdavis.edu/).

Built in Python using `matplotlib`, `pyyaml`, and `numpy`.

---

## Overview

In Next Generation Sequencing (NGS), read **depth** provides a measure of how confidently a region of the genome was sequenced. When depth is visualized alongside gene or product annotations, it allows researchers to:

- Identify under or over-sequenced regions.
- Assess the relationship between functional elements and read depth.
- Validate annotations visually.
- Visualize contiguous coverage intervals and reveal potential breaks or low-support regions.

This tool generates a vertically stacked plot:

- **Top**: Gene product annotations (e.g. RdRp, CP, MP, etc.) as labeled colored boxes on a genome line.
- **Bottom**: Sequencing depth plot over the same genomic range, with optional shading of below-threshold gaps.

Both plots are **aligned on a shared x-axis** and exported as vector-format graphics suitable for publications.

---

## Installation

Requires Python 3.9+.

```bash
pip install .
```

---

## Quick Start

A sample GFF3, depth file, SAM file, and YAML configuration are included in `examples/`:

```bash
# Basic plot
virplot -g examples/sample.gff3 -d examples/sample.dep -y examples/spec.yml
```

```bash
# PNG with smoothing, gap shading, and grid
virplot -g examples/sample.gff3 -d examples/sample.dep -y examples/spec.yml \
  --smooth --shade-breaks --grid -f png
```

```bash
# PDF with custom name and output directory
virplot -g genome.gff3 -d depth.dep -y spec.yml \
  -f pdf --name my_virus -o results/
```

To create your own configuration, copy the template and edit the color mappings:

```bash
cp examples/spec.yml my_project.yml
```

---

## Inputs

Three input files are required:

1. **GFF3 file** — genome annotations containing CDS features and a `region` entry for sequence length
2. **Depth source** — either
   - tab-delimited output from `samtools depth -a`, **or**
   - a **SAM or BAM** file; VirPlot computes per-base depth from the alignments itself (no `samtools`/`pysam` needed)
3. **YAML file** — color mapping and style settings (see [Customization](#customization))

If multiple depth sources are given, VirPlot combines them into a *stacked area chart* showing each sample's depth contribution under a combined depth line. Depth files and SAM/BAM files can be mixed.

### Depth from SAM/BAM

```bash
virplot -g examples/sample.gff3 -d examples/sample.sam -y examples/spec.yml --smooth --shade-breaks
```

The file type is detected from the extension (`.sam`, `.bam`) or, failing that, from the content. Counting follows `samtools depth -a` defaults so results are directly comparable:

- unmapped, secondary, QC-fail and duplicate reads are skipped;
- only aligned bases (CIGAR `M`, `=`, `X`) count — deletions and `N` skips do not;
- no mapping-quality filter unless `--min-mapq` is given.

The reference to plot is taken from the GFF `region` line's sequence id. If the SAM/BAM header has a single `@SQ` entry that one is used regardless; if it has several and none matches, pass `--ref NAME`. BAM is read via the standard library's gzip module; CRAM is not supported.

---

## Usage

```txt
virplot [-h] [-V] -g GFF -d DEPTH [DEPTH ...] [-l LABELS [LABELS ...]]
        [--ref REF] [--min-mapq MIN_MAPQ] -y YAML [-o OUTDIR] [-n] [--grid] [--smooth]
        [--yscale {linear,symlog}] [--linthresh LINTHRESH]
        [--name NAME] [--no-label] [--no-border]
        [-t THRESHOLDS [THRESHOLDS ...]] [-r] [--shade-breaks]
        [--legend] [--title] [-f {svg,pdf,png}] [-v]
```

### Options

| Flag                 | Description                                                    |
| -------------------- | -------------------------------------------------------------- |
| `-g`, `--gff`        | Path to GFF3 annotation file                                   |
| `-d`, `--depth`      | One or more depth sources — `samtools depth` files or SAM/BAM (stacked if multiple) |
| `-l`, `--labels`     | Label(s) for each depth source (same order as `--depth`)       |
| `--ref`              | Reference name to use from SAM/BAM headers (default: GFF sequence id, or the only reference) |
| `--min-mapq`         | Skip SAM/BAM reads with MAPQ below this (default: 0)           |
| `-y`, `--yaml`       | YAML file for color mapping and other specs                    |
| `-o`, `--outdir`     | Output directory (default: `.`)                                |
| `-f`, `--format`     | Output format: `svg`, `pdf`, or `png` (default: `svg`)         |
| `--name`             | Base name for output file (default: `virplot`)                 |
| `-n`, `--normalize`  | Normalize depth values to max=1                                |
| `--smooth`           | Smooth depth plot using moving average                         |
| `--grid`             | Enable background grid on depth plot                           |
| `--yscale`           | Y-axis scale: `linear` or `symlog` (default: `linear`)         |
| `--linthresh`        | Symlog linear threshold around 0 (default: `10.0`)             |
| `--no-label`         | Hide feature labels                                            |
| `--no-border`        | Remove borders around feature rectangles                       |
| `-t`, `--thresholds` | Coverage thresholds for interval/gap analysis (default: `1 5`) |
| `-r`, `--report`     | Write CSV reports of intervals and gaps per threshold          |
| `--shade-breaks`     | Shade coverage gaps on the depth plot                          |
| `--legend`           | Show depth plot legend                                         |
| `--title`            | Show title from YAML                                           |
| `-v`, `--verbose`    | Enable debug-level logging                                     |
| `-V`, `--version`    | Show version and exit                                          |

---

## Customization

All visual elements are customizable via the YAML file:

| Key                   | Description                           |
| --------------------- | ------------------------------------- |
| `color_mapping`       | Maps product names to hex colors      |
| `default_color`       | Color for unmapped products           |
| `shade_color`         | Color for shaded gap regions          |
| `depth_line_color`    | Line color for depth plot             |
| `annotation_fontsize` | Font size for feature labels          |
| `stacked_area_colors` | Color palette for stacked area chart  |
| `legend_location`     | Legend position (e.g. `"upper left"`) |
| `title`               | Plot title text                       |

See `examples/spec.yml` for a complete template.

---

## Supplementary Scripts

Two standalone scripts in `bin/` assist with depth file preparation:

- **`depth_filter.py`** — Filter depth entries by sequence header:

  ```bash
  python3 bin/depth_filter.py -i input.dep -o filtered.dep --seq SEQ_NAME
  ```

- **`depth_merger.py`** — Merge depth counts across multiple files:

  ```bash
  python3 bin/depth_merger.py -i sample1.dep sample2.dep -o combined.dep
  ```

---

## License

This project is released under the MIT License.
