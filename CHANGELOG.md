# Changelog

All notable changes to VirPlot will be documented in this file.

Format based on [Keep a Changelog](https://keepachangelog.com/).

---

## [Unreleased]

### Added

- **ICTV-style circular options.** `circular_labels: horizontal` writes every
  arc label level, outside the ring, as `ORF (product)`; a nested arc's label
  steps one line away from the equator per tier so it does not overprint.
  `circular_arcs: on_circle` sets the innermost arcs astride the genome
  circle instead of just outside it. Both default off; `examples/grbv_ictv.yml`
  and `docs/figures/grbv_circular_ictv.svg` show them.
- **Frameshift and readthrough.** A RefSeq CDS split over several rows with
  `exception=ribosomal slippage` is drawn as one box per segment, each
  continuation flipped across the line and labelled `+1 FS` / `−1 FS` (sign
  from the join geometry); a CDS with `transl_except=` is trimmed to its
  extension and drawn abutting the ORF it reads through, behind a bar
  labelled `RT`. Hand-written files can say the same with `Note=+1
  frameshift` / `Note=readthrough`. A frameshift ORF always flips and a
  readthrough extension never does, whatever the overlap and colour rules
  would say. Linear layout only; docs/GFF_GUIDE.md §8.
- **Map reads with bowtie2.** `-x FASTA` plus `-U` (unpaired) or `-1`/`-2`
  (paired) per sample, `-p` threads and `--bowtie2-args`, mirroring bowtie2's
  own options. VirPlot builds and caches the index under `--outdir` (or uses
  one beside the FASTA), runs bowtie2 once per sample, logs its alignment
  summary, keeps the SAM as `<name>.<label>.sam`, and feeds it into the
  depth pipeline; `-d` sources and read samples can be mixed. Warns when the
  FASTA and GFF disagree on a sequence length. Code in `virplot.analysis`;
  bowtie2 is needed only when reads are given.
- **SAM/BAM input.** `-d/--depth` now accepts SAM or BAM files alongside
  `samtools depth` output; VirPlot computes per-base depth from the
  alignments itself (`virplot.alignments`, standard library only — BAM is
  read through `gzip`; CRAM unsupported). Counting matches `samtools depth -a`
  defaults: unmapped/secondary/QC-fail/duplicate reads skipped, only
  `M`/`=`/`X` bases counted.
- `--ref NAME` to choose the reference from a multi-reference SAM/BAM header
  (default: the GFF sequence id, or the only `@SQ` entry).
- `--min-mapq N` to drop low mapping-quality reads.
- **Multiple RNAs.** A GFF with several `region` lines (segmented genomes,
  satellites, subgenomic RNAs) now produces one figure per RNA, named
  `<name>.<seqid>.<format>`, with per-RNA CSV reports. Depth sources are
  matched to each RNA by sequence id, so one SAM/BAM or multi-sequence depth
  file serves all of them. Figures share a depth y-limit and have widths
  proportional to RNA length so they compare honestly side by side;
  `--free-y` / `--equal-width` opt out, `--rnas` selects a subset.
- `examples/sample.sam` (synthetic, two references; regenerate with
  `examples/make_sample_sam.py`) and `examples/sample_multi.gff3`.

- **Circular genomes.** `Is_circular=true` on a GFF `region` line (as NCBI
  RefSeq writes) marks an RNA circular; `--topology circular|linear` overrides
  it. For a circular genome, alignments that run off the end of the reference
  wrap round to position 1 and positions past the end are taken modulo the
  length, so depth is continuous across the origin for a BAM produced by
  padding the reference and wrapping coordinates back. The annotation track
  drops the 5'/3' marks and shows the backbone continuing past both edges.
- **Circular layout.** `--layout circular` draws a genome as a circle on a
  single polar axes: feature arcs in a two-lane outer ring, depth as a filled
  band inside it, position 1 at the top running clockwise. `--layout auto`
  picks it for genomes marked circular; the default stays `linear`. A feature
  crossing the origin is drawn as the two arcs it occupies
  (`RNA.feature_spans`). `--yscale symlog` is ignored in this layout.
- `examples/byv.gff3` / `byv.sam` — beet yellows virus with real NC_001598.1
  CDS coordinates and ICTV ORF names, the fixture the drawing conventions
  are checked against; `examples/make_reads.py`, a generic flat-coverage
  read simulator. docs/ictv_drawing_conventions.md §7 records the
  comparison with the published figures.
- `examples/grbv.gff3`, `examples/grbv.sam`, `examples/grbv.yml` — grapevine
  red blotch virus (RefSeq NC_022002.1, 3206 nt circular): real RefSeq
  annotation with simulated, origin-spanning reads
  (`examples/make_grbv_sam.py`).

- **ICTV drawing conventions** (docs/ictv_drawing_conventions.md, stages
  0–3 and 6). New `virplot.layout` places features by rule: one-strand genomes
  flip a box across the line only when it overlaps its upstream neighbour;
  two-strand genomes keep + above / − below as arrows and tier same-side
  overlaps outward; circular genomes nest overlapping arcs inward, largest
  outermost, with arrowheads clockwise for + and anticlockwise for − and a
  stem-loop icon at the origin. Colour now encodes function: products not in
  `color_mapping` are classified by name onto the ICTV palette ("putative" →
  lighter tint). spec.yml gains `overlap_mode`, `end_5_label`/`end_3_label`
  (`VPg` draws the grey oval) and `function_palette`. GFF `gene=` is drawn
  outside the glyph as the ORF name. `--title` with no YAML title uses
  `name (length nts)`. Section 6 of the document records what changed and
  what the spec still needs.

- `docs/GFF_GUIDE.md`: how to write a GFF3 for VirPlot — what is read, `product`/`gene`
  naming, strand, `Is_circular` and origin-crossing features, frameshift convention,
  and the function words behind the default colours.
- `docs/ARCHITECTURE.md`: module map, pipeline and class diagrams (Mermaid,
  rendered by GitHub), a step-by-step run, the behaviours that are easy to
  get wrong, extension points and the test suites. `docs/figures/`: example
  renders linked from the README.

### Changed (output)

- Feature glyphs are drawn as flat colour with no outline, as in the
  original VirPlot figure and the ICTV maps. The `--no-border` option is
  removed (there is nothing to turn off). To keep adjacent same-colour boxes
  apart, the one-strand flip rule now also flips a box whose upstream
  neighbour is the same colour and within 1 % of the genome length — the
  choice the ICTV BYV figure makes for CPm/CP. With `examples/byv.yml` (the
  figure's palette) BYV now matches Closteroviridae Fig 2 on every ORF.
- Circular arcs and their arrowheads are one path each, removing the
  hairline seam between body and head.
- Depth trace restyled to match the original figure, whose values were
  measured from its pixels: a 0.3 pt black outline at alpha 0.8 over a fill
  at alpha 0.9 in `depth_line_color`, in both layouts, instead of a heavier
  coloured line over a pale fill; the legend swatch now shows the fill colour.
  The original figure's depth palette, recovered from the same measurement
  (`#567eb0`, `#55b9d9`, `#85dbec`, bottom layer first), is now the built-in
  default for `stacked_area_colors`, and its top layer `#85dbec` the default
  `depth_line_color`, so a single sample looks like sample 1 of that figure.
  The example specs use the same values.
- SVG and PDF output is reproducible: no embedded timestamp, and element
  ids are hashed from the output name, so regenerating an unchanged figure
  gives a byte-identical file. `docs/make_figures.sh` regenerates the
  documentation figures.

- Linear annotation track: boxes no longer alternate above/below by index;
  they follow the ICTV flip/tier rules. Genomes without overlapping ORFs now
  draw every box above the line. Products with no explicit colour that match
  a function class are coloured by that class instead of grey.

### Fixed

- `--report` no longer fails when the output directory does not exist yet.

### Changed

- Internal refactor toward an object model (no user-facing changes; output is
  pixel-identical):
  - New `virplot.models` module with `Feature`, `DepthTrack` and `RNA`
    dataclasses; `parse_gff_rnas()` returns one `RNA` per GFF region. `RNA.circular` is a placeholder flag for future circular-RNA
    support.
  - `parse_gff` now returns `Feature` objects instead of dicts.
  - `plotting.plot()` replaced by a `Plotter` base class and `LinearPlotter`;
    `cli.main()` builds one `RNA` and calls `LinearPlotter(...).render(rna, ...)`.
  - Smoothing/normalisation is computed once (`Plotter._prepared_tracks`)
    instead of separately for the depth panel and the y-limit.
- Added `docs/class_layout.svg` showing the current classes and the seams left
  open for circular RNAs and segmented (multi-RNA) genomes.

## [2.0.0] — 2026-05-03

### Added

- Installable Python package (`pip install .`) with `virplot` entry point
- PEP 621 `pyproject.toml` build metadata with setuptools backend
- `python -m virplot` support via `__main__.py`
- `-V` / `--version` flag
- `-f` / `--format {svg,pdf,png}` flag replacing `--Opdf` / `--Opng`
- `-v` / `--verbose` flag for debug-level logging
- Input validation: depth file format, tab delimiters, cross-file position counts
- YAML schema warnings for unknown keys, info for missing keys using defaults
- `Settings` dataclass for typed configuration
- Type hints on all function signatures
- Logging via `logging` module with `[VirPlot]` prefix (INFO omits level tag)
- Named constants for layout parameters in `plotting.py`
- Unit tests (26 tests across parsers, analysis, settings, CLI formatter)
- Synthetic sample data in `examples/`
- `examples/spec.yml` template configuration
- `CHANGELOG.md`

### Changed

- Restructured from monolithic `bin/virplot` script to `src/virplot/` package
- `parse_depth()` now takes `seq_len` and returns `(np.ndarray, n_entries)` directly, without intermediate dict
- Depth parsing enforces tab-delimited 3-column format with line-number errors
- `load_settings()` returns a `Settings` dataclass instead of a tuple
- Log output cleaned: INFO messages omit redundant level prefix
- `figsize` and `height_ratios` extracted as module-level constants
- `CITATION.cff` updated to v2.0.0
- `README.md` rewritten for new package structure and pip workflow
- `.gitignore` cleaned up, removed irrelevant framework sections

### Removed

- `--Opdf` / `--Opng` flags (replaced by `--format`)
- `bin/virplot` monolithic script
- `bin/env.yml` (superseded by pip workflow)
- `bin/spec.yml` (moved to `examples/spec.yml`)
- `build_depth_array()` function (merged into `parse_depth()`)
- `print()` statements (replaced by `logging`)

---

## [1.0.0] — 2025-08-24

### Added

- Initial release
- Combined annotation + depth plot from GFF3 and samtools depth files
- YAML-based color and style configuration
- Stacked area chart for multiple depth files
- Coverage interval and gap analysis with CSV reports
- Supplementary scripts: `depth_filter.py`, `depth_merger.py`
