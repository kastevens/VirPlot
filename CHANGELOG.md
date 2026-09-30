# Changelog

All notable changes to VirPlot will be documented in this file.

Format based on [Keep a Changelog](https://keepachangelog.com/).

---

## [Unreleased]

### Added

- **SAM/BAM input.** `-d/--depth` now accepts SAM or BAM files alongside
  `samtools depth` output; VirPlot computes per-base depth from the
  alignments itself (`virplot.alignments`, standard library only — BAM is
  read through `gzip`; CRAM unsupported). Counting matches `samtools depth -a`
  defaults: unmapped/secondary/QC-fail/duplicate reads skipped, only
  `M`/`=`/`X` bases counted.
- `--ref NAME` to choose the reference from a multi-reference SAM/BAM header
  (default: the GFF sequence id, or the only `@SQ` entry).
- `--min-mapq N` to drop low mapping-quality reads.
- `examples/sample.sam` (synthetic; regenerate with `examples/make_sample_sam.py`).

### Changed

- Internal refactor toward an object model (no user-facing changes; output is
  pixel-identical):
  - New `virplot.models` module with `Feature`, `DepthTrack` and `RNA`
    dataclasses. `RNA.circular` is a placeholder flag for future circular-RNA
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
