# Changelog

All notable changes to VirPlot will be documented in this file.

Format based on [Keep a Changelog](https://keepachangelog.com/).

---

## [Unreleased]

### Added

- **A real-read example: GLRaV-13 isolate H8881** (`examples/glrav13.gff3`,
  `glrav13.fa`, `glrav13.bam`, `glrav13.yml`). Every other example ships
  simulated reads; this one carries 212,285 real alignments against the real
  reference, mean depth 1,173x with 100% of the genome at >=1x, so the depth
  trace shows genuine coverage structure instead of a flat band.

  The BAM holds only the reads VirPlot counts: the source alignment had 23.0M
  records, of which 212,285 survive the unmapped / secondary / QC-fail /
  duplicate filter, taking it from 899 MB to 4.7 MB. The depth it yields is
  identical to the depth from the unfiltered file. Its `@SQ` name and the
  FASTA header were both changed from the assembler's contig name
  (`GVLRaV13_H8881NODE_4_length_...`) to the GFF3's `GLRaV13_H8881`, so the
  example plots with no `--ref` and no name-mismatch warning — the §7 problem,
  avoided rather than demonstrated.

## [2.2.0] — 2026-10-02

### Changed

- **A circular genome is now drawn as a circle by default.** `--layout`
  defaults to `auto` rather than `linear`, so a GFF whose region line says
  `Is_circular=true` produces the circular figure the ICTV conventions
  describe (§4) without being asked. spec.yml gains a `layout` key
  (`auto|linear|circular`) so a per-example default lives beside the data it
  styles; precedence is `--layout`, then the YAML key, then `auto`.

  `--layout linear` is unchanged and now matters more: it is the way to get
  the linear track of a circular genome. Deleting `Is_circular=true` is not,
  because that attribute is a fact about the molecule — it drives depth
  wrapping at the origin and the splitting of origin-crossing features — so
  removing it would quietly change the data handling too.

  Because the circular layout cannot carry the `±1 FS` / `RT` marks, it now
  warns when a genome reaching it has frameshifted or read-through ORFs,
  rather than dropping the annotation silently.

  `docs/make_figures.sh` passes `--layout linear` explicitly for the
  `grbv_linear` figure, which previously relied on the old default. No shipped
  figure changes.

### Added

- **Stacked figures for segmented genomes, by composing.** `bin/stack_figures.py`
  joins a multipartite genome's per-segment figures into the one stacked figure
  the ICTV conventions describe — largest first, shared scale, left aligned at
  the 5' end (Geminiviridae DNA-A/DNA-B, LIYV RNA-1/RNA-2, TSWV L/M/S, nanovirus
  components). Supporting it: segments now render **largest first** by default
  (`--rna-order length|gff`), and `--bare-x` drops the repeated x-axis label and
  tick numbers from every panel but the bottom one.

  There is deliberately no stacked-figure renderer. Within a run the figure
  width is proportional to genome length and matplotlib's margins are
  fractional, so the separate panels already share an exact x-scale — measured
  across TSWV's 3x length spread, 0.0000% drift — which makes stacking a
  document operation rather than a re-plot. The composer also **aligns the
  panels' plot areas**, correcting the small offset (0.63 pt on LIYV) that
  appears when panels carry y-tick labels of different widths and are cropped
  to a tight bounding box; hand pasting cannot fix that. `tests/test_stacking.py`
  guards the shared-scale property, since the decision holds only while it does.
  See `docs/ictv_drawing_conventions.md` §6 G.

- `docs/ARCHITECTURE.md` §3 class diagram simplified to the objects and
  modules that matter, with the design patterns the code happens to use
  marked (Value Object, Aggregate, Builder, Factory + Iterator, Strategy,
  Rule table); `docs/class_layout.svg` regenerated from it.
- **Subgenomic RNA rows.** The nested sgRNA sets that plus-strand RNA plant
  viruses use for their 3' ORFs are drawn as shorter lines beneath the genome
  line, 5' aligned to their start with an arrowhead there, longest row first
  and each labelled — the ICTV convention. They are read from a transcript row
  (`mRNA`, `transcript`, `ncRNA`, `misc_RNA`, `primary_transcript`,
  `sequence_feature`, `misc_feature`) whose `Note=` says `sgRNA` or
  `subgenomic`, with `gene=` naming what it expresses; the marker is required
  rather than inferred from the type, because RefSeq writes real `mRNA` rows
  (spliced mastrevirus transcripts) that are not sgRNAs. `start` is the 5'
  terminus; an end at or past the genome length means the 3' end, covering the
  usual 3'-coterminal set, while an earlier end is kept for the 5'-proximal
  sgRNAs closteroviruses also make. Row colour is `sgrna_color` in spec.yml.
  `examples/byv.gff3` carries the full closterovirus ladder.

  These rows are annotation and never a depth track: a plant virus sgRNA is
  co-linear with the genome and carries no leader junction, so its reads
  cannot be told from genomic reads and per-sgRNA coverage is not recoverable
  from short reads. The set does leave a step in the aggregate depth at each
  5' end, which is why the ladder shares the depth track's x-axis — a
  closterovirus covered 1x over ORF1a and 150x over CP is showing its
  expression strategy, not a failed assembly.

  Drawn in the linear layout only, deliberately: the ICTV puts no transcript
  rows on circular genomes (geminivirus transcription is bidirectional from
  the IR, nanovirus components carry one ORF each), so `--layout circular`
  warns and skips. See `docs/GFF_GUIDE.md` §11 and
  `docs/ictv_drawing_conventions.md` §6 K.

- **Non-coding features: UTRs, intergenic regions, stem-loops.** GFF
  `five_prime_UTR` / `three_prime_UTR` / `stem_loop` / `origin_of_replication`
  rows are always read; GenBank's catch-all `misc_feature` (NCBI GFF3
  `sequence_feature`), `regulatory_region` and `repeat_region` rows only when
  their name or `Note=` says intergenic / common region / IR / UTR /
  stem-loop, so RefSeq motif annotations stay off the line. A region is a
  grey bar astride the genome line (linear) or a grey arc astride the circle
  (circular, the ICTV begomovirus `CRA`), under the ORFs and outside the
  layout rules; a stem-loop is the hairpin icon at its own position, which on
  a circle replaces the default icon at position 1. A region holding an
  unnamed stem-loop is named once, over the hairpin. Labels: `Name=` first,
  then `product=`/`gene=`/`Note=`, with defaults for the dedicated types.
  `examples/grbv.gff3` now carries its origin-spanning intergenic region and
  hairpin, `examples/tgmv.gff3` its common regions (`CRA`/`CRB`) and
  hairpins — all located from the sequences, since RefSeq annotates none —
  and `examples/fbnyv.gff3` RefSeq's own `stem_loop` rows on every
  component; docs/GFF_GUIDE.md §10.
- Inside labels switch to white on a really dark fill (`color_mapping`
  navy, the non-coding grey); every ICTV palette colour keeps black text.
- **Polyprotein domain dividers.** RefSeq's `mature_protein_region_of_CDS`
  rows (GenBank `mat_peptide`; `Parent=` names the CDS, else containment)
  become `Domain` segments of their CDS. The box is drawn as one glyph split
  by thin vertical lines, each segment in its own colour (same rule as any
  product) and labelled with its domain name — inside when it fits at the
  figure's size, otherwise just outside the box, staggered onto a further
  row when names would overprint. The polyprotein's own `product` is not
  written. Works in both layouts. `examples/pvy.gff3` / `pvy.sam` /
  `pvy.yml` — potato virus Y with its ten mature proteins and PIPO, styled
  to ICTV Potyviridae Fig. 2 — and `docs/figures/pvy_linear.svg`;
  docs/GFF_GUIDE.md §9.
- Frameshift labels carry whatever shift the `Note=` states (`+2 FS` for
  PIPO), not only ±1.
- **ICTV reference examples.** Seven curated RefSeq genomes, each labelled
  and coloured as in its ICTV 9th Report family figure, with a style file and
  simulated reads: BYDV-PAV and PLRV (−1 frameshift, readthrough, VPg), LIYV
  (two RNAs, BAM input), TSWV (three ambisense segments), TGMV (bipartite
  begomovirus), FBNYV (eight-component nanovirus) and PVX (triple gene
  block). Each GFF header lists its accessions and every change from RefSeq.
  `examples/make_example_sam.py` simulates reads for every molecule in a GFF.
- `examples/grbv.gff3` labels curated to the ICTV convention (`gene=V1;
  product=CP`, `gene=C1;product=RepA`, …; RefSeq's `Vn protein` names kept in
  `Note=`), with a header comment recording the curation. `grbv.yml` no longer
  needs a `color_mapping`; the default palette colours the curated products.
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
  original VirPlot figure and the ICTV maps: `--no-border` is now the
  default (and still accepted, so existing commands keep working), and the
  new `--border` restores the black outline. To keep adjacent same-colour boxes
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

- GFF3 rows that precede their sequence's `region` line are no longer dropped
  (with a misleading "no region line" warning); features are attached after
  the whole file is read.
- GFF3 attribute values are percent-decoded (`%2C` → `,`, `%3B` → `;`, as
  NCBI writes inside `product=` and `Note=`), so labels read naturally and the
  function-word classifier sees real text.
- A named or 3-digit colour in `function_palette` (`orange`, `#abc`) no
  longer crashes the "putative = lighter tint" path; `lighter()` accepts any
  matplotlib colour.
- `color_mapping:` left empty in the YAML no longer raises `TypeError`.
- A hand-written join (two rows sharing an `ID`) with `Note=… frameshift` on
  both rows marks only the continuation as the frameshift, not the first
  segment too.
- Join segments that are in frame with each other get no sign (bare `FS`)
  instead of `−0 FS`.
- A frameshift continuation now flips relative to the ORF it continues
  (found at the junction, ±2 nt), not to whichever ORF sorted before it by
  start; a small ORF nested inside ORF1a no longer makes ORF1b land on
  ORF1a's side.
- Origin-crossing coordinates (`start > end`, or `end` past the length) are
  split into their two pieces whether or not the molecule is marked circular,
  so `--topology linear` on a circular genome, or a GFF that forgot
  `Is_circular`, draws the feature where the coordinates say instead of as
  an off-axis box; a `start` past the length (padded reference) is reduced
  modulo the length; the parser warns when a non-circular region has such
  a feature.
- Two `samtools depth` files of different row counts no longer abort for a
  circular genome (one may have been made against a padded reference).
- `Feature.length` counts bases inclusively (`Feature(1, 3).length == 3`).
- Malformed GFF or depth text raises `virplot.parsers.ParseError` (caught by
  the CLI) instead of calling `sys.exit` from library code.
- The linear depth outline is drawn stepped when the fill is (two or more
  tracks), so it sits on the fill's edge.
- All-zero depth no longer trips matplotlib's identical-ylim warning; an
  `RNA` with no depth track plots as flat zero instead of raising; genomes
  under eight bases get a tick step of one instead of a division by zero.
- `CircularPlotter`'s per-figure radii start at sane defaults, so its label
  helpers can be called in any order.
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
