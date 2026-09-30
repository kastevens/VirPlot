# Genome-map drawing conventions (ICTV 9th Report style)

Reference: King, Adams, Carstens & Lefkowitz (eds), *Virus Taxonomy: Ninth Report
of the ICTV* (2012), family figures, e.g. Closteroviridae Figs 2–4, Geminiviridae
Fig 5, Tombusviridae Fig 2, Luteoviridae Fig 2, Betaflexiviridae Fig 2,
Circoviridae Fig 2. All figures are online at https://ictv.global/report_9th.

VirPlot's annotation track should follow these conventions so that output looks
like the figures virologists already read. The rules below were extracted by
inspecting the figures directly; they are a target spec for the plotting code.

---

## 1. Shared rules (all modes)

| Element | Convention |
|---|---|
| Scale | Drawn to nucleotide length. Title: virus name, acronym, length, e.g. `BYV (15,468 nts)`. |
| ORF glyph | Rectangle (linear) or thick arc with arrowhead (circular). |
| Labels | ORF number outside the box (`ORF1a`, `AC1`); product/domain name inside (`RdRp`, `CP`, `39K`). |
| Colour | Encodes **predicted function / homology**, never reading frame or strand. Homologous ORFs share a colour across all panels. |
| Palette | replicase (Mtr/Hel/RdRp, Rep) yellow-orange · CP/CPm magenta-pink · MP teal/blue · HSP70h green · silencing suppressor / small 3′ ORFs purple · unknown grey · non-coding features (IR, UTR, stem-loop) grey or black · "putative" function = lighter tint. |
| Readthrough | One continuous box; vertical bar + small arrow at the stop codon, labelled `RT`. |
| Frameshift | Step between two abutting boxes, labelled `+1 FS` / `−1 FS`. |
| Polyprotein | Domains written as text inside one box, separated by thin vertical lines. |
| sgRNAs | Shorter lines stacked beneath the genome, 5′ aligned to their start, labelled (`sgRNA1`). |
| Segmented genomes | One row (linear) or one circle (circular) per segment, stacked largest first, each labelled (`RNA-1 (3,840 nts)`). |

---

## 2. Linear genomes

Baseline: horizontal black line, 5′ left, 3′ right. 5′ end marked `5′m⁷G` (cap)
or a grey oval `VPg`; 3′ end `3′OH` or `A(n)`.

### Mode L1 — all ORFs on one strand (typical +ssRNA)

- ORF boxes sit **on** the line; reading direction is implicitly left→right.
- First (5′-most) ORF is drawn above the line.
- An ORF that overlaps its upstream neighbour flips to the opposite side.
  Non-overlapping neighbours stay on the same side. Overlapping runs therefore
  alternate above/below (BYV, CTV, ASGV, luteovirids, tombusvirids).
- An ORF that is the **same colour** as its upstream neighbour and within
  ~1 % of the genome length of it also flips, so two flat boxes do not read
  as one (BYV CPm/CP, 70 nt apart, both magenta). Added after checking the
  figures — see §7.
- Reading frame is **not** encoded by vertical position (Caulimoviridae is the
  one family that tiers by frame).

### Mode L2 — ORFs on both strands (ambisense RNA; dsDNA)

- Genomic / + sense ORFs **above** the line with rightward arrows.
- Antigenomic / − sense ORFs **below** the line with leftward arrows.
- Side is fixed by strand, so overlap is resolved by **tiering** further from the
  line within the same strand — never by flipping across it.
- Ambisense segments show the intergenic hairpin between the two oppositely
  oriented ORFs.

---

## 3. Circular genomes

Baseline: a circle with the origin of replication / intergenic region
(stem-loop icon) at 12 o'clock; IR drawn as a grey or black arc segment.
Virus name and length in the centre.

### Mode C1 — all ORFs on one strand

- ORFs are arcs just outside the circle, running **clockwise** from the origin,
  arrowheads showing direction.
- Overlapping ORFs nest **inward** (smaller radius), keeping the same direction.

### Mode C2 — ORFs on both strands (geminiviruses, nanoviruses, circoviruses)

- Virion-sense ORFs run **clockwise** from the origin down the **right** side
  (V1/AV1, V2/AV2; `rep` in PCV).
- Complementary-sense ORFs run **anticlockwise** from the origin down the
  **left** side (C1/AC1 Rep, C2 TrAP, C3 REn, C4; `cap` in PCV).
- Arrowheads point in the direction of translation on each strand.
- Overlapping ORFs on the same strand step inward (AC4 inside AC1), so
  **radius encodes nesting** and **half of the circle encodes strand**.
- Bipartite genomes: two circles side by side (DNA-A, DNA-B), common region
  marked `CRA` / `CRB` at the top of each.

---

## 4. Decision logic

```
circular?       -> C modes, else L modes
both strands?   -> mode 2: side/half fixed by strand; overlap -> tier / nest
                   else mode 1: linear flips side on overlap; circular nests inward
```

Both flags come from the data, not from taxonomy:

- `RNA.circular` — already set from the GFF3 `Is_circular=true` region attribute
  or the `--topology` CLI flag.
- `two_strand` — `any(feature.strand == "-")` over the ORFs of that RNA.

---

## 5. Staged implementation plan

Each stage is independently shippable and testable. File names refer to the
`refactor/oo-models` layout.

| Stage | Scope | Touches |
|---|---|---|
| 0 | Document conventions (this file). Add `two_strand` property to the `RNA` model; add `overlap_mode: auto \| flip \| tier` to the spec.yml schema. | `docs/`, `models.py`, `settings.py`, `examples/spec.yml` |
| 1 | **L1 flip-on-overlap.** Replace current box placement with: first ORF above, flip side only when overlapping upstream neighbour. Regression-test against the BYV/sample example. | `plotting.py`, `tests/` |
| 2 | **Function palette.** Map product names (RdRp, CP, MP, HSP70h, …) to the ICTV palette by default; allow override in spec.yml. | `settings.py`, `plotting.py` |
| 3 | **L2 two-strand.** Carry strand on `Feature`; + above/→, − below/←; same-strand overlap tiers outward. `overlap_mode: auto` picks flip vs tier from `two_strand`. | `models.py`, `parsers.py`, `plotting.py` |
| 4 | **Expression features.** RT bar, FS step, polyprotein domain dividers, sgRNA rows, 5′/3′ end glyphs — driven by GFF3 attributes or spec.yml. | `parsers.py`, `plotting.py` |
| 5 | **Segmented genomes.** Multiple `##sequence-region` entries → stacked rows sharing a scale, largest first (`sample_multi.gff3`). | `plotting.py`, `cli.py` |
| 6 | **Circular C1/C2.** Polar annotation track; origin at top; V-sense clockwise right, C-sense anticlockwise left; inward nesting. Depth track drawn as a ring or as an unrolled linear panel beneath (GRBV example). | new `circular.py`, `plotting.py` |

Stages 1–3 change existing linear output and should each bump the minor
version; stages 4–6 add features.

---

## 6. Implementation status (refactor/oo-models)

What the code does today, stage by stage, and where implementing the spec
showed the spec itself needs a change. Rule numbers refer to the sections above.

| Stage | Status | Where |
|---|---|---|
| 0 | Done. `RNA.two_strand`; `overlap_mode`, `end_5_label`, `end_3_label`, `function_palette` in spec.yml. | `models.py`, `settings.py` |
| 1 | Done. Flip-on-overlap with a collision fallback (see A below). | `layout.py` |
| 2 | Done. Keyword classifier → ICTV palette; `color_mapping` still wins; "putative" → lighter tint. | `settings.py` |
| 3 | Done. `Feature.strand`/`gene` carried; + above →, − below ←; same-side overlap tiers outward; `auto` picks flip/tier from `two_strand`. | `layout.py`, `plotting.py` |
| 4 | **Partly.** 5′/3′ end labels and the VPg oval; ORF name outside the glyph from `gene=`; default title `name (length nts)`; **frameshift step + `±1 FS` label and readthrough bar + `RT`**, read from RefSeq's `exception=ribosomal slippage` / `transl_except=` or from `Note=` (linear layout). Polyprotein dividers, sgRNA rows: **not done** (see F). | `plotting.py`, `parsers.py`, `layout.py` |
| 5 | **Done differently.** One figure per RNA with a shared y-limit and length-proportional width, not stacked rows (see G). | `cli.py`, `plotting.py` |
| 6 | Done. `CircularPlotter`: origin at 12 o'clock with stem-loop icon, arcs just outside the circle with arrowheads, clockwise for + and anticlockwise for −, overlap nests inward; depth as an inner ring. `--layout circular|auto`. | `plotting.py` |

### Changes the document needs

**A′. §2 L1 — the flip rule needs a legibility clause.** Boxes carry no
outline, so two adjacent same-colour boxes merge. The figures flip such a pair
even without overlap (BYV CP below CPm). Implemented as: flip when the
upstream neighbour resolves to the same colour and the gap is under
`NEAR_GAP_FRACTION` (1 %) of the genome. With the figure's palette
(`examples/byv.yml`) this makes BYV agree with Closteroviridae Fig 2 on all
nine ORFs — including ORF1a→RdRp, which flips because both are replicase
yellow and 2 nt apart, not because the code knows about the frameshift.

**A. §2 L1 — the flip rule can overprint.** "Flip when overlapping the upstream
neighbour" only looks one feature back. With a long upstream ORF (a polyprotein
spanning several small downstream ORFs) the flipped side can already be
occupied. The code keeps the rule and adds a fallback: if the wanted side is
taken where the box would land, try the other side; if both are taken, tier
outward on the wanted side. The document should state this so the rule is
complete.

**B. §3 — nesting order must be specified.** "Overlapping ORFs nest inward" is
order-dependent: placed by start coordinate, AC4 (lower start) would take the
outer lane and AC1 nest inside it, the opposite of the figures. The code places
arcs *largest first*, so the small ORF nests inside the big one it overlaps
regardless of coordinates. Suggest adding "largest arc outermost" to C1/C2.

**C. §3 C2 — "half of the circle encodes strand" is a consequence, not a rule.**
With the origin at 12 o'clock and coordinates clockwise, geminivirus V ORFs
fall on the right and C ORFs on the left because of where they sit in the
sequence. The code encodes strand only through arrow direction and nests by
positional overlap irrespective of strand (so opposite-strand ORFs that do
overlap positionally — not the geminivirus case — cannot collide). Suggest
rewording C2 to: direction by strand; nesting by overlap; halves follow.

**D. §1 Labels — say which GFF attribute feeds which label.** ORF number
outside / product inside needs two names per feature. The code uses
`product=` inside and `gene=` outside when present. NCBI RefSeq GFF3 for GRBV
carries neither an ORF name nor a function — only `locus_tag=N761_gp1` and
`product=V1 protein` — so the outside label is empty and function colouring
needs a per-product `color_mapping` (as `examples/grbv.yml` does). Worth
deciding whether `locus_tag` or `Name` should be accepted as the outside label.

**E. §1 Palette — RefSeq product names are often function-free.** The keyword
classifier covers RdRp/Rep/CP/MP/HSP70/p2x-style names; anything else falls to
grey unless mapped. That is the right failure mode, but the document's "map
product names by default" will silently leave many RefSeq genomes grey.

**F. §1 Expression features need an encoding — two now have one.** Frameshift
and readthrough are read from what RefSeq already writes (`exception=ribosomal
slippage` on the rows of a joined CDS; `transl_except=` on a readthrough
product) and from VirPlot's `Note=+1 frameshift` / `Note=readthrough`
convention (docs/GFF_GUIDE.md §8). The frameshift sign is derived from the
join geometry (skip one base = +1, re-read one = −1). Polyprotein domains
(`mat_peptide` rows) and sgRNAs still need an encoding decision, and the
parser still reads only `CDS` rows, so IR/UTR/stem-loop arcs (§3 baseline)
need `parse_gff_rnas` to accept other types.

**G. §1 Segmented genomes — stacked vs separate.** The document specifies
stacked rows, largest first. The branch delivers separate files with a shared
depth y-limit and width proportional to length (so they compose honestly), in
GFF order. Either the document should record that decision or a
`StackedPlotter` remains a to-do; ordering largest-first is a one-line change
either way.

**H. §4 Decision logic — layout default.** "circular? → C modes" implies
circular genomes are drawn as circles by default. The CLI defaults to
`--layout linear` because a linear track carries far more labels legibly;
`--layout auto` gives the document's behaviour. This is a product decision
worth making explicit.

**I. §3 — nesting depth is bounded.** Three lanes fit between the depth ring
and the rim; deeper nesting is clamped to the innermost lane with a warning.
Real cases (nanovirus components, PCV) fit; a very dense circular genome would
need the depth ring made thinner or the arcs unrolled.

**J. Version.** Stages 1–3 change existing linear output (the synthetic sample
now has every box above the line, since none of its ORFs overlap, and its `p7`
takes the small-ORF purple). Per §5 that warrants a minor-version bump; not
applied on the branch — a release decision.

---

## 7. Checked against the published figures

The rules in §1–3 were compared with the figures themselves (viewed at
ictv.global/report_9th on 2026-09-29): Geminiviridae Figs 2 and 5,
Closteroviridae Figs 2, 3 and 5, Luteoviridae Fig 2. Figure images are
Elsevier/ICTV copyright and are not in this repository. The test data was
compared too, and a real fixture added.

### 7.1 The figures are hand-drawn; the rules approximate them

| Rule (§) | Holds in | Breaks in |
|---|---|---|
| L1 flip only on overlap (§2) | BYV: RdRp→p6 (50 nt gap) stay together below; every overlapping pair flips. CTV likewise. | BYV: CP is flipped below CPm despite a 70 nt gap (same colour, close — flipped for legibility; now a rule, see A′). LIYV RNA-1: p31 is above though RdRp below and no overlap. |
| First (5′-most) ORF above (§2) | BYV, LIYV, BYDV | **CTV: ORF1a is below, RdRp above.** Starting side is arbitrary; alternation is what matters. |
| Frameshift = step + `+1 FS` label (§1) | CTV (label present). | BYV shows the same step with no label. The step is simply the flip; the flip happens although the two boxes abut (7997/7999) rather than overlap — so **frameshift continuation flips even without overlap**. *Implemented: a frameshift ORF always flips; label at the junction on the partner's side.* |
| Readthrough = one box + bar (§1) | — | BYDV-PAV: ORF3→ORF5 readthrough is two abutting boxes on the **same** side, no bar. And ORF5 stays above although the rule would put it below with ORF4. So **readthrough continuation does not flip, frameshift does** — opposite behaviours the coordinates alone cannot distinguish. *Implemented: a readthrough extension takes the side of the ORF it abuts and gets a thin bar + `RT`; the bar is VirPlot's addition, the figure has none.* |
| ORF number outside, product inside (§1) | Closteroviridae: `ORF1a`/`ORF1b`/`66K 63K 100K` outside, `L-Pro Mtr Hel RdRp` inside. | **Luteoviridae: ORF numbers are inside the boxes** (`ORF1`…`ORF6`), sizes inside the product boxes beneath. Not universal. |
| Palette (§1) | Closteroviridae only. | **Geminiviridae uses a different palette**: CP green, MP yellow, Rep teal, TrAP pink, REn blue. Luteoviridae: CP pink, MP blue, RTD green, P0 (suppressor) green. Colour is consistent *within a family figure*, not across the Report. "Small 3′ ORFs purple" is BYV p21 only; BYV p6 is grey. |
| Circular: arcs "just outside the circle", nest inward (§3) | Mastrevirus: the arcs *are* the circle — thick coloured arcs replace the line where ORFs lie, the thin black line shows only in the LIR/SIR gaps. | Begomovirus: base arcs sit **on** the circle line; AC2/AC3/AC4 nest **inside** the circle. Nothing is drawn outside the ring except labels and the stem-loop. |
| Circular labels | — | All labels are **horizontal**, outside the ring for base arcs (`AV1 (CP)`), inside the circle for nested ones (`AC2 (TrAP)`); format `ORF (function)`. No rotated text. Centre carries the component name (`DNA-A`) only, no length. **No position ticks or scale.** |
| Circular IR (§3) | Begomovirus: the common region is a thick grey arc straddling 12 o'clock with the stem-loop icon on it, labelled `CRA`/`IR`. Mastrevirus: `LIR` at top, `SIR` at bottom. | — |
| Title `acronym (length nts)` (§1) | BYV, CTV, LIYV: `Beet yellows virus, BYV (15,468 nts)`. Luteoviridae: `Luteovirus, BYDV-PAV (5,677 nts)` (genus first). | — |
| Segmented: stacked rows, largest first, shared scale (§1) | LIYV: `RNA-1 (8,118 nts)` over `RNA-2 (7,193 nts)`, one title, same x-scale, both left-aligned at 5′. | — |
| 5′/3′ ends | `5′m⁷G?` / `3′OH` (Closteroviridae); grey oval `VPg` (PLRV, PEMV); bare `5′` (BYDV). | — |

### 7.2 What VirPlot does about it

* **Flip rule** — overlap, plus the same-colour-and-close clause (A′).
  Together they reproduce all nine BYV placements with the figure's palette;
  with a palette that gives CP and CPm different colours, CP stays beside CPm,
  which is then the legible choice. Since the figures do not agree with each
  other on the starting side, "first ORF above" is treated as a convention of
  ours, not theirs.
* **No outlines** — glyphs are flat colour, as in the figures and the original
  VirPlot output; there is no option to add an outline. Separation between
  neighbours comes from colour and from the flip rule above.
* **Frameshift / readthrough** — read from RefSeq's own attributes and from
  `Note=`; see F above and docs/GFF_GUIDE.md §8. With this, BYV matches
  Closteroviridae Fig 2 on eight ORFs from overlap and mechanism alone, and on
  all nine with the figure's palette (the CPm/CP legibility flip).
* **Palette** — the keyword classifier now colours only names that state a
  function (`RdRp`, `coat protein`, `movement`, `HSP70`, `silencing
  suppressor`); a bare `p6`/`p20` stays `default_color` until mapped. Colour
  should be thought of as a per-family style sheet: a Geminiviridae
  `function_palette` (CP green, MP yellow, Rep teal) is a five-line YAML,
  and the right place for it is a per-example spec, not the code default.
* **Circular ring** — VirPlot keeps its arcs outside the genome circle and
  nests toward it, because the inside of the circle holds the depth ring,
  which the ICTV figures do not have. Base-tier arcs could be drawn *on* the
  circle (Begomovirus style) at no cost; nesting inside the circle cannot
  coexist with the depth band. Rotated arc labels and rim ticks are VirPlot
  additions; a `--labels horizontal` option would match the figures when
  depth is not the point.
* **IR arc and `ORF (function)` labels** — need non-CDS features
  (`misc_feature`/`regulatory` rows) and a second name per feature in the
  GFF; see §6 D–F.

### 7.3 Test data

* `examples/sample.gff3` was written to *evoke* BYV but its nine ORFs are
  spaced 100 nt apart and never overlap, in a different order from the real
  genome, at half the length. It therefore cannot exercise the flip rule at
  all — under the rule every box sits above the line — and it should not be
  read as an ICTV comparison. It remains a fine depth-plotting fixture.
* `examples/byv.gff3` (new) carries the real NC_001598.1 coordinates with
  the ICTV ORF names and short products; `examples/byv.sam` is a flat
  simulated read set from `examples/make_reads.py`; `examples/byv.yml` is
  the Closteroviridae Fig 2 palette. Rendered with that palette it
  reproduces all nine placements of the figure
  (`tests/test_conventions.py::test_byv_with_figure_palette_matches_ictv_on_every_orf`);
  the pure overlap rule alone gets seven
  (`test_byv_fixture_layout_versus_ictv_figure`).
* `examples/grbv.gff3` is a Grablovirus, not in the 9th Report; its nearest
  figure is Begomovirus DNA-A. The C2 conventions (V clockwise right, C
  anticlockwise left, nesting) hold; RefSeq names `V1 protein` etc. carry no
  function, which is why `grbv.yml` maps colours by hand.
