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
| 4 | **Done.** 5′/3′ end labels and the VPg oval; ORF name outside the glyph from `gene=` / `Name=` / `locus_tag=`, omitted when it repeats the product (see D); default title `name (length nts)`; **frameshift step + `±1 FS` label and readthrough bar + `RT`**, read from RefSeq's `exception=ribosomal slippage` / `transl_except=` or from `Note=`; **polyprotein domain dividers** from RefSeq's `mature_protein_region_of_CDS` / `mat_peptide` rows, in both layouts (see F); **sgRNA rows** beneath the line, longest first, from a marked transcript row (see K). The FS/RT marks and the sgRNA ladder are linear-layout only, and the circular layout now warns rather than dropping them silently. | `plotting.py`, `parsers.py`, `layout.py`, `models.py` |
| 5 | **Done by composing.** One figure per segment, largest first, on a shared scale and y-limit; `bin/stack_figures.py` stacks them into the figure §1 describes. A `StackedPlotter` was considered and rejected on a measurement — see G. | `cli.py`, `plotting.py`, `bin/stack_figures.py` |
| 6 | Done. `CircularPlotter`: origin at 12 o'clock with stem-loop icon (moved to a GFF `stem_loop` row when there is one), arcs just outside the circle with arrowheads, clockwise for + and anticlockwise for −, overlap nests inward; **IR / UTR as a grey arc astride the circle** from the GFF's non-coding rows; depth as an inner ring. `--layout circular|auto`. | `plotting.py` |

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
`product=V1 protein`. The shipped `examples/grbv.gff3` is therefore curated
(`gene=V1;product=CP`, RefSeq name kept in `Note=`), which is what users are
expected to do; docs/GFF_GUIDE.md §5 says so.

**Settled, against a real annotation.** The outside label is now taken from
`gene=`, then `Name=`, then `locus_tag=`, so an uncurated RefSeq record labels
its ORFs `N761_gp1` rather than not at all, and a hand annotation that puts the
ORF name in `Name=` is read as intended.

Testing that against the GLRaV-13 H8881 annotation used in the FPS paper
(17,564 nt, 13 ORFs) turned up a second rule the document does not state and
real files need badly: **an outside label equal to the product is dropped.**
That file carries `gene=` on every ORF, but for eleven of thirteen the gene and
the product are the *same string* (`gene=CP;product=CP`, `gene=p53;product=p53`),
because the ORF has one name and no separate known function. Writing it inside
the box and again above it is noise, not the convention's two names. The same
holds in the shipped `examples/grbv.gff3`, where `V3` and `C3` are named twice
for exactly that reason. With the rule, GLRaV-13's one genuine two-name
feature — `Name=polyprotein_1a;product=Methyltransferase/helicase` — is the only
ORF that gets a label above its box, which is what §1 is asking for.

§1 should therefore say: outside label from `gene=`, else `Name=`, else
`locus_tag=`; inside label from `product=`; and the outside label is omitted
when the two would read the same.

**E. §1 Palette — RefSeq product names are often function-free.** The keyword
classifier covers RdRp/Rep/CP/MP/HSP70/p2x-style names; anything else falls to
grey unless mapped. That is the right failure mode, but the document's "map
product names by default" will silently leave many RefSeq genomes grey.

**F. §1 Expression features need an encoding — three now have one.** Frameshift
and readthrough are read from what RefSeq already writes (`exception=ribosomal
slippage` on the rows of a joined CDS; `transl_except=` on a readthrough
product) and from VirPlot's `Note=+1 frameshift` / `Note=readthrough`
convention (docs/GFF_GUIDE.md §8). The frameshift sign is derived from the
join geometry (skip one base = +1, re-read one = −1). Polyprotein domains
are read from RefSeq's `mature_protein_region_of_CDS` rows (GenBank
`mat_peptide`; `Parent=` names the CDS, or containment is used) and drawn as
the §1 rule says — one box, thin lines between segments, a name in each —
with one addition the rule does not state and Potyviridae Fig. 2 shows:
**each segment is coloured on its own** (by `color_mapping`, then function
words), the box's own `product` is not written, and a name too long for its
segment moves outside the box (`6K1`, `6K2` in that figure). The document
should add the per-segment colour and the outside-label fallback to the
Polyprotein row of §1. Non-coding features (§1 palette, §3 baseline) are
now read too: RefSeq's `five_prime_UTR` / `three_prime_UTR` / `stem_loop`
rows always, its catch-all `sequence_feature` (`misc_feature`) /
`regulatory_region` / `repeat_region` rows only when their words say
intergenic / common region / IR / UTR / stem-loop, so motif annotations
inside ORFs stay off the line. A region is a grey bar or arc astride the
genome line; a stem-loop is the hairpin icon at its own position (which on a
circle replaces the default icon at position 1); a region holding an unnamed
stem-loop is named once at the hairpin, as Geminiviridae Fig. 5 writes
`CRA`. The document should say that these are drawn **on** the line and take
no part in the ORF layout. sgRNAs are now encoded and drawn; see K.

**G. §1 Segmented genomes — rendered separately, composed afterwards.** The
document specifies stacked rows, largest first. VirPlot renders one figure per
segment and stacks them with `bin/stack_figures.py` rather than drawing a
stacked figure directly. **This is a decision, not a missing feature**, and the
measurement behind it is this: within one run the figure width is proportional
to genome length and matplotlib's margins are fractional, so the panels already
share an exact x-scale. Across TSWV's L/M/S (8,897 / 4,821 / 2,916 nt — a 3×
spread) the measured scale drift is **0.0000%** and the plot areas start at the
same x to within 0.00 pt. Stacking is therefore a document operation, and a
`StackedPlotter` would re-implement in matplotlib a figure that already
composes exactly.

What the decision buys and costs. The panels stay independently useful (a
single segment is a figure in its own right, which a stacked-only renderer
would lose), the circular case needs no second layout pass, and the composer is
~170 lines of stdlib against a 150–250-line renderer that would have to handle
unequal-width axes in one `GridSpec` and a radial budget for circular
components. The cost is one extra command, which `docs/make_figures.sh` shows.

Three pieces close the gap to the document's wording: segments now render
**largest first** by default (`--rna-order length`), `--bare-x` drops the
repeated x-axis from every panel but the bottom one, and the composer orders
panels largest first and **aligns their plot areas** — a panel whose depth axis
carries wider tick labels starts further right once the figure is cropped to a
tight bounding box (0.63 pt on LIYV), which hand pasting cannot correct.
`tests/test_stacking.py` guards the shared-scale property itself, since the
decision is only sound while it holds.

§1 needs no change: it describes the finished figure, and the composed figure
is that figure — one row per segment, largest first, on a shared scale. It is
§5's stage wording ("→ stacked rows") that should say the rows are composed
from per-segment figures rather than drawn by one renderer.

**H. §4 Decision logic — layout default. Settled: the document wins.** A
circular genome is now drawn as a circle without being asked, as §4 says.
`--layout` defaults to `auto`, and spec.yml gains a `layout` key so a
per-example default lives beside the data it styles; precedence is `--layout`,
then the YAML, then `auto`.

**`--layout linear` stays, and becomes load-bearing.** Three reasons it cannot
be replaced by editing the GFF. First, `Is_circular=true` is a fact about the
*molecule*, not the figure: it drives depth wrapping at the origin and the
splitting of origin-crossing features, so deleting it to get a linear picture
would silently change the data handling as well. Presentation needs its own
control, which is why `--topology` (what the molecule is) and `--layout` (how it
is drawn) are separate axes and should stay separate. Second, the circular
layout cannot carry everything the linear one can — the `±1 FS` and `RT` marks
are linear-only (§6 F), as are sgRNA rows (§6 K), so a circular genome with a
frameshift or a transcript ladder needs a way back to a track; `CircularPlotter`
now warns rather than dropping those marks silently. Third, legibility: a linear
track carries far more labels, and circular nesting is bounded at three lanes
before it clamps (§6 I).

`docs/make_figures.sh` is the immediate proof. Its `grbv_linear` figure asked
for no layout and relied on the old default; under the new one it would have
become a circle, so it now passes `--layout linear` explicitly.

The original note, kept for the record: the CLI defaulted to
`--layout linear` because a linear track carries far more labels legibly;
`--layout auto` gives the document's behaviour. This is a product decision
worth making explicit.

**I. §3 — nesting depth is bounded.** Three lanes fit between the depth ring
and the rim; deeper nesting is clamped to the innermost lane with a warning.
Real cases (nanovirus components, PCV) fit; a very dense circular genome would
need the depth ring made thinner or the arcs unrolled.

**K. §1 sgRNAs — encoded as a curation convention, and linear-only by
design.** The document asks for "shorter lines stacked beneath the genome, 5′
aligned to their start, labelled". The branch draws exactly that, longest row
first, from a transcript row (`mRNA`, `transcript`, `ncRNA`, `misc_RNA`,
`primary_transcript`, `sequence_feature`, `misc_feature`) whose `Note=` says
`sgRNA` / `subgenomic`, with `gene=` naming what it expresses
(docs/GFF_GUIDE.md §12). The marker is required rather than inferred from the
type, because RefSeq writes real `mRNA` rows — spliced mastrevirus transcripts
— that are not sgRNAs. `start` is the 5′ terminus; an `end` at or past the
genome length means the 3′ end, which covers the usual 3′-coterminal set,
while an earlier end is kept for the 5′-proximal sgRNAs closteroviruses also
make. `examples/byv.gff3` carries the ladder, with the two 5′ termini mapped
by Vitushkina et al. (2002, *Virology* 297:299–307) marked as such and the
other five flagged as illustrative placements.

Two things the document should record. First, **these rows are annotation and
never a depth track**: a plant virus sgRNA is co-linear with the genome and
carries no leader junction (unlike *Nidovirales*), so a read from an sgRNA
cannot be told from a genomic read at the same coordinate and per-sgRNA
coverage is not recoverable from short reads. What the set leaves is a step in
the *aggregate* depth at each 5′ end, height proportional to abundance — which
is the reason to draw the ladder on the same x-axis as the depth trace, and
which explains a 5′-to-3′ coverage ramp that would otherwise read as a failed
assembly.

Second, **sgRNA rows are drawn in the linear layout only, deliberately.** The
ICTV draws no transcript rows on circular genomes: geminivirus transcription
is bidirectional from the IR with overlapping transcripts rather than a
3′-coterminal set (Geminiviridae chapter, Fig. 2 shows ORFs and the stem-loop
and nothing else), and nanovirus components carry one ORF each. Nesting seven
near-complete arcs would also be unreadable. `--layout circular` warns and
skips. The genuine circular case is *Caulimoviridae* — CaMV's 35S and 19S —
which is two arcs, not a ladder, and would fit the existing nesting lanes if
it is ever wanted. §1 should say the sgRNA row is a linear-genome convention.

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
| Circular: arcs "just outside the circle", nest inward (§3) | Mastrevirus: the arcs *are* the circle — thick coloured arcs replace the line where ORFs lie, the thin black line shows only in the LIR/SIR gaps. | Begomovirus: base arcs sit **on** the circle line; AC2/AC3/AC4 nest **inside** the circle. Nothing is drawn outside the ring except labels and the stem-loop. *Option: `circular_arcs: on_circle` (innermost tier astride the line; deeper nesting still steps outward).* |
| Circular labels | — | All labels are **horizontal**, outside the ring for base arcs (`AV1 (CP)`), inside the circle for nested ones (`AC2 (TrAP)`); format `ORF (function)`. No rotated text. Centre carries the component name (`DNA-A`) only, no length. **No position ticks or scale.** *Option: `circular_labels: horizontal`.* |
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
  VirPlot output (`--no-border`, the default); `--border` adds a black
  outline for those who want it. Separation between
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
* **Circular ring** — by default VirPlot keeps its arcs outside the genome
  circle and nests toward it, because the inside of the circle holds the
  depth ring, which the ICTV figures do not have. Two spec.yml keys bring it
  to the figures' style: `circular_arcs: on_circle` sets the innermost arcs
  astride the circle line (Begomovirus/Mastrevirus style; outer tiers still
  step outward, since nesting inside cannot coexist with the depth band), and
  `circular_labels: horizontal` writes level `ORF (product)` labels outside
  the ring. The figures label nested arcs inside the circle; VirPlot instead
  steps a nested arc's label one line away from the equator per tier so it
  does not overprint the arc it nests in. Rim ticks remain a VirPlot addition.
* **IR arc** — drawn from the GFF's non-coding rows (docs/GFF_GUIDE.md §10)
  as a grey arc astride the circle, under the ORF arcs, with the hairpin icon
  at the `stem_loop` row's position and the name written once over it, as
  Fig. 5 draws `CRA`; the Mastrevirus `LIR`/`SIR` are the same two rows.
  `examples/grbv.gff3` and `examples/tgmv.gff3` carry their intergenic /
  common regions and hairpins this way (RefSeq annotates neither for those
  records, so the rows were located from the sequences and say so);
  `examples/fbnyv.gff3` keeps RefSeq's own `stem_loop 1..33` rows. The depth
  ceiling label, which stood against the circle at 12 o'clock, hangs inward
  when such an arc straddles the origin. `ORF (function)` labels come from
  `gene=` + `product=` with `circular_labels: horizontal`; see §6 D.

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
  anticlockwise left, nesting) hold. RefSeq's `V1 protein` etc. carry no
  function, so the shipped GFF is curated to `gene=V1;product=CP` and the
  default palette colours it; the RefSeq names are kept in `Note=`.
* `examples/pvy.gff3` (potato virus Y, NC_001616.1) is the polyprotein
  fixture, compared with Potyviridae Fig. 2 (drawn from TEV, same
  organisation). With `pvy.yml` (the figure's per-domain colours, sampled
  from the image) the ten segments, their dividers and names match; `6K1`
  and `6K2` go outside the box as in the figure. Two differences remain:
  the figure draws the polyprotein **on** the line with PIPO above it, where
  L1 puts the polyprotein above and flips PIPO below; and the figure writes
  a second row of names under the box (`35K`, `52K`, … sizes inside,
  functions below), which VirPlot has no second label slot for — the
  curated file keeps the function names. PIPO is really reached by
  polymerase slippage (P3N-PIPO); the file carries the figure's `+2 fs`
  wording as `Note=+2 frameshift`, so the label reads `+2 FS`.
