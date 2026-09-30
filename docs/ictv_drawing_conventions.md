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
