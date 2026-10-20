# To do

What is left on the ICTV work. Tick an item off or delete it in the same
commit that finishes it. The first section matches the ICTV-examples phase of
the JOSS submission checklist.

## ICTV examples

- [ ] **Examples gallery.** Write `docs/EXAMPLES.md` with one row per feature
      and one column per example, so every feature is shown at least once.
      Give each example its figure, the command that made it, and its sources
      (RefSeq accession, ICTV figure). Link the ICTV figures rather than
      copying them; they are Elsevier/ICTV copyright.
- [ ] **Multi-sample depth example.** Two or three samples stacked, coverage
      gaps shaded (`--shade-breaks`), CSV interval and gap reports
      (`--report`), a symlog y-axis, and `samtools depth` file input. None of
      these five features has an example yet. Use real FPS reads if they can
      be shared, simulated reads otherwise.
- [ ] **One script for every figure, with tests.** `docs/make_figures.sh`
      covers BYV, GRBV, the synthetic sample, PVY, TSWV and GLRaV-13. Add
      BYDV, PLRV, LIYV, TGMV, FBNYV and PVX. SVG output is byte-reproducible
      under a fixed matplotlib (3.10.9 made the current set), so a test can
      compare file hashes and catch any drawing change.
- [ ] **Keep example data small.** Everything in `examples/` is under 260 KB
      except `glrav13.bam` (4.9 MB). Subsample it, or move it to Zenodo with
      a fetch script.

## Spec changes not yet written into the rules

`docs/ictv_drawing_conventions.md` §6 records where implementing the spec
showed the spec itself needs changing. The code already does what each item
says; the rules in §1–§5 still carry the old wording. Only A′ (the
same-colour flip, §2) has been written in. E, H and J need no change.

- [ ] **A — §2 L1.** State the overprint fallback: if the side a flip wants is
      taken, use the other side; if both are taken, tier outward.
- [ ] **B — §3.** Overlapping arcs nest largest first, so the largest arc is
      outermost whatever the coordinates.
- [ ] **C — §3 C2.** Direction is set by strand and nesting by overlap; the
      two halves of the circle follow from the coordinates, not a rule.
- [ ] **D — §1 Labels.** Outside label from `gene=`, else `Name=`, else
      `locus_tag=`; inside label from `product=`; the outside label is dropped
      when it repeats the product.
- [ ] **F — §1 Polyprotein.** Each domain segment is coloured on its own, the
      box's own product is not written, and a name too long for its segment
      goes outside. Non-coding features are drawn on the line and take no
      part in the ORF layout.
- [ ] **G — §5 stage 5.** Segment rows are composed from per-segment figures
      (`bin/stack_figures.py`), not drawn by one renderer.
- [ ] **I — §3.** Nesting is bounded at three lanes; deeper nesting clamps to
      the innermost lane with a warning.
- [ ] **K — §1 sgRNAs.** The sgRNA row is a linear-genome convention, and the
      ladder predicts a depth staircase that a given library may or may not
      show.

Known differences from the published figures are recorded in §7 of the same
document and are not on this list.

## Stranded depth and MAF (`feature/stranded-maf`)

Stranded depth is in (`--stranded`; `docs/stranded_depth.md`). Left:

- [ ] **Per-strand reports.** Coverage intervals, gaps and the CSV reports
      use the two strands' sum. Add + and − columns, or per-strand reports.
- [ ] **Re-measure the sgRNA staircase on the + strand.** The − strand reads
      in `glrav13.bam` are concentrated toward the 5′ end (+ : − is 4.4 in
      the first 2 kb, 13–25 elsewhere). Adding them to the 5′ end flattens
      the aggregate 5′→3′ comparison in `ictv_drawing_conventions.md` §6 K
      and `GFF_GUIDE.md` §11.
- [ ] **A stranded circular dataset.** Nothing in the repository checks the
      circular layout against a real stranded library, or measures how
      balanced the strands of a DNA virus are.
- [ ] **Minor-allele frequency.** Decode SEQ and QUAL, count alleles per
      position, and draw MAF against depth; the strand split gives a
      strand-bias filter (`docs/stranded_depth.md` §6).
