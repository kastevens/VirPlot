# Writing a GFF3 file for VirPlot

VirPlot draws the annotation track from a GFF3 file. It reads a deliberately
small subset of GFF3 — two feature types and four attributes — and applies a
few conventions of its own for things GFF3 does not standardise (circularity,
frameshifts, function colours). This guide says exactly what it reads, what it
ignores, and how to write each situation so the figure comes out as intended.
Examples of every case are in [`examples/`](../examples/): `byv.gff3` (linear,
one strand, with a `Note=+1 frameshift`), `grbv.gff3` (circular, both strands),
`sample_multi.gff3` (two segments).

## 1. The minimum

```gff3
##gff-version 3
NC_001598.1	RefSeq	region	1	15480	.	+	.	ID=NC_001598.1;Is_circular=false
NC_001598.1	RefSeq	CDS	108	7997	.	+	0	ID=orf1a;gene=ORF1a;product=L-Pro/Mtr/Hel
NC_001598.1	RefSeq	CDS	9608	11404	.	+	0	ID=orf3;gene=ORF3;product=Hsp70h
```

Nine tab-separated columns per row. VirPlot uses:

| Column | Used for |
|---|---|
| 1 `seqid` | Which molecule the row belongs to; must match the SAM/BAM `@SQ` name, the FASTA header and the depth-file sequence column (see §7). |
| 3 `type` | Only `region` and `CDS` are read. Everything else is skipped silently. |
| 4–5 `start`, `end` | 1-based, inclusive, as GFF3 specifies. |
| 7 `strand` | `+` or `-`. Anything else (`.`, `?`) is treated as `+`. |
| 9 `attributes` | `product`, `gene`, `Is_circular`, `Note` (see below). |

Columns 2 (source), 6 (score) and 8 (phase) are read past and ignored.
Comment lines (`#…`) and blank lines are ignored, including the
`##sequence-region` pragma — the `region` **row** is what defines the molecule.

## 2. The `region` row: one per molecule

Every molecule needs exactly one `region` row. It gives VirPlot the length
(column 5) and, optionally, the topology:

```gff3
NC_022002.1	RefSeq	region	1	3206	.	+	.	ID=NC_022002.1;Is_circular=true
```

* `end` is the genome length. Depth arrays are sized from it, so it must
  match the reference the reads were mapped to; VirPlot warns when a SAM/BAM
  header or the `-x` FASTA disagrees.
* `Is_circular=true` marks the molecule circular (§4). This is the attribute
  GFF3 itself defines for the purpose, and NCBI RefSeq writes it, so a
  downloaded RefSeq GFF3 already carries the right value. Absent or `false`
  means linear. `--topology circular|linear` overrides the file.
* Several `region` rows with different `seqid`s make a segmented genome (§6).
* A CDS whose `seqid` has no `region` row is dropped with a warning; a file
  with no `region` row at all is an error.

The molecule's display name (in titles and the circular hub) is the `seqid`.

## 3. `CDS` rows: one per ORF

Each `CDS` row becomes one glyph. Two attributes name it:

| Attribute | Drawn | Purpose |
|---|---|---|
| `product=` | **Inside** the glyph (or just outside it when the ORF is under 500 bp) | The protein or domain: `RdRp`, `CP`, `Hsp70h`. Also drives the default colour (§5). Missing → `unknown`. |
| `gene=` | **Outside** the glyph, small grey text | The ORF name: `ORF1a`, `AC1`, `V2`. Optional. |

This is the ICTV figure convention: ORF number outside, product inside. Keep
`product` short — it has to fit in the box. The long RefSeq description can go
in `Note=`, which VirPlot ignores (and so never draws).

Other attributes (`ID`, `Parent`, `locus_tag`, `Dbxref`, `protein_id`, …) are
accepted and ignored. `gene` rows, `mRNA` rows, `mat_peptide` rows are ignored
too, so a RefSeq GFF3 can be used as is — with two caveats in §8.

## 4. Strand

Column 7 decides two things.

**Which side of the line.** When every CDS is on one strand, VirPlot uses the
one-strand layout: boxes sit on the line, the 5′-most above, and a box flips
to the other side only when it overlaps its upstream neighbour (or is the same
colour and within 1 % of the genome of it — §5). When CDS rows appear on
*both* strands, VirPlot switches to the two-strand layout: `+` ORFs above the
line drawn as rightward arrows, `-` ORFs below as leftward arrows, and
same-side overlaps tier outward instead of flipping. `overlap_mode: flip|tier`
in the spec.yml forces either.

**Arrow direction on a circle.** In the circular layout, `+` arcs run
clockwise from the origin and `-` arcs anticlockwise, with the arrowhead at
the reading end.

So for a geminivirus, write the virion-sense ORFs with `+` and the
complementary-sense ones with `-`, exactly as RefSeq does:

```gff3
NC_022002.1	RefSeq	CDS	292	807	.	+	0	ID=v2;product=V2 protein
NC_022002.1	RefSeq	CDS	2250	3044	.	-	0	ID=c1;product=C1 protein
```

Coordinates of a `-` strand feature are still written `start < end` in
genome coordinates; VirPlot works out the reading direction from the strand.

## 5. Colour: name the function

VirPlot colours by **predicted function**, never by strand or reading frame
(the ICTV convention). The colour of a glyph is decided in this order:

1. An explicit entry for the exact `product` string in the spec.yml
   `color_mapping`.
2. Otherwise, a **function class** guessed from words in `product`, drawn in
   that class's palette colour.
3. Otherwise `default_color` (grey).

The classes, their default colours, and the words that select them:

| Class | Colour | `product` contains (case-insensitive) |
|---|---|---|
| replicase | `#f5b041` yellow-orange | `RdRp`, `RNA-dependent`, `polymerase`, `replicase`, `replication-associated`, `Rep`, `RepA`, `methyltransferase`, `Mtr`, `helicase`, `Hel`, `protease`, `papain`, `nsp1`…`nsp9`, `P1a`/`P1ab`, `ORF1`/`ORF1a`/`ORF1b`, `polyprotein` |
| hsp70 | `#7dbb4f` green | `HSP70`, `HSP90`, `heat shock` |
| capsid | `#e75480` magenta | `coat`, `capsid`, `nucleocapsid`, `CP`, `CPm`, `CPh`, `N protein`, `virion protein`, `structural` |
| movement | `#3ea6b5` teal | `movement`, `MP`, `cell-to-cell`, `triple gene block`, `TGB1`… |
| suppressor | `#9b6fc4` purple | `silencing`, `suppressor`, `VSR`, `HC-Pro` |
| noncoding | `#6f6f6f` dark grey | `intergenic`, `IR`, `UTR`, `stem-loop`, `hairpin`, `non-coding` |
| *(none)* | `default_color` `#9f9f9f` | anything else |

Rules are tried in the order listed; the first match wins, so
`coat protein` is capsid even though it contains `protein`.

Two consequences worth knowing:

* **Uninformative names stay grey on purpose.** `p20`, `V1 protein`,
  `hypothetical protein`, `ORF6` carry no function, so VirPlot does not guess
  (BYV's `p6` is a membrane protein; its `p21` is a silencing suppressor —
  the name alone cannot tell). Either name the function in `product`
  (`p21 silencing suppressor`) or pin the colour in `color_mapping`.
  `examples/grbv.yml` does the latter for GRBV, whose RefSeq products are all
  `Vn protein` / `Cn protein`.
* **"Putative" lightens the colour.** If `product` also contains `putative`,
  `probable`, `possible`, `hypothetical`, `predicted`, `proposed` or `-like`,
  the class colour is tinted 35 % towards white — the ICTV figures' convention
  for unconfirmed function. `putative movement protein` is pale teal.

Any class colour can be replaced in the spec.yml under `function_palette:`;
this is how to adopt a family's own scheme (Geminiviridae figures, for
instance, draw CP green and MP yellow).

Colour also affects layout: in the one-strand layout, two **adjacent boxes of
the same colour** closer than 1 % of the genome flip to opposite sides so they
do not read as one box (BYV's CPm and CP). Giving neighbours the same function
word therefore separates them; giving them different words keeps them
together.

## 6. Circular genomes

Set `Is_circular=true` on the `region` row. Three things follow:

* **Depth wraps.** Reads that run off the end continue from position 1, and
  positions past the end are taken modulo the length. This is what makes a
  BAM produced by padding the reference (the usual fix for aligners, which
  know nothing of circles) come out continuous. Without the flag the same
  data shows a false dip at the origin.
* **Features may cross the origin.** Write such a feature either with
  `end` past the genome length or with `start > end`; both occur in real
  annotations and both are understood:

  ```gff3
  # a 400-nt ORF starting at 3000 on a 3206-nt circle, two equivalent spellings
  seq	.	CDS	3000	3399	.	+	0	ID=x;product=Rep
  seq	.	CDS	3000	193	.	+	0	ID=x;product=Rep
  ```

  It is drawn as the two arcs (or two boxes) it occupies, with one label and,
  in arrow modes, one arrowhead at the reading end. On a *linear* molecule
  there is no origin to cross, so these spellings produce a wrong glyph
  (`end` past the length runs off the axis; `start > end` collapses); use
  them only with `Is_circular=true`.
* **Layout.** The linear layout drops the 5′/3′ marks and shows the backbone
  continuing past both edges. `--layout circular` draws the molecule as a
  circle with position 1 at the top; `--layout auto` does so for every
  circular molecule in the file.

Two-strand circular genomes (geminiviruses) need nothing extra: strand from
column 7 gives arrow direction, overlaps nest inward, largest arc outermost.

## 7. Matching names across files

The `seqid` in column 1 is the key that joins the GFF to everything else:

| File | Must match `seqid` |
|---|---|
| SAM/BAM | the `@SQ SN:` name (`RNAME` in records) |
| `samtools depth` output | column 1 |
| `-x` reference FASTA | the header's first word |

If a depth source names exactly **one** sequence, VirPlot uses it whatever it
is called (with a warning); with several and no match, it stops and lists the
candidates, and `--ref NAME` picks one. Using the RefSeq accession *with its
version* (`NC_022002.1`) everywhere is the least error-prone choice, since that
is what the RefSeq FASTA header and therefore the aligner's `@SQ` will carry.

## 8. Frameshifts and readthrough

GFF3 has no field for "this ORF is reached from that one by a frameshift" or
"by reading through a stop codon". VirPlot recognises the two ways such ORFs
are actually written and draws them as the ICTV figures do.

### Frameshift

**RefSeq form — nothing to change.** RefSeq encodes a frameshifted
polyprotein as one CDS in several rows sharing an `ID`, each carrying
`exception=ribosomal slippage`:

```gff3
NC_001598.1	RefSeq	CDS	108	7997	.	+	0	ID=cds-NP_041872.1;exception=ribosomal slippage;product=fusion protein of ...
NC_001598.1	RefSeq	CDS	7999	9393	.	+	0	ID=cds-NP_041872.1;exception=ribosomal slippage;product=fusion protein of ...
```

VirPlot groups the rows by `ID`, orders them in translation direction, draws
each segment as its own box, flips every segment after the first across the
line, and writes the shift at the junction. The sign comes from the
coordinates: a junction that skips one base (`…7997`, `7999…`) is `+1 FS`; one
that re-reads a base (`…13468`, `13468…`, as in coronaviruses) is `−1 FS`.
Because RefSeq gives every row the same long product, only the first segment
is labelled with it; the continuation shows the `FS` mark instead. Rows that
share an `ID` *without* a slippage note (a spliced CDS) are left as separate,
unmarked boxes.

**Hand-written form — name the two ORFs.** To get the figure's `L-Pro/Mtr/Hel`
and `RdRp` labels, write the segments as two CDS rows and put the mechanism in
`Note=` on the downstream one:

```gff3
NC_001598.1	RefSeq	CDS	108	7997	.	+	0	ID=orf1a;gene=ORF1a;product=L-Pro/Mtr/Hel
NC_001598.1	RefSeq	CDS	7999	9393	.	+	0	ID=orf1b;gene=ORF1b;product=RdRp;Note=+1 frameshift from ORF1a
```

Any `Note` containing `frameshift` marks the row; `+1` or `-1` in the note
gives the sign, otherwise it is inferred from the gap to the nearest upstream
ORF on the same strand. `examples/byv.gff3` is written this way.

A frameshift continuation **always** flips across the line, whether or not it
overlaps — that step is how every ICTV figure shows the event — so the
same-colour rule of §5 no longer has to do this job.

### Readthrough

**RefSeq form.** A readthrough product is annotated as a CDS that starts
where the shorter ORF starts and runs past its stop, with a `transl_except=`
attribute recording the read-through codon (TMV's 183K beside its 126K, for
instance):

```gff3
NC_001367.1	RefSeq	CDS	69	3419	.	+	0	ID=p126;product=126 kDa replicase
NC_001367.1	RefSeq	CDS	69	4919	.	+	0	ID=p183;product=183 kDa replicase;transl_except=(pos:3417..3419%2Caa:OTHER)
```

VirPlot trims the readthrough CDS to the part beyond the shorter ORF
(3420–4919 here) and draws it as a second box **on the same side**, abutting
the first, with a thin bar at the read-through stop labelled `RT` — the
BYDV ORF3/ORF5 picture.

**Hand-written form.** Give the extension its own coordinates and
`Note=readthrough`:

```gff3
A	.	CDS	100	1000	.	+	0	ID=orf3;gene=ORF3;product=CP
A	.	CDS	1001	1800	.	+	0	ID=orf5;gene=ORF5;product=readthrough domain;Note=readthrough of ORF3 stop
```

The extension takes the side of the ORF whose end it abuts (not merely the
previous ORF by coordinate — a nested ORF4 can sit between ORF3 and ORF5, as
in luteovirids), and never flips.

Both marks are drawn in the linear layout only for now; the circular layout
places the ORFs correctly but omits the `FS`/`RT` text.

### Not yet drawn

Polyprotein domain dividers (`mat_peptide` rows), subgenomic RNAs and
non-coding features (`misc_feature`, `regulatory`, `stem_loop`) are ignored by
the parser; see §6 F and §7.2 of
[`ictv_drawing_conventions.md`](ictv_drawing_conventions.md) for the intended
encodings. The `noncoding` colour class exists for the day they are read.

## 9. Checklist

- [ ] `##gff-version 3` first line; nine tab-separated columns.
- [ ] One `region` row per molecule, `end` = genome length,
      `Is_circular=true` if it is a circle.
- [ ] `seqid` identical across GFF, SAM/BAM `@SQ`, FASTA header, depth file.
- [ ] One `CDS` row per ORF; `product=` short and naming the **function**;
      `gene=` for the ORF name; long descriptions in `Note=`.
- [ ] Strand `+`/`-` correct — it chooses the layout and the arrow direction.
- [ ] Frameshifts: either leave RefSeq's `exception=ribosomal slippage` rows as
      they are, or split into named ORFs with `Note=+1 frameshift`; readthrough
      via `transl_except=` or `Note=readthrough`.
- [ ] Origin-crossing features only on `Is_circular=true` molecules.
- [ ] Products with no function word either renamed or pinned in
      `color_mapping`.
