# Writing a GFF3 file for VirPlot

VirPlot draws the annotation track from a GFF3 file. It reads a deliberately
small subset of GFF3 — a few feature types and a handful of attributes — and
applies a few conventions of its own for things GFF3 does not standardise
(circularity, frameshifts, function colours). This guide says exactly what it
reads, what it ignores, and how to write each situation so the figure comes
out as intended. Examples of every case are in [`examples/`](../examples/):
`byv.gff3` (linear, one strand, with a `Note=+1 frameshift`), `grbv.gff3`
(circular, both strands, RefSeq labels curated), `pvy.gff3` (a polyprotein
with its mature-protein rows), `sample_multi.gff3` (two segments).

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
| `gene=` | **Outside** the glyph, small grey text | The ORF name: `ORF1a`, `AC1`, `V2`. Optional. `Name=` is used when there is no `gene=`, then `locus_tag=`. |

This is the ICTV figure convention: ORF number outside, product inside. Keep
`product` short — it has to fit in the box. The long RefSeq description can go
in `Note=`, which VirPlot ignores (and so never draws).

**One name, written once.** If the outside label would read the same as
`product`, it is dropped and only the box is labelled. Real annotations very
often carry `gene=CP;product=CP` or `gene=p53;product=p53`, because the ORF has
a name and no separately known function — `examples/grbv.gff3` does it for `V3`
and `C3`. Writing the name inside the box *and* above it is noise. So give the
two attributes two different things to say, or just one of them:

```gff3
# two names: ORF name outside, function inside — what the convention wants
A	.	CDS	1071	7655	.	+	0	ID=a;Name=polyprotein_1a;product=Methyltransferase/helicase
# one name: labelled once, inside the box
A	.	CDS	9248	9439	.	+	0	ID=b;gene=p7;product=p7
```

Other attributes (`ID`, `Parent`, `locus_tag`, `Dbxref`, `protein_id`, …) are
accepted and ignored. `gene` rows and `mRNA` rows are ignored too, so a RefSeq
GFF3 can be used as is — with two caveats in §8. Mature-protein rows
(`mature_protein_region_of_CDS` / `mat_peptide`) are read: they split a
polyprotein's box into domains (§9). Non-coding rows — UTRs, intergenic
regions, stem-loops — are drawn on the genome line (§10).

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
NC_022002.1	RefSeq	CDS	292	807	.	+	0	ID=v2;gene=V2;product=putative MP
NC_022002.1	RefSeq	CDS	2250	3044	.	-	0	ID=c1;gene=C1;product=RepA
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
  `examples/grbv.gff3` shows the former: RefSeq calls its products
  `V1 protein` … `C3 protein`, and the shipped file rewrites them to
  `gene=V1;product=CP`, `gene=C1;product=RepA` and so on, keeping the RefSeq
  name in `Note=`. **Expect to do this curation for your own genomes** — a
  header comment in the GFF recording what was changed and on what authority
  is the right place to say so.
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

Set `Is_circular=true` on the `region` row. Four things follow:

* **The figure is drawn as a circle**, unless something says otherwise.
  `--layout` (`auto` by default) decides, falling back to a `layout:` key in
  spec.yml. Use `--layout linear` — or `layout: linear` in the YAML — for the
  linear track of a circular genome; do **not** delete `Is_circular=true` to
  get one, because that would also turn off the two behaviours below. The
  molecule's shape and the figure's shape are separate settings on purpose.
  Note the `±1 FS` / `RT` marks (§8) and sgRNA rows (§12) are drawn in the
  linear layout only, and a circular figure warns when it has to leave them out.

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

Any `Note` containing the word `frameshift` marks the row. **Write the sign
(`+1` or `-1`) in the note.** Without it VirPlot infers the sign from the gap
between this row and the nearest upstream ORF on the same strand, which is
only right when the two coordinates meet at the slippage site as they do in a
RefSeq join; if ORF1b is written from its first full codon, or with any other
convenient boundary, the inferred sign will be wrong or missing (a bare `FS`).
An explicit sign in the note always wins. `examples/byv.gff3` is written this
way.

On the `-` strand write the coordinates as usual (`start < end`); VirPlot
orders the segments in translation direction itself.

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

Any `Note` containing `readthrough` or `read-through` marks the row. **Make
the extension start exactly one base after its partner's end** (`1000` →
`1001`; on the `-` strand, end exactly one base before the partner's start).
That abutment is how VirPlot finds the partner: the extension takes the side
of the ORF it abuts — not merely the previous ORF by coordinate, since a nested
ORF4 can sit between ORF3 and ORF5, as in luteovirids — and never flips. If
nothing abuts, it falls back to the side of the previous ORF, which may be the
wrong one. If instead you give the extension the *full* span of the readthrough
product (starting where the partner starts), that is the RefSeq form above
and is trimmed automatically, provided the `Note` or a `transl_except=` is
present.

### What VirPlot needs, in one table

| To draw | Write | VirPlot then |
|---|---|---|
| a frameshift, from RefSeq | leave the multi-row CDS with `exception=ribosomal slippage` as is | one box per segment, continuation flipped, sign from the coordinates, product labelled once |
| a frameshift, by hand | two CDS rows; on the downstream one `Note=+1 frameshift …` or `Note=-1 frameshift …` | continuation flipped and labelled with the sign you wrote |
| a readthrough, from RefSeq | leave the full-span CDS with `transl_except=` as is | trimmed to the extension, same side as its partner, bar + `RT` |
| a readthrough, by hand | extension row starting at partner end + 1, `Note=readthrough …` | same side as the ORF it abuts, bar + `RT` |

**Layout limit.** The `FS`/`RT` marks are drawn in the linear layout only.
`--layout circular` places these ORFs correctly (a frameshift or readthrough
arc nests like any other) but writes no mark yet.

Subgenomic RNAs are the third expression mechanism; they get their own
section, §12.

## 9. Polyproteins: domains inside one box

Potyviruses, picornaviruses, comoviruses, flaviviruses, coronaviruses — many
genomes translate one long ORF into a polyprotein that is cut into mature
proteins. The ICTV figures draw this as **one box split by thin vertical
lines**, each segment named (Potyviridae Fig. 2: `P1-Pro | HC-Pro | P3 | 6K1 |
CI | …`). VirPlot draws it the same way from the rows RefSeq already provides.

**RefSeq form — nothing to change.** RefSeq lists the mature proteins under
the CDS as `mature_protein_region_of_CDS` rows (the GFF3 spelling of GenBank's
`mat_peptide`), each with `Parent=` naming the CDS and its own `product=`:

```
NC_001616.1	RefSeq	CDS	185	9376	.	+	0	ID=cds-NP_056759.1;product=polyprotein
NC_001616.1	RefSeq	mature_protein_region_of_CDS	185	1036	.	+	.	Parent=cds-NP_056759.1;product=P1 protein
NC_001616.1	RefSeq	mature_protein_region_of_CDS	1037	2404	.	+	.	Parent=cds-NP_056759.1;product=HC-Pro protein
…
```

What VirPlot does with them:

- The CDS is still one feature — one box, one position in the flip/tier
  layout, labelled outside by its `gene=` if it has one. Its own `product=`
  (`polyprotein`) is **not** written: the domains label the box instead.
- Each domain is a segment of the box in its own colour, chosen by the same
  rule as any product (§5: `color_mapping`, then function words, then grey),
  with a thin line at every boundary. Any stretch of the CDS with no domain
  keeps the polyprotein's colour.
- A domain's `product=` is written inside its segment when it fits at the
  figure's size, otherwise just outside the box — the figure's `6K1`, `6K2`.
  Outside labels that would overprint step out onto a second row.
- The last domain usually stops three bases short of the CDS (RefSeq leaves
  the stop codon out); that is not a boundary and no line is drawn there.

**Hand-written form.** Use the `mat_peptide` type if you prefer; `Parent=` is
optional — a row without one is attached to the smallest CDS on its strand
that contains it. Rows that match no CDS are dropped with a warning.

```
A	.	CDS	100	5000	.	+	0	ID=pp;gene=ORF1;product=polyprotein
A	.	mat_peptide	100	2000	.	+	.	product=Pro
A	.	mat_peptide	2001	4997	.	+	.	product=RdRp
```

**Curate the names.** RefSeq's `P1 protein`, `coat protein`, `NIa-VPg
protein` are long for a segment and function-free for the colour rule;
`examples/pvy.gff3` shortens them to the figure's `P1-Pro`, `CP`, `VPg` and
keeps the RefSeq names in `Note=`, and `examples/pvy.yml` pins the figure's
per-domain colours. Do the same for your own genomes.

**Limits.** A domain that straddles a frameshift junction (coronavirus
`nsp12`, which starts in ORF1a and ends in ORF1b) is attached to the segment
holding most of it and drawn clipped to that segment. Domains are drawn in
both layouts; in the circular layout each domain is labelled by the same
rule as an arc (along it, or radially outside when narrow).

## 10. Non-coding features: UTRs, intergenic regions, stem-loops

The ICTV figures mark a few things that are not ORFs: the begomovirus common
region as a thick grey arc with its stem-loop icon (`CRA`, `IR`), the
mastrevirus `LIR`/`SIR`, the hairpin between the two ORFs of an ambisense
segment. VirPlot draws these **on the genome line itself** — they take no
part in the flip/tier/nest layout — in two shapes:

| Shape | Drawn as | For |
|---|---|---|
| region | a grey bar astride the line (linear) or a grey arc astride the circle (circular), named in small text beside it | UTRs, intergenic / common regions, origins |
| stem-loop | a hairpin icon standing on the line (linear) or on the outside of the ring (circular), named at its tip | stem-loops, hairpins, nick sites |

**Which rows are read.** Two tiers, so that a RefSeq file works as is
without its motif annotations cluttering the line:

- Always: `five_prime_UTR`, `three_prime_UTR`, `UTR`, `intergenic_region`,
  `origin_of_replication`, `stem_loop`.
- Only when their words say so: GenBank's catch-all `misc_feature` (which
  NCBI's GFF3 spells `sequence_feature`), `regulatory_region` /
  `regulatory`, `repeat_region`, `sequence_secondary_structure`. These are
  read when the name or `Note=` contains a non-coding word — `intergenic`,
  `common region`, `IR`, `LIR`, `SIR`, `CR`/`CRA`/`CRB`, `UTR`,
  `untranslated`, `stem-loop`, `hairpin`, `non-coding`, `ori` — and skipped
  otherwise (a `misc_feature` reading `RdRp motif` or a `polyA_signal_sequence`
  is not drawn). A short region (under 100 nt) whose words say stem-loop /
  hairpin is drawn as the icon rather than a bar.

**Label.** The first of `Name=`, `product=`, `gene=`, `standard_name=`,
`regulatory_class=`, then `Note=`. The dedicated types have defaults when
nothing is given (`5′ UTR`, `3′ UTR`, `IR`, `ori`); a `stem_loop` row is
labelled only when actually named, since its `Note=` is usually a
description. A region that contains an unnamed stem-loop is named **once, at
the hairpin's tip** — the figures write `CRA` over the icon, not along the
arc. Keep labels short with `Name=`.

**Colour.** Always the `noncoding` grey of the function palette (`#6f6f6f`),
unless the label is pinned in `color_mapping`.

**Circular origin.** By default the hairpin icon marks position 1. A
`stem_loop` row moves it to where the stem-loop actually is — the nick site
is rarely exactly at position 1 in a RefSeq record. Without one, the icon
stays at the origin, unlabelled.

```
# RefSeq form (NCBI GFF3), used as is
NC_001507.1	RefSeq	sequence_feature	2501	326	.	+	.	ID=id-…;Note=common region;gbkey=misc_feature
NC_001507.1	RefSeq	stem_loop	110	142	.	+	.	ID=id-…;Note=conserved stem-loop structure;gbkey=stem_loop

# hand-written, naming things as the figure does
A	.	misc_feature	2501	2926	.	+	.	Name=CRA;Note=common region, spans the origin (2600 + 326)
A	.	stem_loop	110	142	.	+	.	Name=IR
A	.	five_prime_UTR	1	68	.	+	.	ID=utr5
```

An origin-crossing region on a circle is written with `end` past the genome
length (§6), as `examples/grbv.gff3` does for its intergenic region.

## 11. Checklist

- [ ] `##gff-version 3` first line; nine tab-separated columns.
- [ ] One `region` row per molecule, `end` = genome length,
      `Is_circular=true` if it is a circle.
- [ ] `seqid` identical across GFF, SAM/BAM `@SQ`, FASTA header, depth file.
- [ ] One `CDS` row per ORF; `product=` short and naming the **function**;
      `gene=` for the ORF name; long descriptions in `Note=`.
- [ ] Strand `+`/`-` correct — it chooses the layout and the arrow direction.
- [ ] Frameshifts: either leave RefSeq's `exception=ribosomal slippage` rows as
      they are, or split into named ORFs with `Note=+1 frameshift` / `-1` —
      **with the sign written**. Readthrough via `transl_except=`, or an
      extension row that starts at partner end + 1 with `Note=readthrough`.
- [ ] Origin-crossing features only on `Is_circular=true` molecules.
- [ ] Polyproteins: keep RefSeq's `mature_protein_region_of_CDS` rows (or
      write `mat_peptide` rows) and shorten their `product=` to the names you
      want inside the segments.
- [ ] Non-coding landmarks: `five_prime_UTR` / `three_prime_UTR` /
      `stem_loop` rows as RefSeq gives them; an intergenic or common region as
      a `misc_feature` with `Name=IR` (or `CRA`, `LIR`…); a `stem_loop` row
      where the nick site is, so the hairpin icon sits there.
- [ ] Subgenomic RNAs: a `mRNA` (or `transcript`) row per sgRNA with
      `Note=sgRNA` and `gene=` naming what it expresses; `start` is the 5'
      terminus, and the end may be left at the genome length for the usual
      3'-coterminal set.
- [ ] Products with no function word either renamed or pinned in
      `color_mapping`.

## 12. Subgenomic RNAs

Many plus-strand RNA plant viruses express their 3' ORFs from a nested set of
**3'-coterminal subgenomic RNAs** — *Closteroviridae*, *Alphaflexiviridae*,
*Tombusviridae*, *Virgaviridae*. The ICTV figures draw them as shorter lines
stacked beneath the genome, 5' aligned to their start. VirPlot does the same,
longest row first.

GFF3 has no feature type meaning "this is a subgenomic RNA", and RefSeq does
not annotate them for plant viruses, so — exactly as with frameshift and
readthrough (§8) — **this is a curation convention: a transcript row that says
it is one.**

```gff3
# 5' terminus, 3' end at the genome end, named by what it expresses
NC_001598.1	VirPlot	mRNA	13570	15480	.	+	.	ID=sg_cp;gene=CP;Note=sgRNA
```

| Column | Write |
|---|---|
| type | `mRNA`, `transcript`, `ncRNA`, `misc_RNA`, `primary_transcript`, `sequence_feature` or `misc_feature` |
| start | the sgRNA's **5' terminus** — the only coordinate that carries information in a coterminal set |
| end | the genome length (or anything past it) for a 3'-coterminal sgRNA; an earlier end is kept as written, for the 5'-proximal sgRNAs closteroviruses also make |
| `Note=` | must match `sgRNA`, `sg RNA`, `subgenomic` or `sub-genomic` — this marker is **required** |
| `gene=` | what the sgRNA expresses (`CP`, `Hsp70h`); a leading `sgRNA` is stripped, so `gene=sgRNA CP` also labels the row `CP` |

The marker is required rather than inferred from the type, because RefSeq
writes real `mRNA` rows — spliced mastrevirus transcripts, for instance — that
are **not** sgRNAs and must not become ladder rows. A row with no name at all
falls back to `sgRNA1`, `sgRNA2`… numbered 5' to 3'.

`examples/byv.gff3` carries the full closterovirus ladder.

**Why these are annotation and never a depth track.** A plant virus sgRNA is
co-linear with the genome and carries no leader junction (unlike
*Nidovirales*, where `periscope` and `LeTRS` count leader-spanning reads), so
a read from an sgRNA is indistinguishable from a genomic read at the same
coordinate. Per-sgRNA coverage cannot be recovered from short reads. What the
set does leave is a **step in the aggregate depth at each 5' end**, with step
height proportional to that sgRNA's abundance — which is why the rows are
drawn on the same x-axis as the depth trace. A closterovirus at 1× over ORF1a
and 150× over CP is not under-sequenced; it is the expression strategy, and
the ladder is what says so.

**Layout limit.** sgRNA rows are drawn in the **linear layout only**. This is
deliberate, not missing: the ICTV draws no transcript rows on circular
genomes. Geminivirus transcription is bidirectional from the intergenic region
with overlapping transcripts rather than a 3'-coterminal set, and nanovirus
components carry one ORF each. `--layout circular` warns and skips them.
