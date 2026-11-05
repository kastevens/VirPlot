# Within-sample variation and divergence from the reference — design notes

Branch `feature/stranded-maf`. Status: **design**. Nothing here changes VirPlot
yet. This note replaces the MAF sketch in `stranded_depth.md` §6.

Every GLRaV-13 number below comes from `docs/design/variation_measurements.py`,
run on the shipped `examples/glrav13.bam` and `examples/glrav13.fa`. The same
script takes any single-reference BAM and the FASTA it was mapped to:

```bash
python3 docs/design/variation_measurements.py                        # GLRaV-13
python3 docs/design/variation_measurements.py --bam S.bam --fasta ref.fa --around 12239
```

---

## 1. Two questions, two comparisons

Mapped reads can disagree with two different things, and the two answer
different questions.

| | Within-sample variation | Divergence from the reference |
|---|---|---|
| Compares | the reads at a site with **each other** | the reads with **the sequence they were mapped to** |
| Question | Is the virus population one isolate, or several? | Is the reference right for this sample? |
| Needs a reference? | No | Yes |
| Typical causes | sequencing and library error; in vivo replication errors; a co-infecting isolate; minority structural variants | a consensus error or misassembly in the reference; a reference from another isolate; a co-infecting isolate that maps only where it is conserved |

Per site, with *n* bases at base quality ≥ 30, allele fractions *f*, the
read majority *M* and the reference base *R*:

* **Non-major fraction** = 1 − *f*(*M*), and **MAF** = the second allele's
  fraction. Both measure within-sample variation.
* **Non-reference fraction** = 1 − *f*(*R*). It measures divergence.
* Where *M* = *R*, the two fractions are equal: the site shows variation,
  not divergence. Where *M* ≠ *R*, the read consensus disagrees with the
  reference, and the non-reference fraction is at least 0.5.

The same per-read evidence (a mismatch, a soft clip) answers either question
depending on **how many reads at the site show it**. A minority means
variation within the sample. A majority means the reference does not
describe this sample. §3.3 has an example of each from one dataset.

---

## 2. Within-sample variation

### 2.1 What makes reads at one site disagree

1. **Sequencing and library error.** About 0.1% per base at Q30. It is
   spread over many sites, it is unlinked, and it piles up at read ends and
   in low-complexity sequence.
2. **In vivo replication errors.** RNA virus replication errors are real
   molecules, the mutant cloud. Each arises on its own genome, so these
   variants are mostly rare, spread across frequencies, and unlinked.
   Replicate libraries confirm them; replicates do not remove them.
3. **A co-infecting isolate.** It puts about (divergence × genome length)
   sites at **one shared frequency**, its share of the population, and
   those sites are **linked** on reads.
4. **Minority structural variants.** Defective RNAs, duplications and
   deletions carried by some molecules. They show as soft-clip pileups
   rather than base changes (GLRaV-13 has one, §3.3).

The single-versus-mixed question is 3 against 1 + 2.

### 2.2 Statistics

| Statistic | What it measures | Verdict for this question | Source |
|---|---|---|---|
| **π** (nucleotide diversity), per site [*N*(*N*−1) − Σ*n*ᵢ(*n*ᵢ−1)] / [*N*(*N*−1)], averaged over sites | total diversity | The right diversity statistic, unbiased with respect to read depth, unlike entropy or variant counts. But it cannot tell a diverse cloud from two clean isolates, and most of it is error (§2.4). | Zhao & Illingworth 2019 |
| **Mean MAF** | total diversity | About π/2, since π ≈ 2·MAF·(1−MAF) per site. Averaged over all sites it is the error rate. Redundant. | — |
| **% of reads with an alternate allele** | anywhere on the read: error rate × read length | Not useful as is. Restricted to variant sites and asked per pair of sites, it becomes the linkage test below. | — |
| **Watterson's θ, Tajima's D** | the shape of the frequency spectrum | The right idea: a mixture gives an excess of intermediate-frequency variants, so D > 0. But D from pooled reads is unreliable, so read the spectrum directly. | Czech et al. (grenedalf) |
| **Fws** = 1 − *H*w/*H*s | within-sample against population heterozygosity | Built for clonal versus mixed (≥ 0.95 "largely clonal"). Needs population allele frequencies from many samples of the same virus. | Auburn et al. 2012 |
| **Frequency clustering** | one cluster of variant frequencies or two (Gaussian mixture, BIC) | The most direct single-sample test. The second cluster's mean is the minor isolate's share. | Sobkowiak et al. 2025 (MixInfect2) |
| **Pooling across defined sites** | the fraction of an isolate's defining mutations present | Detects a minor isolate below the per-site error floor, because errors do not land on dozens of specified sites. Needs known isolates. | Pipek et al. 2024 |
| **Linkage** (*r*² between variant sites on the same reads) | whether minor alleles travel together | Separates an isolate (linked) from errors and replication variants (unlinked). Limited to sites one read or read pair spans. | Cai & Sun 2022 (RVHaplo); Bangratz et al. 2025 |
| **Haplotype reconstruction** | whole haplotypes from long reads | Confirmation. Does not capture the mutant cloud. | Bangratz et al. 2025 |

### 2.3 How the most relevant papers separate an isolate from errors

* **Sobkowiak et al. 2025 (MixInfect2), TB.** Fits one versus two clusters to
  the frequencies at variant sites and compares the fits by BIC. It masks
  sites that vary in pure samples too (within-host evolution), and merges
  variant sites within 100 bp into one data point, by median frequency, so a
  cluster of nearby artefacts cannot dominate the fit. This merge uses
  position only; it does not check whether the alleles share reads.
  Reconstruction of the minor strain worked at 10% and 30% but failed at 5%.
  A clonal bacterium has a quieter background than an RNA virus, so the null
  is noisier for us.
* **Bangratz et al. 2025 (rice yellow mottle virus, nanopore), via RVHaplo
  (Cai & Sun 2022).** The only one of these that puts linkage inside the
  test, twice:
  1. Candidate variant sites come from two binomial tests against a
     sequencing-error rate. A candidate is then kept only if its minor
     allele predicts the allele at some other site on the same reads with
     conditional probability ≥ 0.65. A real haplotype's variants accompany
     each other; an error accompanies nothing.
  2. Reads are joined in a graph by agreement at variant sites, clustered
     (Markov clustering, then hierarchical merging), and clusters holding
     ≤ 0.5% of reads are dropped. More than one surviving cluster means a
     mixture.

  The authors then checked haplotype divergence: real isolates were
  110–192 nt apart, while artefactual haplotypes differed by ≤ 7 nt.
  Mixtures were recovered down to 99/1 for isolates more than 2% apart,
  and the within-host cloud was not captured.
* **Pipek et al. 2024 (SARS-CoV-2).** Calls co-infection when most defining
  mutations of each lineage are present (0.8 of them), with no allele
  frequency filter. Reads are used only afterwards, to look for
  recombination (a recombinant-read ratio threshold of 0.1, PCR chimeras
  treated as a confounder).
* **Roder et al.** The error side only. Most false-positive variants fall
  below 3%. Use 3% for one library, or 1% when replicate libraries agree,
  and ≥ 200× depth.

Combined, the test for one isolate against two is: a cluster of variant
frequencies (MixInfect2), pooled across many sites so it can sit below the
per-site error floor (Pipek), linked on reads where reads allow (RVHaplo),
above an error floor set by replicates (Roder).

What a second isolate would add to GLRaV-13 (17,564 nt; observed π 0.0016):

| Divergence | Sites at one shared frequency | Extra π at a 2% share | at 10% | at 30% |
|---:|---:|---:|---:|---:|
| 0.5% | ~90 | 0.0002 | 0.0009 | 0.0021 |
| 2% | ~350 | 0.0008 | 0.0036 | 0.0084 |
| 5% | ~880 | 0.0020 | 0.0090 | 0.0210 |

π misses a small share: a 2%-divergent isolate at 2% adds half the observed
π. Its 350 sites at one frequency would be unmistakable in the spectrum.

### 2.4 GLRaV-13 H8881

211,830 primary reads; 17,429 sites with ≥ 100 Q30 bases.

| Statistic | Value |
|---|---|
| π | 0.00162, of which 69% comes from sites with MAF < 0.5% |
| Mean MAF | 0.00076, about the Q30 error rate |
| Reads with ≥ 1 non-major base | 7.4%, against 9.5% predicted by 0.1% error over 100 bp |
| Variant sites ≥ 2% | 22 (1.3 per kb) |
| Folded spectrum | 48 sites at 1–2%, none at 2–3%, 16 at 3–10%, 3 at 10–20%, 3 near 40–45% (39.8, 43.7, 45.1%) |
| Linkage, 7 variant pairs one read spans | 4 unlinked (*r*² ≤ 0.09), 2 unmeasurable, 1 linked (552/558, 6 nt apart, *r*² 0.84) |

**Reading.** A single isolate with a few high-frequency within-host variants.
Even a 0.5%-divergent second isolate would put about 90 sites at one
frequency; there are at most three near a shared frequency (1,364, 6,986,
16,372), too far apart for 100 bp single-end reads to phase, and nearby
variants are independent.

### 2.5 Pitfalls

* **The usual strand-bias filter fails on a stranded library.** Read
  orientation is the RNA strand, so "seen on both orientations" means "seen
  on both RNA strands". GLRaV-13 is 94/6, leaving 40–90 reads of the minority
  orientation at most sites. The linked pair 552/558 is about 10% on forward
  reads (− strand RNA) and absent from reverse reads; so is 567. Filter on
  position in the read, base quality and replicate libraries instead.
* **Read-end artefacts.** All 5 reads with the 3.6% variant at 567 have it
  within 10 nt of a read end.
* **The mapping blind spot.** These statistics see only reads that mapped.
  An isolate too divergent to map is among the unmapped reads. This is
  where within-sample variation meets divergence (§3).

---

## 3. Divergence from the reference

### 3.1 What makes reads disagree with the reference

1. **A consensus error in the reference**, if the reference is this
   sample's own assembly: a wrong base or an indel the assembler got wrong.
2. **A misassembly**: a chimeric join of two isolates or two regions, a
   collapsed repeat, an inversion.
3. **A reference from another isolate** ("reference bias"). Reads map less
   well where the sample diverges, so depth and identity track local
   divergence rather than the virus. shiver builds a sample-specific
   reference from contigs for this reason. Across 115 HIV genomes, mapping
   to it instead of the closest of 3,249 existing references recovered a
   median of 205 bases of missing sequence (Wymant et al. 2018).
4. **A co-infecting divergent isolate** whose reads map only where it is
   conserved, adding depth there and mismatches with it.
5. **Biology that changes depth without any disagreement.** Subgenomic RNA
   5′ ends, defective-RNA junctions, and regions that sequence poorly
   (AT-rich, structured).

### 3.2 Signals

| Signal | Points to | Notes |
|---|---|---|
| Sites where the read majority ≠ reference | 1, 3, or 4 where the second isolate is the majority | The decisive per-site test. A base error is one site; a wrong isolate is many. |
| Non-reference fraction along the genome | 3, 4 | Rises across divergent regions; equal to the non-major fraction elsewhere. |
| Mismatches per read along the genome | 2, 3, 4 | From MD/NM tags or the reference. Flat when the reference fits. |
| Soft-clip pileup at one site, with split (SA) alignments | 2, or a minority junction (§2.1, 4) | The clipped parts' alignment says where the junction goes. Its share of the reads at the site says minority or majority. |
| Depth step or dip | any of 1–5 | Never decisive alone; read it with the signals above. |
| Discordant pairs, insert size (paired-end) | 2 | GLRaV-13's reads are single-end. |

Two practical points:

* **Soft clips appear only with a local aligner** (`bwa mem`, or
  `bowtie2 --local`). In bowtie2's default end-to-end mode, which VirPlot's
  `-x` mapping uses, a read across a junction is either rejected or forced
  through with mismatches. Pass `--bowtie2-args "--local"` to see clips.
* **Prior art.** Pilon (Walker et al. 2014) corrects bases, fixes
  misassemblies and fills gaps from read alignments. REAPR (Hunt et al. 2013)
  scores each base by perfect, unique short-insert coverage and finds
  misassemblies from the fragment-coverage distribution of read pairs. Both
  target bacterial and eukaryotic assemblies; the per-site and per-read
  signals carry over to viruses.

### 3.3 GLRaV-13: no divergence, one minority junction, one depth dip

Against its own reference (`glrav13.fa`), GLRaV-13 has:

* **0 sites where the read majority differs from the reference.** The
  non-reference fraction equals the non-major fraction everywhere (0.00085).
* **Mismatches per read: mean 0.14**, flat along the genome (0.11–0.18 per
  2 kb window of read starts). 0.4% of reads have ≥ 3.
* **Soft clips: 1,593 reads (0.8%) clip ≥ 5 nt, but 72 of them clip at one
  site**, right after 9,324. 53 have split alignments, and in 49 of those the
  clipped tail aligns on the same strand from 9,258 (from 9,257 or 9,259 in
  the other 4). Those reads carry a **67-nt tandem duplication of
  9,258–9,324**: an AT-rich repeated motif (27% GC, against 45% for the
  genome) inside the p7 ORF (9,248–9,439). 303 reads align straight
  through, so the duplication is in about **19% of molecules**. That makes
  it within-sample structural variation, not a reference error. Whether it
  is an in vivo variant or a reverse-transcription template switch, the
  reads cannot say; 67 nt is not a multiple of 3, so it would shift p7's
  frame.
* **A depth dip at the same place: 1,364× at 9,241–9,265, 365× at
  9,341–9,365**, recovering to about 650× by 9,400. The read majority
  matches the reference at every site through the dip.

So in one dataset, a depth change, a clip pileup and a correct reference
coexist. A depth change alone does not indicate an assembly error.

---

## 4. Worked case: ENV024, GLRaV-4 strain 6-like

The figure (VirPlot linear layout, 13,834 nt) shows depth near 400–700×
along most of the genome and a sharp rise at the CP ORF. Read from the SVG's
coordinates (to within a few bases):

| Region | Mean depth |
|---|---:|
| ORF1a (241–7,269) | 503× |
| RdRp (7,262–8,815) | 361× |
| HSP70h (8,959–10,560) | 395× |
| p60 (10,541–12,160) | 456× |
| p60–CP gap (12,160–12,239) | 666× |
| CP (12,239–13,048) | 1,709× |
| p23 (13,083–13,706) | 1,781× |

**The step**: 562× → 1,903× (3.4×). It starts (10% of the rise) at 12,216,
about 20 nt upstream of the CP start codon, and finishes (90%) at 12,361, so
it takes about 145 nt, about one read length. (VirPlot's `--smooth` window
is 15 nt, too short to widen it much.) Read 3′→5′, this is the sudden
decrease in depth.

### 4.1 What could make the step, and what each predicts in the reads

| Explanation | Majority ≠ reference | Mismatches per read | Clip pileup at the step | MAF across the step | Stranded library |
|---|---|---|---|---|---|
| **CP subgenomic RNA 5′ end** (biology) | none | flat | none | unchanged | step on the + strand |
| **Consensus error or misassembly upstream** of 12.2 kb | sites on the low side | raised on the low side | likely, at the boundary | unchanged | step on both strands |
| **Reference from another isolate** | many sites, wherever divergent | raised in divergent regions | at indels or junctions | unchanged | both strands |
| **Co-infecting isolate mapping only downstream** | likely on the high side: at 3.4×, the extra reads are ~70% of the depth there | raised on the high side | possibly | **a cluster of linked sites at one shared frequency** on the high side only | both strands |
| **Chimeric join of two isolates** at 12.2 kb | sites on one side | raised on one side | yes, with split alignments across the join | unchanged | both strands |

**The subgenomic RNA explanation is consistent with the shape.**
Closterovirid subgenomic RNAs are 3′-coterminal. In GLRaV-3, a related
ampelovirus, the CP subgenomic RNA is among the abundant ones, and its 5′
end maps to 48 nt upstream of the CP ORF (Jarugula et al. 2010). No read
from a subgenomic RNA can start upstream of its first base, so its coverage
ramps up over about one read length from there: a step of this shape, at
this place. A mapping boundary also ramps over about one read length, so
the shape alone cannot decide. GLRaV-3 does not follow the
3′-genes-are-higher rule (Jarugula et al. 2010), so the step's height is no
guide either.

### 4.2 What would decide it

With the BAM and the FASTA it was mapped to, aligned locally:

```bash
python3 docs/design/variation_measurements.py --bam ENV024.bam --fasta ENV024_ref.fa --around 12239
```

prints depth, the reverse-read share, the non-major fraction, majority ≠
reference sites, clipped reads and mismatches per read in windows either
side of the CP start.

* **No majority ≠ reference sites, flat mismatches, no clip pileup** at the
  step means the reference is right there, and the step is biology: the CP
  subgenomic RNA. If the library is stranded, the step should be on the +
  strand.
* **Disagreement sites or raised mismatches on one side** means the
  reference does not describe the sample there. The side, and whether the
  sites form a cluster at one shared frequency, separate a misassembly from
  a co-infecting isolate.
* **Assemble or BLAST the unmapped reads** in either case, to check for an
  isolate too divergent to map at all.

---

## 5. Design for VirPlot

### 5.1 One pileup pass

`alignments.py` decodes SEQ and QUAL (4-bit BAM sequence; SAM column 10)
alongside the depth it already counts. One pass collects per-site A/C/G/T
counts at base quality ≥ 30, split by read orientation (the stranded-depth
machinery), plus soft-clip starts and ends with their SA partners, and
mismatches per read. Within-sample statistics need no reference.
Divergence needs one: the FASTA given to `-x`, a new `--fasta`, or the
reference bases recovered from MD tags.

### 5.2 Tracks, under the depth track on the same x-axis

* **Variation**: the non-major fraction per site, depth-gated, with the
  threshold line, and variant sites ticked.
* **Divergence**: the non-reference fraction, sites where the majority ≠
  reference marked, clip pileups as ticks labelled with their share of
  reads.

Where the majority matches the reference the two tracks coincide, so the
divergence track shows only what differs. Circular layouts: later, as rings
inside the depth band.

### 5.3 Report (CSV), per sample

* **Variation**: π over all sites and over variant sites only; variant
  sites per kb at the threshold; the folded spectrum in fixed bins; a
  one- versus two-cluster fit with the minor share and the number of sites
  that agree with it; *r*² for variant pairs a read spans.
* **Divergence**: disagreement sites (position, reference → majority,
  fraction); clip pileups (position, reads, share, SA partner); mismatches
  per read by window.

### 5.4 Options (sketch)

| Option | Default |
|---|---|
| `--variation` | off: draw the variation track and write its report |
| `--divergence` | off: draw the divergence track; needs `--fasta`, `-x`, or MD tags |
| `--min-bq` | 30 |
| `--min-depth` | 200 (Roder et al.) |
| `--maf-min` | 0.03; 0.01 when replicate libraries agree (Roder et al.) |

### 5.5 Open questions

* Indels as alleles: count them, or leave them to the divergence track?
* Replicate libraries: a concordance option, or leave replicates to the
  user?
* A population panel for Fws, once FPS has many samples of one virus.
* Paired-end linkage, to phase sites farther apart than one read.

---

## 6. Order of work

1. [ ] Decode SEQ and QUAL; per-site counts by orientation; clips and SA
   partners; mismatches per read. Test against `variation_measurements.py`.
2. [ ] Variation statistics and report.
3. [ ] Variation track.
4. [ ] Divergence statistics, track and report (`--fasta`, MD fallback).
5. [ ] Cluster fit and linkage in the report.
6. [ ] Circular rings.

---

## 7. Sources

* Zhao L, Illingworth CJR (2019). Measurements of intrahost viral diversity
  require an unbiased diversity metric. *Virus Evolution*.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6354029>
* Auburn S et al. (2012). Characterization of within-host *Plasmodium
  falciparum* diversity using next-generation sequence data. *PLoS ONE*.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC3290604/>
* Sobkowiak et al. (2025). MixInfect2, a method for detecting mixed
  *Mycobacterium tuberculosis* infection and reconstructing the strains.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11771024>
* Czech L, Spence JP, Expósito-Alonso M. grenedalf: population genetic
  statistics for the next generation of pool sequencing. *arXiv* 2306.11622.
  <https://arxiv.org/html/2306.11622v3>
* Roder A et al. (2022). Optimized quantification of intrahost viral
  diversity in SARS-CoV-2 and influenza virus sequence data. *bioRxiv*
  10.1101/2021.05.05.442873.
  <https://www.biorxiv.org/content/10.1101/2021.05.05.442873.full.pdf>
* Bangratz et al. (2025). Mixed infections of rice yellow mottle virus
  resolved by nanopore sequencing and RVHaplo.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11731864/>
* Cai D, Sun Y (2022). Reconstructing viral haplotypes using long reads.
  *Bioinformatics* 38(8):2127–2134.
  <https://academic.oup.com/bioinformatics/article/38/8/2127/6528317>
* Pipek et al. (2024). Systematic detection of co-infection and
  intra-host recombination in more than 2 million global SARS-CoV-2
  samples. *Nature Communications*.
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10789779/>
* Walker BJ et al. (2014). Pilon: an integrated tool for comprehensive
  microbial variant detection and genome assembly improvement. *PLoS ONE*
  9(11):e112963. <https://pubmed.ncbi.nlm.nih.gov/25409509/>
* Hunt M et al. (2013). REAPR: a universal tool for genome assembly
  evaluation. *Genome Biology* 14(5):R47.
  <https://www.sanger.ac.uk/tool/reapr/>
* Wymant C et al. (2018). Easy and accurate reconstruction of whole HIV
  genomes from short-read sequence data with shiver. *Virus Evolution*
  4(1). <https://edoc.rki.de/handle/176904/5755>
* Jarugula S, Gowda S, Dawson WO, Naidu RA (2010). 3′-coterminal subgenomic
  RNAs and putative cis-acting elements of *Grapevine leafroll-associated
  virus 3* reveals "unique" features of gene expression strategy in the
  genus *Ampelovirus*. *Virology Journal* 7:180.
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922190/>
