"""Core data model: Feature, Domain, Noncoding and RNA.

These are deliberately thin. ``Feature`` is an immutable record for one CDS
(with its polyprotein ``Domain`` segments, if any), ``Noncoding`` one
non-coding landmark (UTR, intergenic region, stem-loop),
``DepthTrack`` one sample's depth array; ``RNA`` bundles everything VirPlot knows about one molecule (its length, its
features, and one depth track per sample) so the rest of the code passes a
single object around instead of parallel lists.

``RNA.circular`` comes from the GFF3 ``Is_circular=true`` attribute on the
region line. It currently drives depth wrapping at the origin; origin-spanning
features and circular plotting will hang off it too.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Domain:
    """One mature protein cut from a polyprotein, in 1-based inclusive bp.

    Comes from a RefSeq ``mature_protein_region_of_CDS`` (GenBank
    ``mat_peptide``) row. Drawn as a segment of its parent CDS's box: its own
    colour, its ``product`` inside, a thin divider at each boundary — the ICTV
    "domains written inside one box, separated by thin vertical lines".
    """

    start: int
    end: int
    product: str


@dataclass(frozen=True)
class Feature:
    """One annotated feature (currently always a CDS) in 1-based, inclusive bp.

    ``product`` is the name written inside the glyph (``RdRp``, ``CP``);
    ``gene`` is the optional ORF name written outside it (``ORF1a``, ``AC1``),
    following the ICTV figure convention.

    ``mechanism`` says the ORF is translated by continuing from its upstream
    neighbour — ``"frameshift"`` (drawn flipped across the line, labelled
    ``+1 FS``/``−1 FS``) or ``"readthrough"`` (kept on the same side behind a
    bar labelled ``RT``). The parser sets it from RefSeq's
    ``exception=ribosomal slippage`` / ``transl_except=`` or from a ``Note=``.

    ``domains`` are the mature proteins of a polyprotein. When present they
    take over the inside of the box (segment colours and labels replace the
    single ``product`` label); ``gene`` is still written outside.
    """

    start: int
    end: int
    strand: str
    product: str
    gene: str | None = None
    mechanism: str | None = None    # "frameshift" | "readthrough": how this ORF is reached
    shift: int | None = None        # +1 / -1 for a frameshift, when known
    show_label: bool = True         # False for a RefSeq join segment that repeats its product
    domains: tuple[Domain, ...] = ()  # mature proteins of a polyprotein, genome order

    @property
    def forward(self) -> bool:
        """True unless the feature is on the complementary strand."""
        return self.strand != "-"

    @property
    def length(self) -> int:
        """Bases covered, for 1-based inclusive coordinates (``Feature(1, 3)`` is 3)."""
        return self.end - self.start + 1

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2

    # --- polyprotein domains ------------------------------------------------

    def domains_within(self, start: int, end: int) -> list[Domain]:
        """Domains clipped to the drawable span ``[start, end]``, in genome order.

        A span is one piece of the feature as drawn (a readthrough extension
        is trimmed, an origin-crossing feature is split), so a domain may be
        cut or fall outside it entirely. A domain edge within one codon of the
        span's end is snapped to it: RefSeq's last mature protein stops short
        of the stop codon, and those three bases are not a segment.
        """
        out = []
        for d in self.domains:
            s, e = max(d.start, start), min(d.end, end)
            if s <= e:
                s = start if s - start <= 3 else s
                e = end if end - e <= 3 else e
                out.append(Domain(s, e, d.product))
        return sorted(out, key=lambda d: d.start)

    def dividers_within(self, start: int, end: int) -> list[float]:
        """Positions of the thin lines between domains inside ``[start, end]``.

        One line wherever a domain begins or ends strictly inside the span —
        so a line also separates a domain from an unannotated stretch — placed
        between bases (``x.5``) and deduplicated.
        """
        cuts: set[float] = set()
        for d in self.domains_within(start, end):
            if d.start > start:
                cuts.add(d.start - 0.5)
            if d.end < end:
                cuts.add(d.end + 0.5)
        return sorted(cuts)


@dataclass(frozen=True)
class Noncoding:
    """A non-coding landmark: UTR, intergenic / common region, stem-loop, origin.

    Not an ORF, so it takes no part in the flip/tier/nest layout; it is drawn
    on the genome line itself, as the ICTV figures do — a grey bar or arc for
    a ``"region"``, a hairpin icon for a ``"stem_loop"`` (Begomovirus ``CRA``
    with its stem-loop, Mastrevirus ``LIR``/``SIR``).
    """

    start: int
    end: int
    strand: str
    label: str
    kind: str = "region"            # "region" | "stem_loop"

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2


@dataclass
class DepthTrack:
    """Per-base read depth for one sample over one RNA."""

    label: str
    values: np.ndarray


@dataclass
class RNA:
    """A single RNA molecule with its annotation and per-sample depth tracks."""

    name: str                       # display / output name
    length: int
    features: list[Feature] = field(default_factory=list)
    depth: list[DepthTrack] = field(default_factory=list)
    seqid: str | None = None        # sequence id in GFF / SAM / depth files
    circular: bool = False
    noncoding: list[Noncoding] = field(default_factory=list)  # UTRs, IRs, stem-loops

    # --- depth tracks -------------------------------------------------------

    def add_depth(self, label: str, y: np.ndarray) -> None:
        """Attach a depth track; the array must have one value per position."""
        if len(y) != self.length:
            raise ValueError(
                f"Depth track {label!r} has {len(y)} positions, expected {self.length}"
            )
        self.depth.append(DepthTrack(label, y))

    @property
    def labels(self) -> list[str]:
        return [t.label for t in self.depth]

    @property
    def tracks(self) -> list[np.ndarray]:
        return [t.values for t in self.depth]

    def total_depth(self) -> np.ndarray:
        """Sum of all depth tracks (the single track if there is only one)."""
        tracks = self.tracks
        if not tracks:
            return np.zeros(self.length, dtype=int)
        return np.sum(tracks, axis=0) if len(tracks) > 1 else tracks[0]

    @property
    def two_strand(self) -> bool:
        """True when ORFs sit on both strands (ambisense RNA, geminiviruses...)."""
        return any(not f.forward for f in self.features) and any(
            f.forward for f in self.features)

    # --- coordinates --------------------------------------------------------

    @property
    def positions(self) -> np.ndarray:
        """1-based genome positions, ``[1, 2, ..., length]``."""
        return np.arange(1, self.length + 1, dtype=int)

    def feature_spans(self, feature) -> list[tuple[int, int]]:
        """Drawable 1-based inclusive spans for a feature (or any ``start``/``end`` record).

        Normally one span. A feature crossing the origin — written either as
        ``start > end`` or with an ``end`` past the genome length, both of
        which occur in the wild — becomes two spans so it can be drawn
        continuously. Coordinates past the length (a padded reference) are
        taken modulo the length. The same arithmetic is applied whether or
        not the molecule is marked circular: on a linear one the two pieces
        are at least drawn where the coordinates say (the parser warns).
        """
        L = self.length
        start, end = feature.start, feature.end
        if start > L:
            start, end = ((start - 1) % L) + 1, end - (start - ((start - 1) % L) - 1)
        if end > L:
            wrapped = end - L
            if wrapped >= start:                   # wraps right round
                return [(1, L)]
            return [(start, L), (1, wrapped)]
        if start > end:
            return [(start, L), (1, end)]
        return [(start, end)]

    def wraps(self, feature) -> bool:
        """True when ``feature`` is written across the origin (two drawable spans)."""
        return len(self.feature_spans(feature)) > 1
