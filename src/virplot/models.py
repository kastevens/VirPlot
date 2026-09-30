"""Core data model: Feature and RNA.

These are deliberately thin. ``Feature`` is an immutable record for one CDS,
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
    """

    start: int
    end: int
    strand: str
    product: str
    gene: str | None = None
    mechanism: str | None = None    # "frameshift" | "readthrough": how this ORF is reached
    shift: int | None = None        # +1 / -1 for a frameshift, when known
    show_label: bool = True         # False for a RefSeq join segment that repeats its product

    @property
    def forward(self) -> bool:
        """True unless the feature is on the complementary strand."""
        return self.strand != "-"

    @property
    def length(self) -> int:
        return self.end - self.start

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

    def feature_spans(self, feature: Feature) -> list[tuple[int, int]]:
        """Drawable 1-based inclusive spans for a feature.

        Normally one span. On a circular molecule a feature crossing the
        origin — written either as ``start > end`` or with an ``end`` past the
        genome length, both of which occur in the wild — becomes two spans so
        it can be drawn continuously.
        """
        start, end = feature.start, feature.end
        if not self.circular:
            return [(start, end)]
        if end > self.length:
            wrapped = end - self.length
            if wrapped >= start:                   # wraps right round
                return [(1, self.length)]
            return [(start, self.length), (1, wrapped)]
        if start > end:
            return [(start, self.length), (1, end)]
        return [(start, end)]
