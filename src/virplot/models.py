"""Core data model: Feature and RNA.

These are deliberately thin. ``Feature`` is an immutable record for one CDS,
``DepthTrack`` one sample's depth array; ``RNA`` bundles everything VirPlot knows about one molecule (its length, its
features, and one depth track per sample) so the rest of the code passes a
single object around instead of parallel lists.

``RNA.circular`` is a data-only flag for now: no behaviour hangs off it yet.
It marks where topology-dependent logic (coordinate wrapping, origin-spanning
features, circular plotting) will attach later.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Feature:
    """One annotated feature (currently always a CDS) in 1-based, inclusive bp."""

    start: int
    end: int
    strand: str
    product: str

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

    name: str
    length: int
    features: list[Feature] = field(default_factory=list)
    depth: list[DepthTrack] = field(default_factory=list)
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

    # --- coordinates --------------------------------------------------------

    @property
    def positions(self) -> np.ndarray:
        """1-based genome positions, ``[1, 2, ..., length]``."""
        return np.arange(1, self.length + 1, dtype=int)
