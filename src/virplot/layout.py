"""Where each feature glyph goes: side of the backbone and distance from it.

Implements the placement rules in docs/ictv_drawing_conventions.md:

* ``flip`` (mode L1, one strand): the 5'-most ORF sits above the line; an ORF
  that overlaps its upstream neighbour flips to the other side, otherwise it
  stays on the side of that neighbour.  A caller may add a second trigger
  through ``flip_if(prev, feat)`` — the plotter uses it to flip a glyph that
  is close to, and the same colour as, its neighbour, which is what the ICTV
  figures do for BYV's CPm/CP so two flat boxes do not read as one.  If the
  chosen side is already taken where it would land — a long upstream ORF can
  do this — the other side is tried, and failing that the glyph tiers
  outward.  The document's rule never leaves two boxes on top of each other
  this way.
* ``tier`` (mode L2, both strands): side is fixed by strand (+ above, − below)
  and same-side overlap is resolved by tiering further from the line.
* ``nest`` (modes C1/C2, circular): arcs are placed largest first on the
  outermost lane, and an arc overlapping one already placed steps one lane
  inward — so a small ORF sits inside the large one it overlaps.

Pure functions over ``Feature`` objects, so the rules are testable without a
figure.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

from virplot.models import Feature

Spans = list[tuple[int, int]]
SpanFn = Callable[[Feature], Spans]

ABOVE, BELOW = 1, -1


@dataclass(frozen=True)
class Placement:
    feature: Feature
    side: int       # ABOVE / BELOW the backbone (linear); always ABOVE for circular
    tier: int       # 0 touches the backbone; higher = further from it (or further inward)


def spans_overlap(a: Spans, b: Spans) -> bool:
    """True if any interval of ``a`` intersects any interval of ``b`` (inclusive)."""
    return any(s1 <= e2 and s2 <= e1 for s1, e1 in a for s2, e2 in b)


def resolve_mode(mode: str, two_strand: bool) -> str:
    """Turn the spec.yml setting into a concrete rule for a linear track."""
    if mode == "auto":
        return "tier" if two_strand else "flip"
    return mode


FlipIf = Callable[[Feature, Feature], bool]


def place_features(features: Iterable[Feature], spans: SpanFn, mode: str,
                   flip_if: FlipIf | None = None) -> list[Placement]:
    """Assign side and tier to every feature under ``mode`` (flip / tier / nest).

    ``flip_if(prev, feat)`` is an extra reason to flip in ``flip`` mode, judged
    between each feature and its upstream neighbour; ignored by other modes.
    """
    ordered = sorted(features, key=lambda f: (f.start, f.end))
    if mode == "flip":
        return _place_flip(ordered, spans, flip_if)
    if mode == "tier":
        return _place_tier(ordered, spans, lambda f: ABOVE if f.forward else BELOW)
    if mode == "nest":
        # largest arcs claim the outer lane so a small ORF nests inside the big
        # one it overlaps (AC4 inside AC1), whatever their start coordinates
        covered = lambda f: sum(e - s + 1 for s, e in spans(f))
        by_size = sorted(features, key=lambda f: (-covered(f), f.start))
        return _place_tier(by_size, spans, lambda f: ABOVE)
    raise ValueError(f"unknown layout mode {mode!r}")


# --------------------------------------------------------------------------

def _free_tier(placed: list[Placement], spans: SpanFn, feat: Feature, side: int) -> int:
    """Lowest tier on ``side`` where ``feat`` overlaps nothing already there."""
    tier = 0
    while any(p.side == side and p.tier == tier
              and spans_overlap(spans(p.feature), spans(feat)) for p in placed):
        tier += 1
    return tier


def _place_tier(ordered: list[Feature], spans: SpanFn,
                side_of: Callable[[Feature], int]) -> list[Placement]:
    placed: list[Placement] = []
    for feat in ordered:
        side = side_of(feat)
        placed.append(Placement(feat, side, _free_tier(placed, spans, feat, side)))
    return placed


def _place_flip(ordered: list[Feature], spans: SpanFn,
                flip_if: FlipIf | None = None) -> list[Placement]:
    placed: list[Placement] = []
    current = ABOVE
    prev: Feature | None = None
    for feat in ordered:
        wanted = current
        if prev is not None and (spans_overlap(spans(prev), spans(feat))
                                 or (flip_if is not None and flip_if(prev, feat))):
            wanted = -current
        # the document's rule, with a safety net for a long upstream ORF
        if _free_tier(placed, spans, feat, wanted) == 0:
            side, tier = wanted, 0
        elif _free_tier(placed, spans, feat, -wanted) == 0:
            side, tier = -wanted, 0
        else:
            side, tier = wanted, _free_tier(placed, spans, feat, wanted)
        placed.append(Placement(feat, side, tier))
        current, prev = side, feat
    return placed
