"""YAML settings loader and the default function palette."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

import yaml
from matplotlib.colors import to_hex, to_rgb

log = logging.getLogger(__name__)

_KNOWN_KEYS = {
    "color_mapping", "default_color", "depth_line_color", "shade_color",
    "annotation_fontsize", "stacked_area_colors", "legend_location", "title",
    "overlap_mode", "end_5_label", "end_3_label", "function_palette",
    "circular_labels", "circular_arcs", "sgrna_color",
}

CIRCULAR_LABELS = ("tangential", "horizontal")
CIRCULAR_ARCS = ("outside", "on_circle")
LAYOUTS = ("auto", "linear", "circular")

OVERLAP_MODES = ("auto", "flip", "tier")

# ICTV 9th Report family-figure palette: colour encodes predicted function,
# never strand or reading frame. Keys are function classes; the classifier
# below maps product names onto them. Override any class in spec.yml under
# ``function_palette``, or pin a single product under ``color_mapping``.
DEFAULT_FUNCTION_PALETTE = {
    "replicase": "#f5b041",      # Mtr / Hel / RdRp / Rep — yellow-orange
    "capsid": "#e75480",         # CP, CPm — magenta-pink
    "movement": "#3ea6b5",       # MP — teal/blue
    "hsp70": "#7dbb4f",          # HSP70h — green
    "suppressor": "#9b6fc4",     # silencing suppressors / small 3' ORFs — purple
    "noncoding": "#6f6f6f",      # IR / UTR / stem-loop — grey or black
    "unknown": "#9f9f9f",
}

# ordered: first match wins, so "coat protein" beats the "protein" catch-all
_FUNCTION_RULES: list[tuple[str, str]] = [
    (r"\b(rdrp|rna[- ]dependent|polymerase|replicase|replication[- ]associated|"
     r"\brep\b|repa\b|methyltransferase|mtr\b|helicase|\bhel\b|protease|"
     r"papain|nsp\d|\bp1ab?\b|orf1[ab]?\b|polyprotein)", "replicase"),
    (r"\b(hsp70|hsp90|heat[- ]shock)", "hsp70"),
    (r"\b(coat|capsid|nucleocapsid|cp\b|cpm\b|cph\b|\bn protein|virion protein|"
     r"structural)", "capsid"),
    (r"\b(movement|\bmp\b|cell[- ]to[- ]cell|triple gene block|tgb\d?)", "movement"),
    # named function only: a bare "p20"/"p6" says nothing about role (BYV p6 is
    # a membrane protein and grey in the figure), so pN names fall to default
    (r"\b(silencing|suppressor|vsr|hc-?pro)", "suppressor"),
    (r"\b(intergenic|\bir\b|\b[ls]ir\b|common[- ]region|\bcr[ab]?\b|utr|untranslated|"
     r"stem[- ]loop|hairpin|non-?coding|origin of replication|\bori\b)", "noncoding"),
]

_PUTATIVE = re.compile(r"\b(putative|probable|possible|hypothetical|predicted|"
                       r"proposed|-like)\b", re.I)


def classify_function(product: str) -> str | None:
    """Map a product name onto a function class, or None if nothing matches."""
    text = product.lower()
    for pattern, cls in _FUNCTION_RULES:
        if re.search(pattern, text):
            return cls
    return None


def lighter(color: str, amount: float = 0.35) -> str:
    """Tint any matplotlib colour (``#rrggbb``, ``#abc``, ``tomato``) towards
    white (ICTV: 'putative' = lighter); an unparsable colour is returned as is
    so matplotlib reports it at draw time."""
    try:
        r, g, b = to_rgb(color)
    except ValueError:
        return color
    mix = lambda c: c + (1.0 - c) * amount
    return to_hex((mix(r), mix(g), mix(b)))


@dataclass
class Settings:
    color_mapping: dict[str, str] = field(default_factory=dict)
    default_color: str = "#9F9F9F"
    depth_line_color: str = "#85dbec"        # original figure's top layer
    shade_color: str = "tomato"
    annotation_fontsize: int = 8
    stacked_area_colors: list[str] = field(  # original figure's layers, bottom first
        default_factory=lambda: ["#567eb0", "#55b9d9", "#85dbec"])
    legend_location: str = "upper left"
    title: str = ""
    overlap_mode: str = "auto"           # auto | flip | tier
    end_5_label: str = "5'"              # e.g. "5' m7G", "VPg"
    end_3_label: str = "3'"              # e.g. "3' OH", "A(n)"
    function_palette: dict[str, str] = field(
        default_factory=lambda: dict(DEFAULT_FUNCTION_PALETTE))
    # circular layout style (ICTV figures: horizontal labels, arcs on the circle)
    circular_labels: str = "tangential"      # tangential | horizontal
    circular_arcs: str = "outside"           # outside | on_circle
    sgrna_color: str = "#8C2F39"             # subgenomic RNA rows under the line
    # how this genome should be drawn; --layout overrides it. "auto" draws a
    # circular molecule as a circle and a linear one as a track.
    layout: str = "auto"                     # auto | linear | circular

    def feature_color(self, product: str) -> str:
        """Colour for a product: explicit mapping, then function class, then default.

        A product whose name marks it as putative gets a lighter tint of its
        class colour, as in the ICTV figures.
        """
        if product in self.color_mapping:
            return self.color_mapping[product]
        cls = classify_function(product)
        if cls is None:
            return self.default_color
        color = self.function_palette.get(cls, DEFAULT_FUNCTION_PALETTE[cls])
        return lighter(color) if _PUTATIVE.search(product) else color


def load_settings(yaml_path: str) -> Settings:
    """Load a YAML spec file and return a Settings instance."""
    with open(yaml_path) as fp:
        raw = yaml.safe_load(fp) or {}

    unknown = set(raw) - _KNOWN_KEYS
    if unknown:
        log.warning("Unknown YAML keys ignored: %s", ", ".join(sorted(unknown)))

    missing = _KNOWN_KEYS - set(raw)
    if missing:
        log.info("Using defaults for missing YAML keys: %s", ", ".join(sorted(missing)))

    overlap_mode = str(raw.get("overlap_mode", Settings.overlap_mode)).lower()
    if overlap_mode not in OVERLAP_MODES:
        log.warning("overlap_mode %r not one of %s; using 'auto'",
                    overlap_mode, "/".join(OVERLAP_MODES))
        overlap_mode = "auto"

    palette = dict(DEFAULT_FUNCTION_PALETTE)
    palette.update(raw.get("function_palette") or {})

    def choice(key: str, options: tuple, default: str) -> str:
        value = str(raw.get(key, default)).lower()
        if value not in options:
            log.warning("%s %r not one of %s; using %r", key, value, "/".join(options), default)
            return default
        return value

    return Settings(
        color_mapping=dict(raw.get("color_mapping") or {}),
        default_color=raw.get("default_color", Settings.default_color),
        depth_line_color=raw.get("depth_line_color", Settings.depth_line_color),
        shade_color=raw.get("shade_color", Settings.shade_color),
        annotation_fontsize=raw.get("annotation_fontsize", Settings.annotation_fontsize),
        stacked_area_colors=list(raw.get("stacked_area_colors") or Settings().stacked_area_colors),
        legend_location=raw.get("legend_location", Settings.legend_location),
        title=raw.get("title", Settings.title),
        overlap_mode=overlap_mode,
        end_5_label=str(raw.get("end_5_label", Settings.end_5_label)),
        end_3_label=str(raw.get("end_3_label", Settings.end_3_label)),
        function_palette=palette,
        circular_labels=choice("circular_labels", CIRCULAR_LABELS, Settings.circular_labels),
        circular_arcs=choice("circular_arcs", CIRCULAR_ARCS, Settings.circular_arcs),
        sgrna_color=str(raw.get("sgrna_color", Settings.sgrna_color)),
        layout=choice("layout", LAYOUTS, Settings.layout),
    )
