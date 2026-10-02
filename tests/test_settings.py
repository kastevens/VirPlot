"""Tests for virplot.settings."""

import textwrap

from virplot.settings import Settings, load_settings


def test_load_settings_full(tmp_path):
    yml = tmp_path / "spec.yml"
    yml.write_text(textwrap.dedent("""\
        color_mapping:
          RdRp: '#ff0000'
          CP: '#00ff00'
        default_color: '#aaaaaa'
        depth_line_color: '#0000ff'
        shade_color: '#ff6666'
        annotation_fontsize: 12
        stacked_area_colors:
          - '#111111'
          - '#222222'
        legend_location: "lower right"
        title: "My Plot"
    """))
    s = load_settings(str(yml))
    assert s.color_mapping == {"RdRp": "#ff0000", "CP": "#00ff00"}
    assert s.default_color == "#aaaaaa"
    assert s.depth_line_color == "#0000ff"
    assert s.annotation_fontsize == 12
    assert s.legend_location == "lower right"
    assert s.title == "My Plot"
    assert len(s.stacked_area_colors) == 2


def test_load_settings_defaults(tmp_path):
    yml = tmp_path / "empty.yml"
    yml.write_text("")
    s = load_settings(str(yml))
    assert s == Settings()


def test_load_settings_partial(tmp_path):
    yml = tmp_path / "partial.yml"
    yml.write_text("annotation_fontsize: 14\n")
    s = load_settings(str(yml))
    assert s.annotation_fontsize == 14
    assert s.default_color == Settings.default_color


def test_lighter_accepts_any_matplotlib_colour():
    from virplot.settings import lighter
    assert lighter("#000000") == "#595959"
    assert lighter("#abc").startswith("#") and len(lighter("#abc")) == 7
    assert lighter("tomato").startswith("#")          # named colour, no traceback
    assert lighter("not-a-colour") == "not-a-colour"   # left for matplotlib to report


def test_putative_product_with_named_palette_colour_does_not_crash():
    s = Settings(function_palette={"capsid": "orange"})
    assert s.feature_color("putative coat protein").startswith("#")


def test_empty_color_mapping_in_yaml_means_no_mapping(tmp_path):
    p = tmp_path / "s.yml"
    p.write_text("color_mapping:\ntitle: x\n")        # key present, value null
    s = load_settings(str(p))
    assert s.color_mapping == {}
    assert s.feature_color("CP")                      # no TypeError
