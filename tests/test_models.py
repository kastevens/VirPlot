"""Tests for virplot.models."""

import numpy as np
import pytest

from virplot.models import DepthTrack, Feature, RNA


def test_feature_length_and_midpoint():
    f = Feature(start=100, end=400, strand="+", product="CP")
    assert f.length == 301          # 1-based inclusive
    assert f.midpoint == 250.0


def test_feature_is_immutable():
    f = Feature(start=1, end=2, strand="+", product="x")
    with pytest.raises(Exception):
        f.start = 5  # type: ignore[misc]


def test_rna_positions_are_one_based():
    rna = RNA(name="r", length=5)
    assert rna.positions.tolist() == [1, 2, 3, 4, 5]


def test_rna_add_depth_and_total():
    rna = RNA(name="r", length=4)
    rna.add_depth("a", np.array([1, 2, 3, 4]))
    rna.add_depth("b", np.array([10, 10, 10, 10]))
    assert rna.labels == ["a", "b"]
    assert rna.total_depth().tolist() == [11, 12, 13, 14]
    assert isinstance(rna.depth[0], DepthTrack)


def test_rna_single_track_total_is_the_track():
    rna = RNA(name="r", length=3)
    y = np.array([5, 0, 5])
    rna.add_depth("only", y)
    assert rna.total_depth() is y


def test_rna_duplicate_labels_allowed():
    rna = RNA(name="r", length=2)
    rna.add_depth("s", np.array([1, 1]))
    rna.add_depth("s", np.array([2, 2]))
    assert rna.labels == ["s", "s"]
    assert len(rna.tracks) == 2


def test_rna_add_depth_rejects_wrong_length():
    rna = RNA(name="r", length=3)
    with pytest.raises(ValueError):
        rna.add_depth("bad", np.array([1, 2]))


def test_rna_no_tracks_total_is_zeros():
    rna = RNA(name="r", length=3)
    assert rna.total_depth().tolist() == [0, 0, 0]


def test_rna_defaults_to_linear():
    assert RNA(name="r", length=1).circular is False
