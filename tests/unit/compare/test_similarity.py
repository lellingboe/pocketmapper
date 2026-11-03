"""
Unit tests for `pocketmapper.compare.similarity`: X in the normalising sequence is skipped, and scores stay in -1..1.
"""

import math
from itertools import product

import pytest

from pocketmapper.compare.similarity import BLOSUM62, binary_similarity, full_similarity

NAN = math.nan

FULL_TABLE = [
    pytest.param("ACDW", "ACDW", 1.0, id="identical"),
    pytest.param("X", "W", NAN, id="x-in-a"),
    pytest.param("XW", "WW", 1.0, id="x-in-a-skipped"),
    pytest.param("WW", "WX", (1 + BLOSUM62["W", "X"] / BLOSUM62["W", "W"]) / 2, id="x-in-b-scored"),
    pytest.param("acdw", "ACDW", 1.0, id="lowercase"),
    pytest.param("OW", "WW", 1.0, id="o-as-x"),
    pytest.param("UW", "WW", 1.0, id="u-as-x"),
    pytest.param("-W", "WW", 1.0, id="gap-as-x"),
    pytest.param("*W", "WW", 1.0, id="stop-as-x"),
    pytest.param("WW", "WO", (1 + BLOSUM62["W", "X"] / BLOSUM62["W", "W"]) / 2, id="o-in-b-as-x"),
    pytest.param("XXX", "ACD", NAN, id="all-x"),
    pytest.param("B", "Z", BLOSUM62["B", "Z"] / BLOSUM62["B", "B"], id="b-z-biopython"),
]

BINARY_TABLE = [
    pytest.param("ACDW", "ACDW", 1.0, id="identical"),
    pytest.param("X", "W", NAN, id="x-in-a"),
    pytest.param("XW", "WW", 1.0, id="x-in-a-skipped"),
    pytest.param("WW", "WX", 0.5, id="x-in-b-scored"),
    pytest.param("acdw", "ACDW", 1.0, id="lowercase"),
    pytest.param("OW", "WW", 1.0, id="o-as-x"),
    pytest.param("UW", "WW", 1.0, id="u-as-x"),
    pytest.param("-W", "WW", 1.0, id="gap-as-x"),
    pytest.param("XXX", "ACD", NAN, id="all-x"),
    # Biopython scores B-Z as 1; the old bundled matrix scored it 0
    pytest.param("B", "Z", 1.0, id="b-z-biopython"),
]


def assert_score(actual, expected):
    """Compare two scores, treating NaN as equal to NaN."""
    if math.isnan(expected):
        assert math.isnan(actual)
    else:
        assert actual == pytest.approx(expected)


@pytest.mark.parametrize(("seq_a", "seq_b", "expected"), FULL_TABLE)
def test_full_similarity(seq_a, seq_b, expected):
    """`full_similarity` normalises by seq_a, skipping its X positions."""
    assert_score(full_similarity(seq_a, seq_b), expected)


@pytest.mark.parametrize(("seq_a", "seq_b", "expected"), BINARY_TABLE)
def test_binary_similarity(seq_a, seq_b, expected):
    """`binary_similarity` skips the same positions as `full_similarity`."""
    assert_score(binary_similarity(seq_a, seq_b), expected)


def test_full_similarity_bounded():
    """No letter pair in the matrix alphabet scores outside -1..1."""
    for a, b in product(BLOSUM62.alphabet, repeat=2):
        score = full_similarity(a, b)
        assert math.isnan(score) or -1 <= score <= 1, (a, b, score)
