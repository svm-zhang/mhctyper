"""Small numerical examples with independently calculated probabilities."""

import math

import pytest

from mhctyper.score_alleles import score_log_liklihood


def test_scores_matches_from_each_base_quality():
    """Every matched base contributes its own probability and scale term."""
    actual = score_log_liklihood({"bqs": [10, 20, 30], "mds": ["3"]})
    expected = 69 + math.log(0.9) + math.log(0.99) + math.log(0.999)
    assert actual == pytest.approx(expected, abs=1e-10, rel=0)


@pytest.mark.parametrize(
    "mds, probabilities",
    [
        (["0", "C", "2"], [0.1 / 3, 0.99, 0.999]),
        (["1", "C", "1"], [0.9, 0.01 / 3, 0.999]),
        (["2", "C", "0"], [0.9, 0.99, 0.001 / 3]),
    ],
    ids=["first-base", "middle-base", "last-base"],
)
def test_mismatch_uses_quality_at_its_aligned_position(mds, probabilities):
    """Unequal qualities expose a shifted MD-to-quality index."""
    actual = score_log_liklihood({"bqs": [10, 20, 30], "mds": mds})
    expected = 69 + sum(math.log(p) for p in probabilities)
    assert actual == pytest.approx(expected, abs=1e-10, rel=0)


def test_zero_length_match_block_does_not_advance_quality_index():
    """Adjacent mismatches separated by MD zero use consecutive qualities."""
    actual = score_log_liklihood(
        {"bqs": [10, 20, 30], "mds": ["0", "C", "0", "G", "1"]}
    )
    expected = 69 + math.log(0.1 / 3) + math.log(0.01 / 3) + math.log(0.999)
    assert actual == pytest.approx(expected, abs=1e-10, rel=0)
