"""First-round selection depends on total evidence within each locus."""

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from mhctyper.score_alleles import get_winners


def test_selects_by_total_pair_evidence_independently_per_gene():
    """Two modest pairs beat one stronger pair; class II is scored separately."""
    scores = pl.DataFrame(
        {
            "qnames": ["p1", "p2", "p1", "p1", "p1"],
            "gene": ["hla_a"] * 3 + ["hla_drb1"] * 2,
            "allele": [
                "hla_a_01_01", "hla_a_01_01", "hla_a_02_01",
                "hla_drb1_01_01", "hla_drb1_02_01",
            ],
            "scores": [6.00003, 6.00003, 10.0, 30.0, 40.0],
        }
    )
    expected = pl.DataFrame(
        {
            "allele": ["hla_a_01_01", "hla_drb1_02_01"],
            "gene": ["hla_a", "hla_drb1"],
            # Rounding each pair before summing would incorrectly give 12.0.
            "tot_scores": [12.0001, 40.0],
        }
    )
    assert_frame_equal(get_winners(scores).sort("gene"), expected, check_exact=True)


@pytest.mark.parametrize(
    "second_score, expected_allele, expected_total",
    [
        (10.00001, "hla_a_01_01", 10.0),
        (10.00004, "hla_a_01_01", 10.0),
        (10.00016, "hla_a_02_01", 10.0002),
    ],
    ids=["exact-tie", "tie-after-rounding", "distinct-after-rounding"],
)
def test_rounds_totals_before_lexical_tie_break(
    second_score, expected_allele, expected_total
):
    """Four-decimal ties are deterministic even when input order reverses."""
    scores = pl.DataFrame(
        {
            "allele": ["hla_a_02_01", "hla_a_01_01"],
            "gene": ["hla_a", "hla_a"],
            "scores": [second_score, 10.00001],
        }
    )
    for table in (scores, scores.reverse()):
        assert get_winners(table).to_dicts() == [
            {"allele": expected_allele, "gene": "hla_a", "tot_scores": expected_total}
        ]
