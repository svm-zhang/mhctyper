"""Shared and private pair evidence determine the second allele."""

import polars as pl
import pytest

from mhctyper.score_alleles import get_winners, score_second_by_gene


def test_private_evidence_can_change_the_second_winner():
    """A wins 20:19, but B's private pair makes it win after reweighting."""
    scores = pl.DataFrame(
        {
            "qnames": ["shared", "shared", "private", "shared"],
            "gene": ["hla_a", "hla_a", "hla_a", "hla_b"],
            "allele": ["hla_a_01_01", "hla_a_02_01", "hla_a_02_01", "hla_b_01_01"],
            "scores": [20.0, 12.0, 7.0, 100.0],
        }
    )
    winners = get_winners(scores)
    assert winners.filter(pl.col("gene") == "hla_a")["allele"].to_list() == [
        "hla_a_01_01"
    ]
    winner_pairs = scores.join(winners, on=["gene", "allele"])
    second = score_second_by_gene("hla_a", scores, winner_pairs)

    # The B-locus row shares a query name, but must not join into A evidence.
    assert second["gene"].to_list() == ["hla_a"] * 3
    actual = {
        (row["allele"], row["qnames"]): row["scores"]
        for row in second.to_dicts()
    }
    assert actual == pytest.approx(
        {
            ("hla_a_01_01", "shared"): 10.0,
            ("hla_a_02_01", "shared"): 4.5,
            ("hla_a_02_01", "private"): 7.0,
        }
    )
    assert get_winners(second).to_dicts() == [
        {"allele": "hla_a_02_01", "gene": "hla_a", "tot_scores": 11.5}
    ]


def test_single_supported_candidate_can_win_both_rounds():
    """The first allele stays eligible with half its original pair evidence."""
    scores = pl.DataFrame(
        {
            "qnames": ["p1", "p2"],
            "gene": ["hla_drb1"] * 2,
            "allele": ["hla_drb1_01_01"] * 2,
            "scores": [12.0, 8.0],
        }
    )
    first = get_winners(scores)
    winner_pairs = scores.join(first, on=["gene", "allele"])
    second = score_second_by_gene("hla_drb1", scores, winner_pairs)
    assert second.sort("qnames")["scores"].to_list() == [6.0, 4.0]
    calls = pl.concat([first, get_winners(second)])
    assert calls["allele"].to_list() == ["hla_drb1_01_01"] * 2
    assert calls["tot_scores"].to_list() == [20.0, 10.0]
