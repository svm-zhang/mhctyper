"""A miniature BAM exercises both scoring rounds and the public file contract."""

import math

import polars as pl
from polars.testing import assert_frame_equal

from mhctyper import run_mhctyper


def test_types_class_i_and_ii_and_writes_pair_scores(typing_inputs, tmp_path):
    """Real workers retain eligible evidence and preserve homozygous call rows."""
    bam, freq = typing_inputs
    outdir = tmp_path / "results"
    calls, result_path = run_mhctyper(bam, freq, outdir, min_ecnt=1, nproc=1)
    pair_score = 8 * (23 + math.log(0.999))
    expected = pl.DataFrame(
        {
            "allele": [
                "hla_a_01_01_01", "hla_a_02_01_01",
                "hla_drb1_01_01_01", "hla_drb1_01_01_01",
            ],
            "gene": ["hla_a", "hla_a", "hla_drb1", "hla_drb1"],
            # A ties in round one; A2 wins round two on private evidence.
            "tot_scores": [
                round(2 * pair_score, 4), round(1.5 * pair_score, 4),
                round(pair_score, 4), round(pair_score / 2, 4),
            ],
            "sample": ["sample_1"] * 4,
        }
    )
    keys = ["allele", "tot_scores"]
    assert_frame_equal(
        calls.sort(keys), expected.sort(keys), abs_tol=1e-8, rel_tol=0
    )
    assert result_path == outdir / "sample_1.hlatyping.res.tsv"
    assert_frame_equal(pl.read_csv(result_path, separator="\t"), calls)

    a1 = pl.read_csv(outdir / "sample_1.a1.tsv", separator="\t")
    expected_pairs = pl.DataFrame(
        {
            "qnames": ["shared", "a1_private", "shared", "a2_private", "drb1_private"],
            "scores": [pair_score] * 5,
            "allele": (
                ["hla_a_01_01_01"] * 2 + ["hla_a_02_01_01"] * 2
                + ["hla_drb1_01_01_01"]
            ),
            "gene": ["hla_a"] * 4 + ["hla_drb1"],
        }
    )
    pair_keys = ["allele", "qnames"]
    assert_frame_equal(
        a1.sort(pair_keys), expected_pairs.sort(pair_keys),
        abs_tol=1e-10, rel_tol=0,
    )

    a2 = pl.read_csv(outdir / "sample_1.a2.tsv", separator="\t")
    expected_second = expected_pairs.with_columns(
        pl.Series(
            "scores",
            [pair_score / 2, pair_score / 2, pair_score / 2, pair_score, pair_score / 2],
        )
    )
    assert_frame_equal(
        a2.select(expected_second.columns).sort(pair_keys),
        expected_second.sort(pair_keys), abs_tol=1e-10, rel_tol=0,
    )
