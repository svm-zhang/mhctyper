"""Allele names retain their fields, spelling, and supported resolution."""

import pytest

from mhctyper.hla_allele import (
    HLAllele,
    HLAllelePattern,
    decompose,
    reduce_resolution,
)


@pytest.mark.parametrize(
    "name, expected",
    [
        ("A*01:01:01", HLAllele("", "A", "01:01:01", "*")),
        ("HLA-c*01:03", HLAllele("HLA-", "c", "01:03", "*")),
        ("HLA-B*04:16N", HLAllele("HLA-", "B", "04:16N", "*")),
        ("hla-c*02:376n", HLAllele("hla-", "c", "02:376n", "*")),
        ("hla_c_07", HLAllele("hla_", "c", "07", "_")),
        ("hla_dqb1_01_01_01", HLAllele("hla_", "dqb1", "01_01_01", "_")),
        (
            "hla_drb1_11_01_01_01",
            HLAllele("hla_", "drb1", "11_01_01_01", "_"),
        ),
    ],
)
def test_decomposes_and_reconstructs_allele_without_changing_identity(name, expected):
    """Class I/II fields, case, separators, and expression suffixes survive parsing."""
    allele = decompose(name, HLAllelePattern())
    assert allele == expected
    assert str(allele) == name


@pytest.mark.parametrize(
    "name, resolution, expected",
    [
        ("A*02:06", 1, "A*02"),
        ("hla_a_01_01_01", 2, "hla_a_01_01"),
        ("hla_drb1*11_01_01:347N", 3, "hla_drb1*11_01_01"),
        ("HLA-dqb1*06:02", 1, "HLA-dqb1*06"),
        ("hla_drb1*11_01_01:347N", -1, "hla_drb1*11"),
    ],
)
def test_reduces_only_fields_beyond_requested_resolution(
    name, resolution, expected, recwarn
):
    """Reduction preserves retained spelling and discards later fields/suffixes."""
    assert reduce_resolution(name, HLAllelePattern(resolution=resolution)) == expected
    assert not recwarn


@pytest.mark.parametrize(
    "name, resolution",
    [
        ("A*02", 1),
        ("hla_a_01_01N", 2),
        ("hla_a_01_01_01", 3),
        ("hla_drb1*11_01_01:347N", 4),
        ("HLA-dqb1*06:02", 3),
        ("hla_drb1_11_01_01_01", 7),
    ],
)
def test_unchanged_resolution_warns_and_preserves_the_complete_name(name, resolution):
    """Already-short names are not padded; resolution is capped at four fields."""
    with pytest.warns(RuntimeWarning, match="No resolution reduced"):
        assert reduce_resolution(name, HLAllelePattern(resolution=resolution)) == name


@pytest.mark.parametrize(
    "name, resolution, expected",
    [
        ("hla_a_01-01-01-01", 2, "hla_a_01-01"),
        ("DQB1*11-01-01-01", 2, "DQB1*11-01"),
        ("HLA-B*11:01-01_01", 3, "HLA-B*11:01-01"),
    ],
)
def test_custom_field_separators_work_for_parsing_and_reduction(
    name, resolution, expected
):
    """A custom delimiter is interpreted without rewriting the name."""
    pattern = HLAllelePattern(digit_field_sep=["-", ":", "_"])
    assert str(decompose(name, pattern)) == name
    reduced_pattern = HLAllelePattern(
        digit_field_sep=["-", ":", "_"], resolution=resolution
    )
    assert reduce_resolution(name, reduced_pattern) == expected


def test_whitespace_around_a_custom_pattern_does_not_change_parsing():
    """Pattern formatting is trimmed before matching a valid class II locus."""
    pattern = HLAllelePattern(locus="\t([A-Z]+[0-9])\n")
    assert decompose("DQB1*06:02", pattern) == HLAllele("", "DQB1", "06:02", "*")


@pytest.mark.parametrize("locus", ["", " \t\n"])
def test_rejects_empty_pattern_fields(locus):
    """A missing locus pattern cannot define an allele parser."""
    with pytest.raises(ValueError):
        HLAllelePattern(locus=locus)


@pytest.mark.parametrize("name", ["", " \t\n", "A*", "11:01:01", "hla_drb1_"])
def test_rejects_names_without_locus_and_numeric_fields(name):
    """Both operations reject inputs lacking the required allele components."""
    pattern = HLAllelePattern()
    with pytest.raises(ValueError):
        decompose(name, pattern)
    with pytest.raises(ValueError):
        reduce_resolution(name, pattern)


@pytest.mark.parametrize("name", ["QQQ*01:01:01", "HLA-ivy*01:03", "hla_drq1_01_01"])
def test_decomposition_rejects_unknown_loci(name):
    """Preserve existing locus validation without adding a new reduction policy."""
    with pytest.raises(ValueError):
        decompose(name, HLAllelePattern())


@pytest.mark.parametrize("name", ["hla_dqb1_02_01_01_", "hla_c*01-01-01"])
def test_decomposition_rejects_incomplete_or_unsupported_field_separators(name):
    """Default parsing requires complete fields using its supported separators."""
    with pytest.raises(ValueError):
        decompose(name, HLAllelePattern())
