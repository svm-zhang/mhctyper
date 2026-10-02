"""Resource selection preserves allele identity and sample provenance."""

from types import SimpleNamespace

import pytest

from mhctyper.utils import (
    collect_alleles_to_type,
    load_allele_pop_freq,
    load_rg_sm_from_bam,
)


def test_frequency_filter_keeps_fractional_support_and_full_allele_names(tmp_path):
    """Only eligible families survive; distinct references are never collapsed."""
    freq = tmp_path / "freq.tsv"
    freq.write_text(
        "Allele\tlocal_population\tother_population\n"
        "hla_a_01_01\t0.0001\t0\n"
        "hla_a_01_01\t0.2\t0\n"  # Repeated family must not multiply candidates.
        "hla_a_02_01\t0\t0\n"
        "hla_drb1_01_01\t0\t0.3\n"
    )
    eligible = load_allele_pop_freq(freq)
    assert eligible["Allele"].to_list() == [
        "hla_a_01_01", "hla_a_01_01", "hla_drb1_01_01"
    ]
    names = [
        "hla_a_01_01_01", "hla_a_01_01_02N", "hla_a_02_01_01",
        "hla_a_03_01_01", "hla_drb1_01_01_01",
    ]
    metadata = SimpleNamespace(seqnames=lambda: names)
    assert collect_alleles_to_type(metadata, eligible["Allele"].to_list()) == [
        "hla_a_01_01_01", "hla_a_01_01_02N", "hla_drb1_01_01_01"
    ]


def test_reduced_names_match_frequency_entries_exactly():
    """A suffix on a retained field remains part of the eligibility key."""
    metadata = SimpleNamespace(
        seqnames=lambda: ["hla_a_01_01_01N", "hla_a_01_01N"]
    )
    assert collect_alleles_to_type(metadata, ["hla_a_01_01"]) == ["hla_a_01_01_01N"]
    assert collect_alleles_to_type(metadata, ["hla_a_01_01N"]) == ["hla_a_01_01N"]


def test_incompatible_resources_fail_instead_of_typing_unlisted_alleles():
    """A nonempty BAM reference list still needs an eligible family."""
    metadata = SimpleNamespace(seqnames=lambda: ["hla_a_01_01_01"])
    with pytest.raises(SystemExit) as error:
        collect_alleles_to_type(metadata, ["hla_b_01_01"])
    assert error.value.code == 1


def test_reads_sample_name_from_single_read_group():
    """Output provenance uses SM rather than the read-group ID."""
    metadata = SimpleNamespace(read_groups=[{"ID": "lane_1", "SM": "sample_1"}])
    assert load_rg_sm_from_bam(metadata) == "sample_1"


@pytest.mark.parametrize(
    "groups",
    [
        [{"ID": "lane_1"}],
        [{"ID": "lane_1", "SM": ""}],
        [{"ID": "lane_1", "SM": "sample"}, {"ID": "lane_2", "SM": "sample"}],
    ],
    ids=["missing-sample", "empty-sample", "multiple-groups-same-sample"],
)
def test_rejects_unusable_sample_metadata(groups):
    """Missing SM and multiple read groups follow the existing handled error path."""
    with pytest.raises(SystemExit) as error:
        load_rg_sm_from_bam(SimpleNamespace(read_groups=groups))
    assert error.value.code == 1
