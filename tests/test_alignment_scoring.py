"""Real BAM records exercise mhctyper's reader and pair-filtering boundary."""

import math

import pytest

from mhctyper.score_alleles import score_per_allele

from .bam_helpers import make_pair, write_bam


ALLELE = "hla_a_01_01_01"
PAIR_SCORE = 8 * (23 + math.log(0.999))


def test_sums_mates_only_for_the_requested_allele(tmp_path):
    """An alternative-reference pair cannot enter the target's evidence table."""
    bam = write_bam(
        tmp_path / "reads.bam", [ALLELE, "hla_b_01_01_01"],
        make_pair("target", 0) + make_pair("other_reference", 1),
    )
    result = score_per_allele(ALLELE, bam, min_ecnt=1)
    assert result is not None
    assert result.select("qnames", "allele", "gene").to_dicts() == [
        {"qnames": "target", "allele": ALLELE, "gene": "hla_a"}
    ]
    assert result["scores"].item() == pytest.approx(PAIR_SCORE, abs=1e-10, rel=0)


def test_keeps_exact_mismatch_threshold_and_drops_pair_above_it(tmp_path):
    """A failing mate removes its otherwise-good partner, not just itself."""
    boundary = make_pair("boundary")
    boundary[0].set_tag("MD", "1C2")
    boundary[0].set_tag("NM", 1)
    # A C at reference position 11 makes the boundary read's MD consistent.
    # Put the second pair elsewhere, against two C reference positions.
    above = make_pair("above")
    above[0].reference_start = 50
    above[1].reference_start = 70
    above[0].next_reference_start = 70
    above[1].next_reference_start = 50
    above[0].set_tag("MD", "0C0C2")
    above[0].set_tag("NM", 2)
    bam = write_bam(tmp_path / "reads.bam", [ALLELE], boundary + above)
    result = score_per_allele(ALLELE, bam, min_ecnt=1)
    assert result is not None
    assert result["qnames"].to_list() == ["boundary"]
    expected = 8 * 23 + 7 * math.log(0.999) + math.log(0.001 / 3)
    assert result["scores"].item() == pytest.approx(expected, abs=1e-10, rel=0)


@pytest.mark.parametrize(
    "cigar, sequence, md",
    [("2M1I2M", "AAAAA", "4"), ("2M1D2M", "AAAA", "2^A2")],
    ids=["insertion", "deletion"],
)
def test_indel_in_one_mate_excludes_the_whole_pair(tmp_path, cigar, sequence, md):
    """An indel-free partner alone must not contribute evidence."""
    pair = make_pair("indel")
    pair[0].query_sequence = sequence
    pair[0].query_qualities = [30] * len(sequence)
    pair[0].cigarstring = cigar
    pair[0].set_tag("MD", md)
    pair[0].set_tag("NM", 1)
    bam = write_bam(tmp_path / "reads.bam", [ALLELE], pair)
    assert score_per_allele(ALLELE, bam, min_ecnt=1) is None


@pytest.mark.parametrize(
    "flag", [512, 1024, 2048], ids=["qc-fail", "duplicate", "supplementary"]
)
def test_excluded_flag_cannot_leave_a_scored_orphan(tmp_path, flag):
    """Flag exclusions apply before the surviving-pair count."""
    excluded = make_pair("excluded")
    excluded[0].flag |= flag
    bam = write_bam(tmp_path / "reads.bam", [ALLELE], excluded + make_pair("retained"))
    result = score_per_allele(ALLELE, bam, min_ecnt=1)
    assert result is not None
    assert result["qnames"].to_list() == ["retained"]
    assert result["scores"].item() == pytest.approx(PAIR_SCORE)


def test_improper_extra_alignment_is_removed_before_pair_counting(tmp_path):
    """An extra singleton does not invalidate two properly aligned mates."""
    extra = make_pair("retained")[0]
    extra.flag = 1 | 8 | 64 | 256  # Mapped read1, mate unmapped, secondary, not proper.
    extra.next_reference_id = -1
    extra.next_reference_start = -1
    extra.template_length = 0
    improper = make_pair("improper")
    for read in improper:
        read.flag &= ~2
    bam = write_bam(
        tmp_path / "reads.bam", [ALLELE],
        make_pair("retained") + [extra] + improper + make_pair("orphan")[:1],
    )
    result = score_per_allele(ALLELE, bam, min_ecnt=1)
    assert result is not None
    assert result["qnames"].to_list() == ["retained"]
    assert result["scores"].item() == pytest.approx(PAIR_SCORE)


def test_secondary_pair_remains_evidence_for_an_alternative_allele(tmp_path):
    """Secondary status does not discard a properly paired alternative mapping."""
    other = "hla_a_02_01_01"
    bam = write_bam(
        tmp_path / "reads.bam", [ALLELE, other],
        make_pair("shared", 0) + make_pair("shared", 1, secondary=True),
    )
    result = score_per_allele(other, bam, min_ecnt=1)
    assert result is not None
    assert result["qnames"].to_list() == ["shared"]
    assert result["allele"].to_list() == [other]
    assert result["scores"].item() == pytest.approx(PAIR_SCORE)


def test_reverse_mate_uses_qualities_in_stored_alignment_order(tmp_path):
    """A first-position mismatch on read2 uses its first stored quality."""
    pair = make_pair("reverse")
    pair[1].query_qualities = [10, 20, 30, 40]
    pair[1].set_tag("MD", "0C3")
    pair[1].set_tag("NM", 1)
    bam = write_bam(tmp_path / "reads.bam", [ALLELE], pair)
    result = score_per_allele(ALLELE, bam, min_ecnt=1)
    assert result is not None
    expected = (
        8 * 23 + 4 * math.log(0.999) + math.log(0.1 / 3)
        + math.log(0.99) + math.log(0.999) + math.log(0.9999)
    )
    assert result["scores"].item() == pytest.approx(expected, abs=1e-10, rel=0)


def test_trailing_soft_clip_does_not_contribute_to_pair_score(tmp_path):
    """Only the four aligned qualities contribute; trailing low qualities do not."""
    pair = make_pair("trailing_clip")
    pair[0].query_sequence = "AAAAAA"
    pair[0].query_qualities = [10, 20, 30, 40, 1, 2]
    pair[0].cigarstring = "4M2S"
    bam = write_bam(tmp_path / "reads.bam", [ALLELE], pair)
    result = score_per_allele(ALLELE, bam, min_ecnt=1)
    assert result is not None
    expected = (
        8 * 23 + math.log(0.9) + math.log(0.99) + math.log(0.999)
        + math.log(0.9999) + 4 * math.log(0.999)
    )
    assert result["scores"].item() == pytest.approx(expected, abs=1e-10, rel=0)
