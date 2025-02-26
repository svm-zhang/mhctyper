import pytest

from mhctyper.hla_allele import (
    HLAllele,
    HLAllelePattern,
    decompose,
    reduce_resolution,
)


@pytest.fixture(autouse=True, scope="module")
def hla_allele():
    return HLAllele


@pytest.fixture(autouse=True, scope="module")
def hla_allele_default_pattern():
    return HLAllelePattern()


def test_hlallele_default_pattern():
    ap = HLAllelePattern()
    assert ap.prefix == "(?P<prefix>(:?HLA|hla)[-_])"
    assert ap.locus == "(?P<locus>[a-zA-Z]+[0-9]?)"
    assert (
        ap.digit_fields
        == "(?P<digit_fields>(?:[0-9]+[:_]){0,3}[0-9NLSCAQnlscaq]+)"
    )


@pytest.mark.parametrize("resolution, expect", [(-1, 1), (5, 4)])
def test_hlaallele_resolution_value(resolution, expect):
    """
    Test if HLAllelePattern.resolution is set correctly,
    when its value is either less than 1 or greater than 4
    """
    ap = HLAllelePattern(resolution=resolution)
    assert ap.resolution == expect


def test_hlaallele_digit_field_with_custom_resolution():
    ap = HLAllelePattern(resolution=2)
    assert (
        ap.digit_fields
        == "(?P<digit_fields>(?:[0-9]+[:_]){0,1}[0-9NLSCAQnlscaq]+)"
    )
    ap = HLAllelePattern(resolution=0)
    assert (
        ap.digit_fields
        == "(?P<digit_fields>(?:[0-9]+[:_]){0,0}[0-9NLSCAQnlscaq]+)"
    )
    ap = HLAllelePattern(resolution=7)
    assert (
        ap.digit_fields
        == "(?P<digit_fields>(?:[0-9]+[:_]){0,3}[0-9NLSCAQnlscaq]+)"
    )


def test_empty_pattern():
    with pytest.raises(ValueError):
        HLAllelePattern(locus="", prefix="   ")


def test_pattern_with_empty_whitespaces():
    """
    Test if whitespaces are trimmed off from ends when there are
    whitespaces in the given pattern
    """
    ap = HLAllelePattern(locus="    ([A-Z]+[0-9]) ")
    assert ap.locus == "(?P<locus>([A-Z]+[0-9]))"
    ap = HLAllelePattern(locus="\t([A-Z]+[0-9])\n")
    assert ap.locus == "(?P<locus>([A-Z]+[0-9]))"


@pytest.mark.parametrize(
    "allele, nth_field, expect",
    [
        ("A*02:06", 1, "A*02"),
        ("hla_a_01_01_01", 2, "hla_a_01_01"),
        ("hla_drb1*11_01_01:347N", 3, "hla_drb1*11_01_01"),
        ("HLA-dqb1*06:02", 1, "HLA-dqb1*06"),
        ("hla_drb1*11_01_01:347N", -1, "hla_drb1*11"),
    ],
)
def test_reduce_resolution(allele, nth_field, expect, recwarn):
    ap = HLAllelePattern(resolution=nth_field)
    reduced_allele = reduce_resolution(allele, ap)
    assert len(recwarn) == 0
    assert reduced_allele == expect


@pytest.mark.parametrize(
    "allele, nth_field",
    [
        ("A*02", 1),
        ("hla_a_01_01_01", 3),
        ("hla_drb1*11_01_01:347N", 4),
        ("HLA-dqb1*06:02", 3),
    ],
)
def test_reduce_resolution_warning(allele, nth_field):
    with pytest.warns(
        RuntimeWarning, match="No resolution reduced for input allele="
    ):
        ap = HLAllelePattern(resolution=nth_field)
        reduce_resolution(allele, ap)


@pytest.mark.parametrize(
    "allele, nth_field, expect",
    [
        ("hla_a_01-01-01-01", 2, "hla_a_01-01"),
        ("DQB1*11-01-01-01", 2, "DQB1*11-01"),
        ("HLA-B*11:01-01_01", 3, "HLA-B*11:01-01"),
    ],
)
def test_reduce_resolution_with_custom_pattern(allele, nth_field, expect):
    """
    Test reduce resolution with custom pattern
    """
    # Allow digit field using either -, :, _, or a mixture of 3 as separator
    ap = HLAllelePattern(digit_field_sep=["-", ":", "_"], resolution=nth_field)
    assert reduce_resolution(allele, ap) == expect


def test_decompose_empty_allele(hla_allele_default_pattern):
    with pytest.raises(ValueError):
        decompose("", hla_allele_default_pattern)


@pytest.mark.parametrize(
    "allele",
    [
        "jsldkfjlasdjwe",
        "A*",
        "HLA",
        "HLA-C*",
        "hla_drb1_",
        "hla_dqb1_02_01_01_",
        "hla_c*01-01-01",
        "hla_c-01-01-01",
        "_01_01_01",
        "*11:01:01",
        "01-01-01",
        "11:01:01",
        "11_01_01",
    ],
)
def test_decompose_fail_to_find_given_pattern(
    allele, hla_allele_default_pattern
):
    with pytest.raises(ValueError):
        reduce_resolution(allele, hla_allele_default_pattern)
        decompose(allele, hla_allele_default_pattern)


@pytest.mark.parametrize(
    "allele, expect",
    [
        ("A*01:01:01", "A"),
        ("HLA-c*01:03", "c"),
        ("HLA-B*04:16N", "B"),
        ("HLA-A*02:376N", "A"),
        ("hla_dqb1_01_01_01", "dqb1"),
        ("hla_drb1_11_01_01_01", "drb1"),
    ],
)
def test_decompose_locus(allele, expect, hla_allele_default_pattern):
    assert decompose(allele, hla_allele_default_pattern).locus == expect


@pytest.mark.parametrize(
    "allele",
    [
        ("QQQ*01:01:01"),
        ("HLA-ivy*01:03"),
        ("hla_drq1_01_01_01"),
    ],
)
def test_decompose_bad_locus(allele, hla_allele_default_pattern):
    with pytest.raises(ValueError):
        decompose(allele, hla_allele_default_pattern)


@pytest.mark.parametrize(
    "allele, expect",
    [
        ("A*01:01:01", "01:01:01"),
        ("HLA-B*04:16N", "04:16N"),
        ("hla-c*02:376n", "02:376n"),
        ("hla_c_07", "07"),
    ],
)
def test_decompose_digit_field(allele, expect, hla_allele_default_pattern):
    assert decompose(allele, hla_allele_default_pattern).digit_field == expect


def test_hlallele_str():
    """Test str(HLAllele())"""
    assert (
        str(
            HLAllele(prefix="", locus="C", digit_field="01:01:01:02n", sep="*")
        )
        == "C*01:01:01:02n"
    )
    assert (
        str(
            HLAllele(
                prefix="hla_",
                locus="DRB1",
                digit_field="01_01_01:02n",
                sep="_",
            )
        )
        == "hla_DRB1_01_01_01:02n"
    )
