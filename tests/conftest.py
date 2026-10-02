"""Shared state isolation and one small public-workflow input."""

import pytest

from mhctyper.logger import logger

from .bam_helpers import make_pair, write_bam


@pytest.fixture(autouse=True)
def isolated_logger():
    """Restore shared logging state so in-process tests do not affect each other."""
    handlers, level, configured = logger.handlers[:], logger.level, logger._configured
    logger.handlers = []
    logger._configured = False
    try:
        yield
    finally:
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            handler.close()
        logger.handlers = handlers
        logger.setLevel(level)
        logger._configured = configured


@pytest.fixture
def typing_inputs(tmp_path):
    """Shared/private pairs yield A heterozygosity and a DRB1 homozygous call."""
    references = [
        "hla_a_01_01_01", "hla_a_02_01_01", "hla_a_03_01_01",
        "hla_a_04_01_01", "hla_drb1_01_01_01",
    ]
    records = (
        make_pair("shared", 0)
        + make_pair("a1_private", 0)
        + make_pair("shared", 1, secondary=True)
        + make_pair("a2_private", 1)
        + make_pair("drb1_private", 4)
    )
    # A3 would win on volume at the default threshold, but each pair has a
    # mate with two mismatches. min_ecnt=1 must eliminate all three pairs.
    for i in range(3):
        pair = make_pair(f"above_threshold_{i}", 2)
        pair[0].set_tag("MD", "0C0C2")
        pair[0].set_tag("NM", 2)
        records += pair
    # A4 has abundant evidence but zero frequency, so it must never compete.
    for i in range(3):
        records += make_pair(f"excluded_{i}", 3)
    bam = write_bam(tmp_path / "input.bam", references, records)
    freq = tmp_path / "freq.tsv"
    freq.write_text(
        "Allele\tlocal_population\n"
        "hla_a_01_01\t0.1\n"
        "hla_a_02_01\t0.1\n"
        "hla_a_03_01\t0.1\n"
        "hla_a_04_01\t0\n"
        "hla_drb1_01_01\t0.2\n"
    )
    return bam, freq
