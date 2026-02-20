from __future__ import annotations

import sys
import warnings
from collections.abc import Sequence
from typing import TYPE_CHECKING

import polars as pl
import polars.selectors as cs
from tinyscibio import BAMetadata

from .hla_allele import HLAllelePattern, reduce_resolution
from .logger import logger

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Optional

    from tinyscibio import _PathLike


def load_rg_sm_from_bam(bam_metadata: BAMetadata) -> str:
    logger.info("Load SM value from read group in the given BAM.")
    try:
        rg = bam_metadata.read_groups
        logger.debug(f"bam_metadata.read_groups: {rg=}")
        if len(rg) > 1:
            raise ValueError(f"Found more than 1 read groups: {rg}")
        rg_sm = next(iter(rg)).get("SM", "")
        logger.debug(f"{rg_sm=}")
        if not rg_sm:
            raise ValueError(
                "Failed to get SM from read group. "
                f"Please check parsed read group: {rg}"
            )
        return rg_sm
    except ValueError as e:
        logger.error(e)
        sys.exit(1)


def load_allele_pop_freq(freq_fspath: _PathLike) -> pl.DataFrame:
    """Load allele population frequency data"""
    logger.info("Load HLA alleles from population frequency file.")
    freq_df = pl.read_csv(freq_fspath, separator="\t")
    # remove supertype has all zero pop frequency
    # https://github.com/svm-zhang/mhctyper/issues/4
    # pl.fold returns 0 of dtype 0 when either acc initializes as i32 or
    # return_dtype is not set. The supercast rule does not seem to apply
    # https://docs.pola.rs/api/python/stable/reference/dataframe/api/polars.DataFrame.fold.html
    # Below is old solution:
    # return freq_df.filter(
    #     pl.fold(0, lambda acc, s: acc + s, pl.all().exclude(pl.String)) > 0.0
    # )
    # To fix:
    # pl.fold(acc=0.0, lambda acc, s: acc+s, pl.exclude(pl.String))
    # pl.fold(
    #   acc=0,
    #   lambda acc, s: acc+s,
    #   pl.exclude(pl.String),
    #   return_dtype=pl.Float64
    # )
    # Should be better off using pl.sum_horizontal instead
    return freq_df.filter(pl.sum_horizontal(cs.numeric()) > 0.0)


def collect_alleles_to_type(
    bam_metadata: BAMetadata,
    kept: Optional[Sequence[str]] = None,
) -> list[str]:
    """Collect HLA alleles to type"""
    logger.info("Load HLA alleles to type.")
    try:
        ap = HLAllelePattern(resolution=2)
        alleles_df = pl.DataFrame({"Allele": bam_metadata.seqnames()})
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if kept is not None:
                alleles_df = alleles_df.with_columns(
                    pl.col("Allele")
                    .map_elements(
                        lambda c: reduce_resolution(c, ap),
                        return_dtype=pl.String,
                    )
                    .alias("D4")
                ).filter(pl.col("D4").is_in(kept))
        alleles = alleles_df["Allele"].to_list()
        if not alleles:
            raise ValueError("Failed to collect any alleles for typing")
        logger.info(f"Loaded {len(alleles)} alleles to type.")
        return alleles
    except ValueError as e:
        logger.error(e)
        sys.exit(1)
