#!/usr/bin/env python

from __future__ import annotations

import polars as pl
from tinyscibio import BAMetadata, get_parent_dir, make_dir

from .cli import parse_cmd
from .logger import logger
from .score_alleles import get_winners, score_a_one, score_a_two
from .utils import (
    collect_alleles_to_type,
    load_allele_pop_freq,
    load_rg_sm_from_bam,
)


def run_mhctyper() -> int:
    parser = parse_cmd()
    args = parser.parse_args()

    make_dir(args.outdir, exist_ok=True, parents=True)

    # set up logger accordingly
    if args.debug:
        debug_log_fspath = args.outdir / f"{__name__}.debug.log"
        if debug_log_fspath.exists():
            debug_log_fspath.unlink()
        logger.initialize(args.debug, debug_log_fspath)
    else:
        logger.initialize(args.debug)

    logger.info(f"Start HLA typing from given BAM file: {args.bam}")

    allele_pop_freq = load_allele_pop_freq(freq_fspath=args.freq)

    bam_metadata = BAMetadata(args.bam)
    # collect all alleles to type from BAM header
    alleles_to_type = collect_alleles_to_type(
        bam_metadata, kept=allele_pop_freq["Allele"].to_list()
    )

    rg_sm = load_rg_sm_from_bam(bam_metadata)

    out_a1 = args.outdir / f"{rg_sm}.a1.tsv"
    out_a2 = args.outdir / f"{rg_sm}.a2.tsv"
    hla_res = args.outdir / f"{rg_sm}.hlatyping.res.tsv"
    if args.overwrite:
        logger.info("Overwrite specified. Delete results previously computed.")
        out_a1.unlink(missing_ok=True)
        out_a2.unlink(missing_ok=True)
        hla_res.unlink(missing_ok=True)

    a1_scores = pl.DataFrame()
    if not out_a1.exists():
        a1_scores = score_a_one(
            alleles_to_score=alleles_to_type,
            bam=args.bam,
            min_ecnt=args.min_ecnt,
            nproc=args.nproc,
            debug=args.debug,
        )
        a1_scores.write_csv(out_a1, separator="\t")
    else:
        logger.info("Found scores of first alleles previously computed.")
        a1_scores = pl.read_csv(out_a1, separator="\t")

    logger.info("Get winner for the first typed allele.")
    a1_winners = get_winners(allele_scores=a1_scores)
    winner_scores = a1_scores.join(
        a1_winners, on=["gene", "allele"], how="inner"
    )

    a2_scores = score_a_two(
        a1_scores=a1_scores,
        a1_winners=winner_scores,
        nproc=args.nproc,
    )
    a2_scores.write_csv(out_a2, separator="\t")
    logger.info("Get winner for the second typed allele.")
    a2_winners = get_winners(allele_scores=a2_scores)

    logger.info("Combine winnes for both first and second alleles.")
    hla_res_df = pl.concat([a1_winners, a2_winners])
    hla_res_df = hla_res_df.with_columns(sample=pl.lit(rg_sm)).sort(
        by="allele"
    )
    logger.info(f"Final HLA typing result: {hla_res_df}")
    hla_res_df.write_csv(hla_res, separator="\t")
    return 0
