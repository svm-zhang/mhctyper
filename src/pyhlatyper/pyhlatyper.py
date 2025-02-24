#!/usr/bin/env python

from __future__ import annotations

import warnings
from pprint import pprint
from typing import TYPE_CHECKING

import polars as pl
from tinyscibio import BAMetadata, get_parent_dir, make_dir

from .cli import parse_cmd
from .hla_allele import HLAllelePattern, reduce_resolution
from .score_alleles import score_a_one, score_a_two, score_per_allele

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Optional

    from tinyscibio import _PathLike


def get_winners(allele_scores: pl.DataFrame) -> pl.DataFrame:
    # round scores to 4 decimal places to avoid precision
    # problem when getting alleles whose scores equal to max scores
    tot_scores = allele_scores.group_by(["allele", "gene"]).agg(
        pl.col("scores").sum().round(4)
    )
    winners = tot_scores.filter(
        pl.col("scores") == pl.col("scores").max().over("gene")
    )
    # when there is a tie in scores for each gene group,
    # select allele with least string value lexicographically
    # e.g. hla_a_26_01_24 and hla_a_26_01_01 (latter selected)
    winners = winners.sort(by=["allele"]).unique(subset="gene", keep="first")
    winners = winners.rename({"scores": "tot_scores"})
    return winners


def load_allele_pop_freq(freq_fspath: _PathLike) -> pl.DataFrame:
    """Load allele population frequency data"""
    freq_df = pl.read_csv(freq_fspath, separator="\t")
    # remove supertype has all zero pop frequency
    return freq_df.filter(
        pl.fold(0, lambda acc, s: acc + s, pl.all().exclude(pl.String)) > 0.0
    )


def collect_alleles_to_type(
    bam_metadata: BAMetadata,
    kept: Optional[Sequence[str]] = None,
    # bam: _PathLike, kept: Optional[Sequence[str]] = None
) -> list[str]:
    """Collect HLA alleles to type"""
    ap = HLAllelePattern(resolution=2)
    alleles_df = pl.DataFrame({"Allele": bam_metadata.seqnames()})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # alleles = list(map(lambda a: reduce_resolution(a, ap), alleles))
        # if kept is not None:
        #     alleles = [a for a in alleles if a in kept]
        if kept is not None:
            alleles_df = alleles_df.with_columns(
                pl.col("Allele")
                .map_elements(
                    lambda c: reduce_resolution(c, ap), return_dtype=pl.String
                )
                .alias("D4")
            ).filter(pl.col("D4").is_in(kept))
    alleles = alleles_df["Allele"].to_list()
    if not alleles:
        raise ValueError("Failed to collect any alleles for typing")
    pprint(f"[INFO]: Collected {len(alleles)} alleles to type")

    return alleles


def run_pyhlatyper() -> int:
    parser = parse_cmd()
    args = parser.parse_args()

    outdir = get_parent_dir(args.out)
    make_dir(outdir, exist_ok=True)

    allele_pop_freq = load_allele_pop_freq(freq_fspath=args.freq)

    bam_metadata = BAMetadata(args.bam)
    # collect all alleles to type from BAM header
    alleles_to_type = collect_alleles_to_type(
        bam_metadata, kept=allele_pop_freq["Allele"].to_list()
    )

    # score_df = score_per_allele_v2(
    #     allele="hla_c_01_02_01", bam_fspath=args.bam, min_ecnt=1
    # )
    # with pl.Config(fmt_str_lengths=1000, set_tbl_cols=-1):
    #     print(score_df.head())
    # raise SystemExit

    rg = bam_metadata.read_groups
    if len(rg) > 1:
        raise ValueError(f"Found more than 1 read groups: {rg}")
    rg_sm = next(iter(rg)).get("SM", None)
    if rg_sm is None:
        raise ValueError(
            "Failed to get SM from read group. "
            f"Please check parsed read group: {rg}"
        )
    out_a1 = outdir / f"{rg_sm}.a1.tsv"
    a1_scores = pl.DataFrame()
    if not out_a1.exists():
        a1_scores = score_a_one(
            alleles_to_score=alleles_to_type,
            bam=args.bam,
            min_ecnt=args.min_ecnt,
            nproc=args.nproc,
        )
        a1_scores.write_csv(out_a1, separator="\t")
    else:
        a1_scores = pl.read_csv(out_a1, separator="\t")

    a1_winners = get_winners(allele_scores=a1_scores)
    winner_scores = a1_scores.join(
        a1_winners, on=["gene", "allele"], how="inner"
    )
    with pl.Config(fmt_str_lengths=1000, set_tbl_cols=-1):
        print(a1_scores.head())
        print(a1_winners)
    raise SystemExit

    out_a2 = outdir / f"{rg_sm}.a2.tsv"
    a2_scores = score_a_two(
        a1_scores=a1_scores,
        a1_winners=winner_scores,
        out=out_a2,
        nproc=args.nproc,
    )
    a2_winners = get_winners(allele_scores=a2_scores)

    hla_res = f"{outdir}/{rg_sm}.hlatyping.res.tsv"
    hla_res_df = pl.concat([a1_winners, a2_winners])
    hla_res_df = hla_res_df.with_columns(sample=pl.lit(rg_sm)).sort(
        by="allele"
    )
    print(hla_res_df)
    hla_res_df.write_csv(hla_res, separator="\t")
    return 0
