import math
import sys
from functools import partial
from multiprocessing import get_context

import numpy as np
import polars as pl
from tinyscibio import BAMetadata, walk_bam
from tqdm import tqdm

from .hla_allele import HLAllelePattern, decompose
from .logger import logger


def score_log_liklihood(
    row: dict[str, np.ndarray], scale=math.exp(23)
) -> float:
    score: float = 0.0
    start, end = 0, 0
    base_qs = row["bqs"]
    md = row["mds"]
    for i in range(len(md)):
        if md[i].startswith("^"):
            continue
        block = np.array([])
        if md[i].isalpha():
            # mismatches
            end = start + len(md[i])
            block = np.array([10 ** (-k / 10) / 3 for k in base_qs[start:end]])
        else:
            # matches
            end = start + int(md[i])
            block = np.array([1 - 10 ** (-k / 10) for k in base_qs[start:end]])
        score += np.sum(np.log(np.multiply(block, scale)))
        start = end
    return score


def score_per_allele(
    allele: str, bam_fspath: str, min_ecnt: int
) -> pl.DataFrame | None:
    ap = HLAllelePattern()
    hla_allele = decompose(allele, ap)
    hla_gene = f"{hla_allele.prefix}{hla_allele.locus}"

    bametadata = BAMetadata(bam_fspath)
    df = walk_bam(
        bam_fspath,
        allele,
        exclude=3584,
        return_ecnt=True,
        return_bq=True,
        return_md=True,
        return_qname=True,
    )
    # walk_bam returns rname as id
    # convert it back to rname
    df = df.with_columns(
        pl.col("rnames").replace_strict(bametadata.idx2seqname()),
    )
    df = df.filter(
        (pl.col("indel_ecnt") == 0)
        & (pl.col("mm_ecnt") <= min_ecnt)  # good aln
        & (pl.col("propers"))
    )
    # we only keep read in pairs after applying above filters
    # do not do the following before the above
    df = (
        df.with_columns(
            pl.col("qnames").count().over("qnames").alias("n"),
        )
        .filter(pl.col("n") == 2)
        .drop("n")
    )

    if df.shape[0] == 0:
        return None
    # we convert dtype of bqs and mds column to get ready for score likelihood
    df = df.with_columns(
        pl.col("bqs").map_elements(
            lambda x: x.tolist(),
            return_dtype=pl.List(pl.UInt8),
        ),
        pl.col("mds").map_elements(
            lambda x: x, return_dtype=pl.List(pl.String)
        ),
    )
    # score the likelihood given bqs and mds
    df = df.with_columns(
        pl.struct(["bqs", "mds"])
        .map_elements(score_log_liklihood, return_dtype=pl.Float64)
        .alias("scores"),
    )
    # Sum up the score per aligned pair
    df = df.group_by("qnames").agg(pl.col("scores").sum())
    df = df.with_columns(allele=pl.lit(allele), gene=pl.lit(hla_gene))
    return df


def score_a_one(
    alleles_to_score: list[str],
    bam: str,
    min_ecnt: int,
    nproc: int = 8,
) -> pl.DataFrame:
    logger.info("Score first allele.")
    try:
        score_tables: list[pl.DataFrame] = []
        logger.debug(f"# alleles to score: {len(alleles_to_score)}.")
        with get_context("spawn").Pool(processes=nproc) as pool:
            allele_iterator = tqdm(
                pool.imap_unordered(
                    partial(
                        score_per_allele, bam_fspath=bam, min_ecnt=min_ecnt
                    ),
                    alleles_to_score,
                ),
                total=len(alleles_to_score),  # shows pbar
                desc="Score first allele: ",
                ncols=100,  # define width
            )
            for res in allele_iterator:
                if res is None:
                    continue
                # display allele being processed
                allele_iterator.set_postfix(
                    {"allele": f"{res['allele'].unique().item()}"}
                )
                score_tables.append(res)
        if not score_tables:
            raise ValueError("Failed to score for any first alleles.")
        scores = pl.concat([s for s in score_tables])
        logger.info(f"Alleles scored: {len(score_tables)}.")
    except ValueError as e:
        logger.error(e)
        sys.exit(1)
    return scores


def score_second_by_gene(
    gene: str, a1_scores: pl.DataFrame, a1_winners: pl.DataFrame
) -> pl.DataFrame:
    a1_winners = a1_winners.filter(pl.col("gene") == gene)
    a1_scores = a1_scores.filter(pl.col("gene") == gene)
    score_table = a1_scores.join(a1_winners, on=["qnames", "gene"], how="left")
    score_table = score_table.with_columns(
        pl.col("scores_right").fill_null(0.0)
    )
    score_table = score_table.with_columns(
        factor=pl.col("scores") / (pl.col("scores") + pl.col("scores_right"))
    )
    score_table = score_table.with_columns(
        scores=pl.col("scores") * pl.col("factor")
    )
    return score_table


def score_a_two(
    a1_scores: pl.DataFrame,
    a1_winners: pl.DataFrame,
    nproc: int = 8,
) -> pl.DataFrame:
    logger.info("Score second allele.")
    try:
        score_tables: list[pl.DataFrame] = []
        genes = a1_winners["gene"].unique().to_list()
        nproc = min(nproc, len(genes))
        with get_context("spawn").Pool(processes=nproc) as pool:
            gene_iterator = tqdm(
                pool.imap_unordered(
                    partial(
                        score_second_by_gene,
                        a1_scores=a1_scores,
                        a1_winners=a1_winners,
                    ),
                    genes,
                ),
                total=len(genes),
                desc="Score second allele: ",
                ncols=100,
            )
            for res in gene_iterator:
                if res is None:
                    continue
                score_tables.append(res)
        if not score_tables:
            raise ValueError("Failed to score for any second alleles.")
        a2_scores = pl.concat(score_tables)
        return a2_scores
    except ValueError as e:
        logger.error(e)
        sys.exit(1)


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
