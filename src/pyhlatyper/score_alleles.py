import math
import os
from collections.abc import Sequence
from functools import partial
from multiprocessing import get_context

import numpy as np
import polars as pl
import pysam
from tqdm import tqdm

from .bam import count_mismatch_events, count_unaligned_events


def score_log_liklihood(
    base_qs: Sequence[int], md: Sequence[str], scale=math.exp(23)
) -> float:
    score: float = 0.0
    start, end = 0, 0
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
) -> pl.DataFrame:
    ids: list[str] = []
    scores: list[float] = []
    with pysam.AlignmentFile(bam_fspath, "rb") as bamf:
        for aln in bamf.fetch(contig=allele):
            # skip alignments marked as QCFAIL and Supplementary
            if aln.is_qcfail or aln.is_supplementary:
                continue
            # skip alignment that are not properly mapped
            if not aln.is_proper_pair:
                continue
            if aln.query_qualities is None:
                raise TypeError(
                    "QUAL should not be NoneType in the alignment record"
                )
            if aln.query_name is None:
                raise TypeError(
                    "QNAME should not be NoneType in the alignment record"
                )
            if aln.cigarstring is None:
                raise TypeError()
            if not aln.has_tag("MD"):
                raise TypeError()

            n_unaligned_events = count_unaligned_events(aln.cigarstring)
            if n_unaligned_events > 0:
                continue
            n_mm_events = count_mismatch_events(str(aln.get_tag("MD")))
            if n_mm_events > min_ecnt:
                continue

            ids += [aln.query_name]
            scores += [
                score_log_liklihood(
                    aln.query_qualities, str(aln.get_tag("MD"))
                )
            ]

        score_df = pl.DataFrame(
            {
                "ids": ids,
                "scores": scores,
            },
        )
        score_df = (
            score_df.with_columns(pl.col("ids").count().over("ids").alias("n"))
            .filter(pl.col("n") == 2)
            .drop("n")
            .group_by("ids")
            .agg(pl.col("scores").sum())
            .with_columns(allele=pl.lit(allele))
        )
        return score_df


def score_a_one(
    alleles_to_score: list[str],
    bam: str,
    min_ecnt: int,
    nproc: int = 8,
) -> pl.DataFrame:
    score_tables: list[pl.DataFrame] = []
    with get_context("spawn").Pool(processes=nproc) as pool:
        with tqdm(
            total=len(alleles_to_score),
            desc="extract_aln",
            ncols=100,
            leave=False,
        ) as pbar:
            for res in tqdm(
                pool.imap_unordered(
                    partial(
                        score_per_allele, bam_fspath=bam, min_ecnt=min_ecnt
                    ),
                    alleles_to_score,
                )
            ):
                pbar.update()
                if res is None:
                    continue
                score_tables.append(res)
    scores = pl.concat([s for s in score_tables])
    return scores


def score_second_by_gene(
    gene: str, a1_scores: pl.DataFrame, a1_winners: pl.DataFrame
) -> pl.DataFrame:
    a1_winners = a1_winners.filter(pl.col("gene") == gene)
    a1_scores = a1_scores.filter(pl.col("gene") == gene)
    score_table = a1_scores.join(a1_winners, on=["ids", "gene"], how="left")
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


def score_second(
    a1_scores: pl.DataFrame, a1_winners: pl.DataFrame, out: str, nproc: int = 8
) -> pl.DataFrame:
    if os.path.exists(out):
        return pl.read_csv(
            out,
            separator="\t",
        )

    score_tables: list[pl.DataFrame] = []
    genes = a1_winners["gene"].unique().to_list()
    # nproc likely more than number of genes, only spawn # procs
    # based on # genes
    nproc = min(nproc, len(genes))
    # I actually dont think need to parallele this
    # but until that day...
    with get_context("spawn").Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(
                score_second_by_gene,
                a1_scores=a1_scores,
                a1_winners=a1_winners,
            ),
            genes,
        ):
            if res is None:
                continue
            score_tables.append(res)
    a2_scores = pl.concat(score_tables)
    a2_scores.write_csv(out, separator="\t")
    return a2_scores
