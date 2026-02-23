# mhctyper

[![PyPI version](https://img.shields.io/pypi/v/mhctyper)](https://pypi.org/project/mhctyper/)
![Python versions](https://img.shields.io/pypi/pyversions/mhctyper)
![License](https://img.shields.io/pypi/l/mhctyper)

Faster MHC class I and II typing based on Polysolver algorithm.

## Features

- Supports both class I and II typing with high
  [accuracy](https://github.com/svm-zhang/hla_benchmark?tab=readme-ov-file).
- Achieves dramatic speedups by eliminating thousands of disk-heavy I/O and
  leveraging Polars for lighting-fast data manipulation
- Features a robust CLI and standardized packaging, ensuring seamless integration
  existing workflow/pipeline.

## Highlights


### Alignment tuning

Not all the alignments are suitable for typing. `mhctyper` provides a single,
customizable parameter, `--min_ecnt`, to filter
the alignments used in the typing process. Lowering this value prioritizes
"high quality" data. On the 1000 genome dataset, a value of `1` yields optimal results.
By default, all alignments are included regardless of the number of mismatch event counts.


### Unified output

`mhctyper` replaces the "thousands of files" approach with a single, structured
tabular output. This unified format eliminates directory clutter and allows for
streamlined searching, querying, and downstream analysis.

### Smart skip
Scoring for the first allele accounts for the majority `mhctyper`'s runtime.
To optimize efficiency, the first-allele score table is cached once
generated; `mhctyper` then automatically uses the cached results for subsequent
steps to avoid redundant BAM traversal.


## Installation

You can install mhctyper from PyPI

=== "pip"

    ```bash
    pip install mhctyper
    ```

## Quick start

`mhctyper` simply requires 2 inputs:

- Alignment to HLA alleles in BAM format: `$bam`.
- Population frequency from the original `polysolver`: `HLA_FREQ.txt`.

```bash
mhctyper --bam "$bam" \
--freq "HLA_FREQ.txt" \
--outdir "$outdir" \
--nproc 8
```

## Output explain

The above `mhctyper` command yields 3 output files:

- `{RG_SM}.a1.tsv`: score table for the first allele.
- `{RG_SM}.a2.tsv`: score table for the second allele.
- `{RG_SM}.hlatyping.res.tsv`: HLA typing result.

`{RG_SM}` represents the value of `SM` field of read group provided in the
given BAM file. `mhctyper` checks the existence of read group information
and terminates if either no or more than one read group value set.

### Score table

Score table has following format with self-explanatory columns:

| qnames               | scores             | allele            | gene   |
| :------------------- | :----------------- | :---------------- | :----- |
| SRR702076.1195861    | 4645.90100859441   | hla_a_01_01_01_01 | hla_a  |
| SRR702076.11670250   | 4644.914284618211  | hla_a_01_01_01_01 | hla_a  |
| SRR702076.3308570    | 4645.911824273519  | hla_a_01_01_01_01 | hla_a  |
| SRR702076.23151566   | 4645.900731687712  | hla_a_01_01_01_01 | hla_a  |

Each row represents the score typed for an allele from one read (not a pair).

### HLA typing result

HLA typing result contains 4 columns:

- allele: typed alleles.
- gene: HLA gene locus.
- tot_scores: total loglikelihood scores for 2 alleles per gene locus.

| allele           | gene     | tot_scores   | sample   |
| :--------------- | :------- | :----------- | :------- |
| hla_a_11_01_01   | hla_a    | 4120184.7405 | NA18740  |
| hla_a_11_01_01   | hla_a    | 2060092.3702 | NA18740  |
| hla_b_13_01_01   | hla_b    | 1054296.1982 | NA18740  |
| hla_b_40_01_02   | hla_b    | 1221557.8978 | NA18740  |
| hla_c_03_04_04   | hla_c    | 1474826.4096 | NA18740  |
| hla_c_07_02_10   | hla_c    | 1913741.4495 | NA18740  |

## Filters applied silently

`mhctyper` quietly applies the following filters during typing:

- Only properly aligned read pairs are used.
- QC-failed, supplementary, and duplicate-marked alignments are removed
(exclude sam flag = 3584).
- Alignments with indels are removed.
- **Alignments with mismatches more than specified value of `--min_ecnt` are
removed (see below)**.
- HLA alleles (their 4-digits representation) who have sum of zero population allele
frequencies across all populations defined in the `HLA_FREQ.txt` file are
excluded from typing.


## Polysolver comparison

While `mhctyper` implements the core Polysolver algorithm, results may not be
identical in every case. However, a high degree of concordance between the two
tools should be expected in the majority runs. 


### Troubleshooting

If `mhctyper` and Polysovler produce different results, follow these
steps to pinpoint the discrepancy:

1. **Reads concordance**: Verify if the "fished" read IDs are identical;
2. **Bam statistics**: Compare alignment stats of the realigned BAM files
   using the `sametoosl flagstat`;
3. **Score difference**: Check scores for the specific alleles in question.
   Is the numerical difference marginal or significant?
4. **Manual inspection**: Generate BAM files for the alleles in question to
   compare alignment stats, CIGAR and MD strings, and/or view alignments in IGV.

## Notes

### Race

Unlike the original `polysovler` algorithm, `mhctyper` does not provide a
`race` option. This is intentional because most of the case we do not have
such demographic information to begin with. Upon testing on 1000 Genome
dataset, this does not affect the typing result.


### Debug mode

If encountering errors when running `mhctyper`,
the debug mode can be toggled to generate a log file under the output folder
you specify. `mhctyper` will automatically ignore the `--nproc` value and
switch to the single process mode. Please share this debug file when opening an issue.

```bash
mhctyper --bam "$bam" \
    --freq "HLA_FREQ.txt" \
    --outdir "$outdir" \
    --debug
```

### Overwrite

Use the `--overwrite` flag to force a full clean `mhctyper` re-run. Cached
results from previous run will be deleted automatically.


## Citation

Please cite the original
[Polysolver](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4747795/) paper.

If you use `mhctyper`, please cite the github repo.
