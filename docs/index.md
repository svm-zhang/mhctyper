# mhctyper

[![PyPI version](https://img.shields.io/pypi/v/mhctyper)](https://pypi.org/project/mhctyper/)
![Python versions](https://img.shields.io/pypi/pyversions/mhctyper)

Polars-accelerated MHC class I and II typing based on Polysolver algorithm.

## Features

- Supports both class I and II typing with good
  [accuracy](https://github.com/svm-zhang/hla_benchmark?tab=readme-ov-file).
- Runtime speedup boosted by polars.
- Minimum I/O operations.
- Easy integration to workflow/pipeline with better CLI and proper packaging.

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

```text
qnames  scores  allele  gene
SRR702076.1195861       4645.90100859441        hla_a_01_01_01_01       hla_a
SRR702076.11670250      4644.914284618211       hla_a_01_01_01_01       hla_a
SRR702076.3308570       4645.911824273519       hla_a_01_01_01_01       hla_a
SRR702076.23151566      4645.900731687712       hla_a_01_01_01_01       hla_a
```

Each row represents the score typed for an allele from one read (not a pair).

### HLA typing result

HLA typing result contains 4 columns:

- allele: typed alleles.
- gene: HLA gene locus.
- tot_scores: total loglikelihood scores for 2 alleles per gene locus.

```text
   allele  gene    tot_scores  sample
   hla_a_11_01_01  hla_a   4120184.7405    NA18740
   hla_a_11_01_01  hla_a   2060092.3702    NA18740
   hla_b_13_01_01  hla_b   1054296.1982    NA18740
   hla_b_40_01_02  hla_b   1221557.8978    NA18740
   hla_c_03_04_04  hla_c   1474826.4096    NA18740
   hla_c_07_02_10  hla_c   1913741.4495    NA18740
```

## Silently applied filters

`mhctyper` quietly applies the following filters during typing:

- Only properly aligned read pairs are used.
- QC-failed, supplementary, and duplicate-marked alignments are removed
(exclude sam flag = 3584).
- Alignments with indels are removed.
- Alignments with mismatches more than specified value of `--min_ecnt` are
removed.
- HLA alleles (their 4-digits representation) who have zero population allele
frequency across all populations defined in the `HLA_FREQ.txt` file are
excluded from typing.

## Notes

### Optimal value for --min_ecnt

`mhctyper` provides a single customizable parameter `--min_ecnt` to control for
the alignments used in the typing. The lower the value, the more "high quality"
data is used (but less in quantity). On the 1000 genome dataset, I found
setting `--min_ecnt` to 1 gives fairly good result. You can tune this parameter
to find the best value for your data. By default, all alignments, regardless
of the number of mismatch events allowed, will be used.

### Race

Unlike the original `polysovler` algorithm, `mhctyper` does not provide a
`race` option. This is intentional because most of the case we do not have
such demographic information to begin with. From testing on 1000 genome
dataset, this does not affect the typing result.

### Break and continue

Scoring for the first allele will take the majority runtime of `mhctyper`.
Therefore, once the score table for the first allele is generated, `mhctyper`
will not repeat it to avoid walking through the BAM file for another time.
You can use `--overwrite` option to force `mhctyper` to re-score the first
allele.

### Run mhctyper in debug mode

If you encounter errors when running `mhctyper` on your data, you can toggle
the debug mode to generate a log file under the output folder you specify.
Sharing the debug file will very useful for quickly identifying the problem(s).
In debug mode, `mhctyper` will automatically switch to single process and
ignore the value (8 in the example below) provided via `--nproc` option.

```bash
mhctyper --bam "$bam" \
    --freq "HLA_FREQ.txt" \
    --outdir "$outdir" \
    --nproc 8 \
    --debug
```

## Citation

Please cite the original
[Polysolver](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4747795/) paper.

If you use `mhctyper`, please cite the github repo.
