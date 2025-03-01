# mhctyper

[![PyPI version](https://img.shields.io/pypi/v/mhctyper)](https://pypi.org/project/mhctyper/)
![Python versions](https://img.shields.io/pypi/pyversions/mhctyper)
[![PyPI Downloads](https://img.shields.io/pypi/dm/mhctyper)](https://pypistats.org/packages/mhctyper)
![License](https://img.shields.io/pypi/l/mhctyper)

Polars-accelerated MHC class I and II typing based on Polysolver algorithm.

## Features

- Supports both class I and II typing with good
  [accuracy](https://github.com/svm-zhang/hla_benchmark?tab=readme-ov-file)
- Runtime speedup boosted by polars
- Minimum I/O operations
- Easy integration to workflow/pipeline with better CLI and proper packaging.

## Installation

mhctyper can be installed from PyPI:

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

Please refer to [documentation](https://svm-zhang.github.io/mhctyper) for more details.
