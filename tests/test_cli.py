"""Smoke-test the installed command path without duplicating scoring cases."""

import subprocess
import sys

import polars as pl


def test_cli_types_a_bam_into_the_requested_directory(typing_inputs, tmp_path):
    """Required paths and scoring options reach the real typing workflow."""
    bam, freq = typing_inputs
    outdir = tmp_path / "cli_results"
    completed = subprocess.run(
        [
            sys.executable, "-m", "mhctyper", "--bam", str(bam),
            "--freq", str(freq), "--outdir", str(outdir),
            "--min_ecnt", "1", "--nproc", "1",
        ],
        capture_output=True, text=True, timeout=60,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    calls = pl.read_csv(outdir / "sample_1.hlatyping.res.tsv", separator="\t")
    assert calls["allele"].to_list() == [
        "hla_a_01_01_01", "hla_a_02_01_01",
        "hla_drb1_01_01_01", "hla_drb1_01_01_01",
    ]
    assert calls["sample"].to_list() == ["sample_1"] * 4
