import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"


def run_cli(args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(
        [sys.executable, "-m", "funtaxfinder", *args],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def test_cli_help():
    result = run_cli(["--help"])
    assert result.returncode == 0
    assert "Screen functional microbial ASVs" in result.stdout


def test_cli_basic_screening(tmp_path):
    outdir = tmp_path / "out"
    result = run_cli(
        [
            "--asv-table",
            str(DATA / "asv_table.tsv"),
            "--ko",
            "K02586",
            "--ko-predicted",
            str(DATA / "combined_KO_predicted.tsv"),
            "--outdir",
            str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    subtable = outdir / "sub_asv_table.tsv"
    assert subtable.exists()
    text = subtable.read_text()
    assert "ASV1\t" in text
    assert "ASV2\t" not in text


def test_cli_taxonomy_split(tmp_path):
    outdir = tmp_path / "out"
    result = run_cli(
        [
            "--asv-table",
            str(DATA / "asv_table.tsv"),
            "--ko",
            "K02586",
            "K02588",
            "--ko-predicted",
            str(DATA / "combined_KO_predicted.tsv"),
            "--tax-table",
            str(DATA / "taxonomy.tsv"),
            "--split-rank",
            "genus",
            "--outdir",
            str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (outdir / "matched_asv_taxonomy.tsv").exists()
    assert (outdir / "tax_split" / "tax_split_summary.tsv").exists()
