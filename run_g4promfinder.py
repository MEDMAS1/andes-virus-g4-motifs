#!/usr/bin/env python3
"""Run G4PromFinder on FASTA files in data/ and save results under results/G4PromFinder."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "G4PromFinder"
TOOL_DIR = ROOT / "library" / "G4PromFinder"
SCRIPT = TOOL_DIR / "G4PromFinder.py"


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def run_command(command, cwd, input_text):
    return subprocess.run(
        command,
        cwd=str(cwd),
        check=True,
        text=True,
        input=input_text,
        capture_output=True,
    )


def main():
    parser = argparse.ArgumentParser(description="Run G4PromFinder for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        for filename in ("Promoter coordinates.txt", "ABOUT PROMOTERS.txt"):
            path = out_dir / filename
            if path.exists():
                path.unlink()

        command = [sys.executable, str(SCRIPT)]
        result = run_command(command, out_dir, f"{str(fasta)}\n")
        if result.stderr:
            print(result.stderr, end="")
        print(result.stdout, end="")

        coords = out_dir / "Promoter coordinates.txt"
        about = out_dir / "ABOUT PROMOTERS.txt"

        coords_written = coords.exists()
        about_written = about.exists()

        if coords_written:
            coords.rename(out_dir / f"{fasta.stem}_g4promfinder_coordinates.txt")
        if about_written:
            about.rename(out_dir / f"{fasta.stem}_g4promfinder_about.txt")

        if coords_written or about_written:
            print(f"Wrote results to {out_dir}")
        else:
            raise SystemExit(f"G4PromFinder did not generate output for {fasta}")


if __name__ == "__main__":
    main()
