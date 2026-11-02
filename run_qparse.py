#!/usr/bin/env python3
"""Run QPARSE on FASTA files in data/ and save one result folder per input."""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "QPARSE"
TOOL_DIR = ROOT / "library" / "qparse"
SCRIPT = TOOL_DIR / "QPARSE_1.10.py"


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def run_command(command, cwd):
    return subprocess.run(
        command,
        cwd=str(cwd),
        check=True,
        text=True,
        capture_output=True,
    )


def main():
    parser = argparse.ArgumentParser(description="Run QPARSE for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--min-len", type=int, default=2, help="Minimum island length (-m).")
    parser.add_argument("--max-loop", type=int, default=12, help="Maximum loop length (-L).")
    parser.add_argument("--island-num", type=int, default=4, help="Number of islands (-n).")
    parser.add_argument("--perfect", type=int, default=1, help="Minimum number of perfect islands (-p).")
    parser.add_argument("--output-format", choices=["txt", "tsv"], default="txt", help="Output file extension.")
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        output = out_dir / f"{fasta.stem}_qparse.{args.output_format}"

        command = [
            sys.executable,
            str(SCRIPT),
            "-i",
            str(fasta),
            "-o",
            str(output),
            "-m",
            str(args.min_len),
            "-L",
            str(args.max_loop),
            "-n",
            str(args.island_num),
            "-p",
            str(args.perfect),
        ]

        result = run_command(command, TOOL_DIR)
        if result.stderr:
            print(result.stderr, end="")
        print(result.stdout, end="")
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
