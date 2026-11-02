#!/usr/bin/env python3
"""Run G4RNA Screener on FASTA files in data/."""

import argparse
import platform
import shlex
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "G4RNA_Screener"
TOOL_DIR = ROOT / "library" / "g4rna_screener"
PYTHON2 = ROOT / "library" / "g4rna_env" / "bin" / "python"
SCREEN = TOOL_DIR / "screen.py"


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def wsl_path(path):
    path = Path(path).resolve()
    drive = path.drive.rstrip(":").lower()
    rest = path.as_posix().split(":", 1)[1].lstrip("/")
    return f"/mnt/{drive}/{rest}"


def run_command(command):
    if platform.system() == "Windows":
        subprocess.run(["wsl", "bash", "-lc", command], check=True)
    else:
        subprocess.run(command, shell=True, check=True)


def main():
    parser = argparse.ArgumentParser(description="Run G4RNA Screener for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--window", type=int, default=60)
    parser.add_argument("--step", type=int, default=10)
    parser.add_argument(
        "--columns",
        nargs="+",
        default=["description", "sequence", "start", "cGcC", "G4H", "G4NN"],
    )
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        output = out_dir / f"{fasta.stem}_g4rna_screener.tsv"

        if platform.system() == "Windows":
            python2 = wsl_path(PYTHON2)
            screen = wsl_path(SCREEN)
            fasta_arg = wsl_path(fasta)
            output_arg = wsl_path(output)
            cwd = wsl_path(TOOL_DIR)
        else:
            python2 = str(PYTHON2)
            screen = str(SCREEN)
            fasta_arg = str(fasta)
            output_arg = str(output)
            cwd = str(TOOL_DIR)

        parts = [
            "cd", shlex.quote(cwd), "&&",
            shlex.quote(python2), shlex.quote(screen),
            shlex.quote(fasta_arg),
            "-w", str(args.window),
            "-s", str(args.step),
            "-c", *[shlex.quote(column) for column in args.columns],
            ">", shlex.quote(output_arg),
        ]
        run_command(" ".join(parts))
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
