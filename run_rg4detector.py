#!/usr/bin/env python3
"""Run rG4detector on FASTA files in data/."""

import argparse
import os
import platform
import shlex
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "rG4detector"
TOOL_DIR = ROOT / "library" / "rG4detector"
SCRIPT = TOOL_DIR / "code" / "predict_fasta.py"
MODEL_DIR = TOOL_DIR / "model"


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


def make_dna_copy(src, dst):
    text = src.read_text()
    lines = []
    for line in text.splitlines():
        if line.startswith(">"):
            lines.append(line)
        else:
            lines.append(line.replace("U", "T").replace("u", "t"))
    dst.write_text("\n".join(lines) + "\n")
    return dst


def main():
    parser = argparse.ArgumentParser(description="Run rG4detector for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--model-dir", type=Path, default=MODEL_DIR)
    parser.add_argument("--python", default=os.environ.get("PYTHON3", "python3"))
    parser.add_argument("--detect", action="store_true", help="Use detection mode instead of prediction mode.")
    parser.add_argument("--plot", action="store_true", help="Save plots in detection mode.")
    parser.add_argument("--ensemble-size", type=int, default=11)
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        dna_fasta = make_dna_copy(fasta, out_dir / f"{fasta.stem}_dna_for_rg4detector.fasta")

        if platform.system() == "Windows":
            script = wsl_path(SCRIPT)
            fasta_arg = wsl_path(dna_fasta)
            out_arg = wsl_path(out_dir)
            model_arg = wsl_path(args.model_dir)
            cwd = wsl_path(TOOL_DIR)
        else:
            script = str(SCRIPT)
            fasta_arg = str(dna_fasta)
            out_arg = str(out_dir)
            model_arg = str(args.model_dir)
            cwd = str(TOOL_DIR)

        parts = [
            "cd", shlex.quote(cwd), "&&",
            shlex.quote(args.python), shlex.quote(script),
            "-f", shlex.quote(fasta_arg),
            "-m", shlex.quote(model_arg),
            "-o", shlex.quote(out_arg),
            "-s", str(args.ensemble_size),
        ]
        if args.detect:
            parts.append("-d")
        if args.plot:
            parts.append("-p")

        run_command(" ".join(parts))
        print(f"Wrote {out_dir}")


if __name__ == "__main__":
    main()
