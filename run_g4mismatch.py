#!/usr/bin/env python3
"""Run G4mismatch scan mode on FASTA files in data/."""

import argparse
import os
import platform
import shlex
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "G4mismatch"
TOOL_DIR = ROOT / "library" / "G4mismatch"
SCRIPT = TOOL_DIR / "G4mismatch_scan.py"
DEFAULT_MODEL = TOOL_DIR / "models" / "model_K.h5"


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
    parser = argparse.ArgumentParser(description="Run G4mismatch scan mode for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--python", default=os.environ.get("PYTHON3", "python3"))
    parser.add_argument("--flank-size", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument(
        "--filter-bins",
        type=float,
        default=None,
        help="Only keep bins at or above this score. Default keeps unfiltered predictions.",
    )
    parser.add_argument("--merge", type=int, choices=[0, 1], default=1)
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        dna_fasta = make_dna_copy(fasta, out_dir / f"{fasta.stem}_dna_for_g4mismatch.fasta")

        if platform.system() == "Windows":
            script = wsl_path(SCRIPT)
            fasta_arg = wsl_path(dna_fasta)
            out_arg = wsl_path(out_dir)
            model_arg = wsl_path(args.model)
            cwd = wsl_path(TOOL_DIR)
        else:
            script = str(SCRIPT)
            fasta_arg = str(dna_fasta)
            out_arg = str(out_dir)
            model_arg = str(args.model)
            cwd = str(TOOL_DIR)

        parts = [
            "cd", shlex.quote(cwd), "&&",
            shlex.quote(args.python), shlex.quote(script),
            "-dp", shlex.quote(fasta_arg),
            "-mp", shlex.quote(model_arg),
            "-of", shlex.quote(out_arg),
            "-obg", shlex.quote(f"{fasta.stem}_g4mismatch"),
            "-f", str(args.flank_size),
            "-bs", str(args.batch_size),
            "-m", str(args.merge),
        ]
        if args.filter_bins is not None:
            parts.extend(["-fb", str(args.filter_bins)])
        run_command(" ".join(parts))
        print(f"Wrote {out_dir}")


if __name__ == "__main__":
    main()
