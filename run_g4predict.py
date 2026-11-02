#!/usr/bin/env python3
"""Run G4Predict on FASTA files in data/ and save results under results/G4Predict."""

import argparse
import os
import platform
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "G4Predict"
TOOL_DIR = ROOT / "library" / "g4predict"
SCRIPT = TOOL_DIR / "g4funcs" / "g4predict.py"


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def run_command(command, cwd, env=None):
    return subprocess.run(
        command,
        cwd=str(cwd),
        check=True,
        text=True,
        capture_output=True,
        env=env,
    )


def wsl_path(path):
    path = Path(path).resolve()
    drive = path.drive.rstrip(":").lower()
    rest = path.as_posix().split(":", 1)[1].lstrip("/")
    return f"/mnt/{drive}/{rest}"


def run_g4predict(command, cwd, env):
    if platform.system() != "Windows":
        return run_command(command, cwd, env=env)

    wsl_tool_dir = wsl_path(TOOL_DIR)
    wsl_script = wsl_path(SCRIPT)
    wsl_pythonpath = wsl_tool_dir
    wsl_command = [
        "cd",
        shlex.quote(wsl_tool_dir),
        "&&",
        f"PYTHONPATH={shlex.quote(wsl_pythonpath)}",
        "python3",
        shlex.quote(wsl_script),
    ]
    wsl_command.extend(shlex.quote(str(part)) for part in command[2:])
    return subprocess.run(
        ["wsl", "-e", "bash", "-lc", " ".join(wsl_command)],
        cwd=str(cwd),
        check=True,
        text=True,
        capture_output=True,
    )


def writable_output(path):
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a"):
            pass
        return path
    except PermissionError:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return path.with_name(f"{path.stem}_{stamp}{path.suffix}")


def main():
    parser = argparse.ArgumentParser(description="Run G4Predict for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--mode", choices=["intra", "inter"], default="intra", help="Prediction mode.")
    parser.add_argument("--min-tetrad", type=int, default=2, help="Minimum tetrads. Default 2 matches G2-rich viral calls.")
    parser.add_argument("--max-tetrad", type=int, default=3, help="Maximum tetrads.")
    parser.add_argument("--min-loop", default="1", help="Minimum loop length(s).")
    parser.add_argument("--max-loop", default="12", help="Maximum loop length(s). Default 12 matches QGRS/Quadron-style windows.")
    parser.add_argument("--write-bed6", action="store_true", help="Write BED6 output instead of BED12.")
    parser.add_argument("--filter-overlapping", action="store_true", help="Filter overlapping predictions.")
    parser.add_argument("--merge-overlapping", action="store_true", help="Merge overlapping predictions.")
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    env = os.environ.copy()
    env["PYTHONPATH"] = str(TOOL_DIR)

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        output = writable_output(out_dir / f"{fasta.stem}_g4predict.bed")

        command = [
            sys.executable,
            str(SCRIPT),
            args.mode,
            "-f",
            wsl_path(fasta) if platform.system() == "Windows" else str(fasta),
            "-b",
            wsl_path(output) if platform.system() == "Windows" else str(output),
            "-tmin",
            str(args.min_tetrad),
            "-tmax",
            str(args.max_tetrad),
            "-lmin",
            str(args.min_loop),
            "-lmax",
            str(args.max_loop),
        ]
        if args.write_bed6:
            command.append("-s")
        if args.filter_overlapping:
            command.append("-F")
        if args.merge_overlapping:
            command.append("-M")

        result = run_g4predict(command, TOOL_DIR, env)
        if result.stderr:
            print(result.stderr, end="")
        print(result.stdout, end="")
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
