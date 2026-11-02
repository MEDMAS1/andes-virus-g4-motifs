#!/usr/bin/env python3
"""Run QGRS on FASTA files in data/ and save one result folder per input."""

import argparse
import os
import platform
import shlex
import subprocess
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "QGRS"
TOOL_DIR = ROOT / "library" / "qgrs-cpp"
QGRS_BIN = TOOL_DIR / "bin" / "qgrs"


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def read_fasta_sequence(path):
    chunks = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            chunks.append(line)
    return "".join(chunks).upper()


def wsl_path(path):
    path = Path(path).resolve()
    drive = path.drive.rstrip(":").lower()
    rest = path.as_posix().split(":", 1)[1].lstrip("/")
    return f"/mnt/{drive}/{rest}"


def run_command(command, cwd):
    if platform.system() == "Windows":
        return subprocess.run(
            ["wsl", "bash", "-lc", command],
            check=True,
            text=True,
            capture_output=True,
        )
    else:
        return subprocess.run(command, cwd=str(cwd), shell=True, check=True, text=True, capture_output=True)


def writable_output(path):
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a"):
            pass
        return path
    except PermissionError:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return path.with_name(f"{path.stem}_{stamp}{path.suffix}")


def write_website_format(raw_output, output):
    rows = []
    for line in raw_output.splitlines():
        line = line.strip()
        if not line or line.startswith("ID") or line.startswith("-"):
            continue
        parts = line.split(None, 7)
        if len(parts) < 8:
            continue
        _, t1, _t2, _t3, _t4, _tetrads, g_score, seq = parts
        rows.append((int(t1) + 1, len(seq), seq, g_score))

    output = writable_output(output)
    with output.open("w", newline="") as handle:
        handle.write("Position\tLength\tQGRS\tG-Score\n")
        for row in rows:
            handle.write("{}\t{}\t{}\t{}\n".format(*row))
    return output


def main():
    parser = argparse.ArgumentParser(description="Run QGRS for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--tetrads", type=int, default=2, help="Minimum tetrads passed to QGRS -t.")
    parser.add_argument("--score", type=int, default=0, help="Minimum G-score passed to QGRS -s.")
    parser.add_argument("--csv", action="store_true", help="Write CSV instead of the default text output.")
    parser.add_argument(
        "--no-overlap",
        action="store_true",
        help="Do not include overlapping QGRS results. By default all overlapping QGRS are included.",
    )
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        suffix = "csv" if args.csv else "txt"
        output = out_dir / f"{fasta.stem}_qgrs.{suffix}"
        website_output = out_dir / f"{fasta.stem}_qgrs_website_format.tsv"
        clean_input = out_dir / f"{fasta.stem}_qgrs_input_sequence.txt"
        clean_input.write_text(read_fasta_sequence(fasta) + "\n")

        if platform.system() == "Windows":
            qgrs = wsl_path(QGRS_BIN)
            fasta_arg = wsl_path(clean_input)
            output_arg = wsl_path(output)
            cwd_arg = wsl_path(TOOL_DIR)
        else:
            qgrs = str(QGRS_BIN)
            fasta_arg = str(clean_input)
            output_arg = str(output)
            cwd_arg = str(TOOL_DIR)

        parts = [
            shlex.quote(qgrs),
            "-i", shlex.quote(fasta_arg),
            "-t", str(args.tetrads),
            "-s", str(args.score),
        ]
        if args.csv:
            parts.append("-csv")
        if not args.no_overlap:
            parts.append("-v")

        result = run_command(" ".join(parts), cwd_arg)
        if result.stderr:
            print(result.stderr, end="")
        output.write_text(result.stdout)
        if not args.csv:
            website_output = write_website_format(result.stdout, website_output)
        print(f"Wrote {output}")
        if not args.csv:
            print(f"Wrote {website_output}")


if __name__ == "__main__":
    main()
