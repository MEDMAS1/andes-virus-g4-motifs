#!/usr/bin/env python3
"""Run G4Catchall on FASTA files in data/ and save one result folder per input."""

import argparse
import platform
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "G4Catchall"
TOOL_DIR = ROOT / "library" / "G4Catchall"
SCRIPT = TOOL_DIR / "G4Catchall.py"


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def run_command(command, cwd):
    if platform.system() == "Windows":
        return subprocess.run(
            command,
            cwd=str(cwd),
            shell=True,
            check=True,
            text=True,
            capture_output=True,
        )
    return subprocess.run(command, cwd=str(cwd), shell=True, check=True, text=True, capture_output=True)


def main():
    parser = argparse.ArgumentParser(description="Run G4Catchall for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--g3l", type=str, default="1..8", help="G3 loop limits passed to --G3L.")
    parser.add_argument("--g2l", type=str, help="G2 loop limits passed to --G2L.")
    parser.add_argument("--g4h-threshold", type=float, default=0.0, help="Minimum absolute G4Hunter score passed to --G4HThreshold.")
    parser.add_argument("--no-overlap", action="store_true", help="Pass --dont_merge_overlapping.")
    parser.add_argument("--no-reverse", action="store_true", help="Pass --no_reverse.")
    parser.add_argument("--include-flanks", action="store_true", help="Pass --include_flanks.")
    parser.add_argument("--extreme-loop", nargs="?", const="1..30", help="Pass --extreme_loop [limits].")
    parser.add_argument("--min-gtract-extreme-loop", type=int, choices=[2, 3], default=3, help="Pass --min_Gtract_for_extreme_loop.")
    parser.add_argument("--bulge-only", action="store_true", help="Pass --bulge_only.")
    parser.add_argument("--max-length", type=int, help="Pass --max.")
    parser.add_argument("--g2gqs-allowed", action="store_true", help="Pass --G2GQs_allowed.")
    parser.add_argument("--max-imperfect-gtacts", type=int, default=1, help="Pass --max_imperfect_Gtracts.")
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    python_exe = shlex.quote(sys.executable)
    script_path = shlex.quote(str(SCRIPT))

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        output = out_dir / f"{fasta.stem}_g4catchall.tsv"

        cmd_parts = [python_exe, script_path, "-f", shlex.quote(str(fasta)), "--G3L", args.g3l]
        if args.g2l:
            cmd_parts.extend(["--G2L", args.g2l])
        if args.g4h_threshold is not None:
            cmd_parts.extend(["--G4HThreshold", str(args.g4h_threshold)])
        if args.no_overlap:
            cmd_parts.append("--dont_merge_overlapping")
        if args.no_reverse:
            cmd_parts.append("--no_reverse")
        if args.include_flanks:
            cmd_parts.append("--include_flanks")
        if args.extreme_loop is not None:
            cmd_parts.append("--extreme_loop")
            if args.extreme_loop != "1..30":
                cmd_parts.append(args.extreme_loop)
        if args.min_gtract_extreme_loop is not None:
            cmd_parts.extend(["--min_Gtract_for_extreme_loop", str(args.min_gtract_extreme_loop)])
        if args.bulge_only:
            cmd_parts.append("--bulge_only")
        if args.max_length is not None:
            cmd_parts.extend(["--max", str(args.max_length)])
        if args.g2gqs_allowed:
            cmd_parts.append("--G2GQs_allowed")
        if args.max_imperfect_gtacts is not None:
            cmd_parts.extend(["--max_imperfect_Gtracts", str(args.max_imperfect_gtacts)])

        command = " ".join(cmd_parts)
        result = run_command(command, TOOL_DIR)
        if result.stderr:
            print(result.stderr, end="")
        output.write_text(result.stdout)
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
