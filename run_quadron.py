#!/usr/bin/env python3
"""Run Quadron on FASTA files in data/ and save results under results/Quadron."""

import argparse
import os
import shutil
import subprocess
import sys
import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "Quadron"
TOOL_DIR = ROOT / "library" / "Quadron"
QUADRON_LIB = TOOL_DIR / "Quadron.lib"
WSL_R_LIB = ROOT / ".r-lib-wsl"
WINDOWS_R_LIB = ROOT / ".r-lib"


def is_wsl():
    if platform.system() != "Linux":
        return False
    try:
        version = Path("/proc/version").read_text(encoding="utf-8").lower()
    except OSError:
        return False
    return "microsoft" in version or "wsl" in version


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def find_rscript():
    if platform.system() == "Windows" and shutil.which("wsl"):
        probe = subprocess.run(
            ["wsl", "-e", "which", "Rscript"],
            text=True,
            capture_output=True,
        )
        if probe.returncode == 0 and probe.stdout.strip():
            return ("wsl", probe.stdout.strip())

    rscript = shutil.which("Rscript")
    if rscript:
        return ("native", rscript)
    system = platform.system()
    if system == "Windows":
        candidates = [
            Path("C:/Program Files/R/R-4.6.0/bin/Rscript.exe"),
            Path("C:/Program Files/R/R-4.5.0/bin/Rscript.exe"),
            Path("C:/Program Files/R/R-4.4.0/bin/Rscript.exe"),
        ]
    else:
        candidates = [
            Path("/usr/bin/Rscript"),
            Path("/usr/local/bin/Rscript"),
            Path("/opt/R/bin/Rscript"),
        ]
    for candidate in candidates:
        if candidate.exists():
            return ("native", str(candidate))
    raise SystemExit(
        "Rscript not found. Install R and make sure Rscript is available on PATH."
    )


def to_wsl_path(path):
    resolved = Path(path).resolve()
    drive = resolved.drive.rstrip(":").lower()
    parts = resolved.parts[1:]
    return "/mnt/" + drive + "/" + "/".join(part.replace("\\", "/") for part in parts)


def r_quote(value):
    return "'" + str(value).replace("\\", "/").replace("'", "\\'") + "'"


def run_command(command, cwd, env=None):
    result = subprocess.run(
        command,
        cwd=str(cwd),
        text=True,
        capture_output=True,
        env=env,
    )
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, end="", file=sys.stderr)
        raise SystemExit(f"Command failed with exit code {result.returncode}: {' '.join(command)}")
    return result


def build_r_command(rscript_kind, rscript, fasta, output, ncpu):
    if rscript_kind == "wsl":
        tool_dir = to_wsl_path(TOOL_DIR)
        fasta_path = to_wsl_path(fasta)
        output_path = to_wsl_path(output)
        r_lib = to_wsl_path(WSL_R_LIB)
        expression = (
            f".libPaths(c({r_quote(r_lib)}, .libPaths()));"
            f"setwd({r_quote(tool_dir)});"
            f"load('Quadron.lib');"
            f"Quadron(FastaFile={r_quote(fasta_path)}, OutFile={r_quote(output_path)}, nCPU={ncpu})"
        )
        return ["wsl", "-e", "env", f"R_LIBS_USER={r_lib}", "Rscript", "-e", expression], None

    r_lib = WSL_R_LIB if is_wsl() else WINDOWS_R_LIB
    expression = (
        f".libPaths(c({r_quote(r_lib.as_posix())}, .libPaths()));"
        f"setwd({r_quote(TOOL_DIR.as_posix())});"
        f"load('Quadron.lib');"
        f"Quadron(FastaFile={r_quote(fasta.as_posix())}, OutFile={r_quote(output.as_posix())}, nCPU={ncpu})"
    )
    env = os.environ.copy()
    env["R_LIBS_USER"] = str(r_lib)
    return [rscript, "-e", expression], env


def main():
    parser = argparse.ArgumentParser(description="Run Quadron for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--ncpu", type=int, default=1, help="Number of CPUs used by Quadron.")
    args = parser.parse_args()

    if not QUADRON_LIB.exists():
        raise SystemExit(
            "Quadron cannot run because library/Quadron/Quadron.lib is missing. "
            "Please build or install Quadron before using run_quadron.py."
        )

    rscript_kind, rscript = find_rscript()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        output = out_dir / f"{fasta.stem}_quadron.txt"

        command, env = build_r_command(rscript_kind, rscript, fasta, output, args.ncpu)

        result = run_command(command, TOOL_DIR, env=env)
        if result.stderr:
            print(result.stderr, end="")
        print(result.stdout, end="")
        if not output.exists():
            raise SystemExit(f"Quadron did not write output file: {output}")
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
