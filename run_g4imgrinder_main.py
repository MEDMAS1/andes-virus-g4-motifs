#!/usr/bin/env python3
"""Core G4iMGrinder runner implementation."""

import argparse
import shutil
import subprocess
import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "G4iMGrinder"


def find_rscript():
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript
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
            return str(candidate)
    raise SystemExit(
        "Rscript not found. Install R and make sure Rscript is available on PATH."
    )


def check_rscript_and_packages(rscript):
    expr = (
        "options(warn=2);"
        "for (pkg in c('stringi', 'stringr', 'plyr')) {"
        "if (!requireNamespace(pkg, quietly = TRUE)) {"
        "stop(paste0('R package ', pkg, ' is required. Install it with install.packages(\\'', pkg, '\\') and rerun.'))}};"
    )
    try:
        run_command([rscript, "-e", expr], cwd=ROOT)
    except subprocess.CalledProcessError as exc:
        if exc.stdout:
            print(exc.stdout, end="")
        if exc.stderr:
            print(exc.stderr, end="")
        raise SystemExit("R package preflight failed.")


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def read_fasta_sequence(path):
    name = path.stem
    chunks = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0] or path.stem
            else:
                chunks.append(line)
    return name, "".join(chunks)


def safe_r_string(value):
    return value.replace("\\", "/").replace("'", "\\'")


def g4imgrinder_loader_expr():
    source_files = [
        "library/G4iMGrinder/R/G4iM.Grinder.SUBFuns.R",
        "library/G4iMGrinder/R/G4iM.Grinder.Funs.R",
    ]
    files_expr = "c(" + ",".join(f"'{path}'" for path in source_files) + ")"
    return (
        "suppressPackageStartupMessages(library(stringr));"
        "suppressPackageStartupMessages(library(stringi));"
        "suppressPackageStartupMessages(library(plyr));"
        "load('library/G4iMGrinder/data/GiG.DB.rda');"
        f"for (src in {files_expr}) {{"
        "code <- readLines(src, warn=FALSE);"
        "code <- gsub('G4iMGrinder:::', '', code, fixed=TRUE);"
        "code <- gsub('G4iMGrinder::GiG.DB', 'GiG.DB', code, fixed=TRUE);"
        "code <- gsub('packageDescription(\"G4iMGrinder\")$Version', '\"1.6.5\"', code, fixed=TRUE);"
        "eval(parse(text=code), envir=.GlobalEnv);"
        "};"
        ".PackageLoading <- function(...) invisible(TRUE);"
        ".G4Hscore <- function(y) {"
        "vals <- strsplit(as.character(y), '')[[1]];"
        "scores <- numeric(length(vals));"
        "i <- 1L;"
        "while (i <= length(vals)) {"
        "j <- i;"
        "while (j <= length(vals) && vals[j] == vals[i]) j <- j + 1L;"
        "run_len <- j - i;"
        "if (vals[i] == 'G') scores[i:(j-1L)] <- min(run_len, 4L);"
        "if (vals[i] == 'C') scores[i:(j-1L)] <- -min(run_len, 4L);"
        "i <- j;"
        "};"
        "mean(scores)"
        "};"
    )


def run_command(command, cwd):
    return subprocess.run(
        command,
        cwd=str(cwd),
        check=True,
        text=True,
        capture_output=True,
    )


def main():
    parser = argparse.ArgumentParser(description="Run G4iMGrinder for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--ncpu", type=int, default=1, help="Number of CPUs used by G4iMGrinder.")
    args = parser.parse_args()

    rscript = find_rscript()
    check_rscript_and_packages(rscript)

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        seqname, sequence = read_fasta_sequence(fasta)
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        rds_path = out_dir / f"{fasta.stem}_g4imgrinder.rds"
        csv_path = out_dir / f"{fasta.stem}_g4imgrinder_PQSM2a.tsv"
        csv_path_m3 = out_dir / f"{fasta.stem}_g4imgrinder_PQSM3a.tsv"

        sequence_str = safe_r_string(sequence)
        name_str = safe_r_string(fasta.stem)
        rds_str = safe_r_string(rds_path.as_posix())
        csv_str = safe_r_string(csv_path.as_posix())
        csv_m3_str = safe_r_string(csv_path_m3.as_posix())

        expression = (
            "options(warn=2);"
            + g4imgrinder_loader_expr()
            + "Sequence <- '" + sequence_str + "';"
            "Rs <- G4iMGrinder(Name='" + name_str + "', Sequence=Sequence, NCores=" + str(args.ncpu) + ", Verborrea=FALSE, KnownQuadruplex=FALSE, KnownNOTQuadruplex=FALSE);"
            "saveRDS(Rs, file='" + rds_str + "');"
            "if ('PQSM2a' %in% names(Rs)) write.table(Rs$PQSM2a, file='" + csv_str + "', sep='\\t', row.names=FALSE, quote=FALSE);"
            "if ('PQSM3a' %in% names(Rs)) write.table(Rs$PQSM3a, file='" + csv_m3_str + "', sep='\\t', row.names=FALSE, quote=FALSE);"
            "cat('done\\n')"
        )

        try:
            result = run_command([rscript, "-e", expression], cwd=ROOT)
        except subprocess.CalledProcessError as exc:
            if exc.stdout:
                print(exc.stdout, end="")
            if exc.stderr:
                print(exc.stderr, end="")
            raise SystemExit(f"G4iMGrinder failed for {fasta.name}")

        if result.stderr:
            print(result.stderr, end="")
        print(result.stdout, end="")

        if not rds_path.exists():
            raise SystemExit(f"G4iMGrinder did not write output file: {rds_path}")
        print(f"Wrote {rds_path}")
        if csv_path.exists():
            print(f"Wrote {csv_path}")
        if csv_path_m3.exists():
            print(f"Wrote {csv_path_m3}")


if __name__ == "__main__":
    main()
