#!/usr/bin/env python3
"""Find canonical potential quadruplex-forming sequences (PQS) in FASTA files.

This lightweight runner scans FASTA files from data/ with the canonical PQS
pattern G-run + loop + G-run + loop + G-run + loop + G-run. It is useful when
R/Bioconductor pqsfinder is not installed. It does not reproduce pqsfinder's
full imperfection-tolerant scoring model.
"""

import argparse
import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "PQS_Finder"


def fasta_files(data_dir):
    files = []
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        files.extend(data_dir.glob(pattern))
    return sorted(files)


def read_fasta(path):
    name = None
    chunks = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip() or path.stem
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def clean_sequence(seq):
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def reverse_complement(seq):
    normalized = seq.upper().replace("U", "T")
    table = str.maketrans("ACGTN", "TGCAN")
    return normalized.translate(table)[::-1]


def pqs_pattern(run_min, run_max, loop_min, loop_max):
    run = f"G{{{run_min},{run_max}}}"
    loop = f"[ACGTUN]{{{loop_min},{loop_max}}}"
    return re.compile(f"(?=(({run})({loop})({run})({loop})({run})({loop})({run})))")


def find_pqs(seq, pattern, strand):
    for match in pattern.finditer(seq):
        pqs = match.group(1)
        start = match.start(1) + 1
        end = start + len(pqs) - 1
        runs = [match.group(i) for i in (2, 4, 6, 8)]
        loops = [match.group(i) for i in (3, 5, 7)]
        tetrads = min(len(run) for run in runs)
        yield {
            "start": start,
            "end": end,
            "strand": strand,
            "length": len(pqs),
            "tetrads": tetrads,
            "run1_len": len(runs[0]),
            "loop1_len": len(loops[0]),
            "run2_len": len(runs[1]),
            "loop2_len": len(loops[1]),
            "run3_len": len(runs[2]),
            "loop3_len": len(loops[2]),
            "run4_len": len(runs[3]),
            "sequence": pqs,
        }


def main():
    parser = argparse.ArgumentParser(description="Find canonical PQS motifs in FASTA files.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--run-min", type=int, default=2, help="Minimum G-run length.")
    parser.add_argument("--run-max", type=int, default=11, help="Maximum G-run length.")
    parser.add_argument("--loop-min", type=int, default=0, help="Minimum loop length.")
    parser.add_argument("--loop-max", type=int, default=30, help="Maximum loop length.")
    parser.add_argument(
        "--strand",
        choices=["+", "-", "both"],
        default="+",
        help="Scan forward strand, reverse-complement strand, or both.",
    )
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    pattern = pqs_pattern(args.run_min, args.run_max, args.loop_min, args.loop_max)
    columns = [
        "record",
        "start",
        "end",
        "strand",
        "length",
        "tetrads",
        "run1_len",
        "loop1_len",
        "run2_len",
        "loop2_len",
        "run3_len",
        "loop3_len",
        "run4_len",
        "sequence",
    ]

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        output = out_dir / f"{fasta.stem}_pqs_finder.tsv"

        count = 0
        with output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            for record, raw_seq in read_fasta(fasta):
                seq = clean_sequence(raw_seq)
                strands = []
                if args.strand in ("+", "both"):
                    strands.append(("+", seq))
                if args.strand in ("-", "both"):
                    strands.append(("-", reverse_complement(seq)))
                for strand, scan_seq in strands:
                    for row in find_pqs(scan_seq, pattern, strand):
                        row["record"] = record
                        writer.writerow(row)
                        count += 1

        print(f"Wrote {output} ({count} PQS)")


if __name__ == "__main__":
    main()
