#!/usr/bin/env python3
"""Run G4Hunter on FASTA files in data/."""

import argparse
import csv
import platform
import shlex
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results" / "G4Hunter"
TOOL_DIR = ROOT / "library" / "G4-Hunter-master"
PYTHON2 = ROOT / "library" / "g4rna_env" / "bin" / "python"
G4HUNTER = TOOL_DIR / "G4Hunter.py"


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
    return name, "".join(chunks).upper()


def reverse_complement(seq):
    table = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return seq.translate(table)[::-1].upper()


def base_scores(seq):
    scores = []
    i = 0
    while i < len(seq):
        base = seq[i].upper()
        if base == "G":
            j = i
            while j < len(seq) and seq[j].upper() == "G":
                j += 1
            value = min(j - i, 4)
            scores.extend([value] * (j - i))
            i = j
        elif base == "C":
            j = i
            while j < len(seq) and seq[j].upper() == "C":
                j += 1
            value = -min(j - i, 4)
            scores.extend([value] * (j - i))
            i = j
        else:
            scores.append(0)
            i += 1
    return scores


def g4hunter_sequence_score(seq):
    scores = base_scores(seq)
    return sum(scores) / float(len(scores)) if scores else 0.0


def parse_window_file(path):
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">") or line.startswith("Start"):
                continue
            parts = [part.strip() for part in line.split("\t")]
            parts = [part for part in parts if part]
            if len(parts) < 5:
                continue
            start, end, seq, length, score = parts[:5]
            rows.append(
                {
                    "start0": int(start),
                    "end0": int(end),
                    "sequence": seq,
                    "length": int(length),
                    "score": float(score),
                }
            )
    return rows


def make_website_rows(seqname, full_sequence, windows, threshold, window_size):
    hits = []
    cluster = []
    for row in windows:
        if abs(row["score"]) >= threshold:
            if cluster and row["start0"] > cluster[-1]["end0"]:
                hits.append(cluster)
                cluster = []
            cluster.append(row)
        elif cluster:
            hits.append(cluster)
            cluster = []
    if cluster:
        hits.append(cluster)

    rows = []
    for cluster in hits:
        union_start0 = cluster[0]["start0"]
        union_end0 = cluster[-1]["end0"]
        raw_seq = full_sequence[union_start0:union_end0]
        scores = [row["score"] for row in cluster]
        max_score = max(scores, key=abs)
        strand = "+" if max_score >= 0 else "-"
        oriented_seq = raw_seq if strand == "+" else reverse_complement(raw_seq)
        first_g = oriented_seq.find("G")
        last_g = oriented_seq.rfind("G")
        if first_g == -1 or last_g == -1 or first_g > last_g:
            continue
        reported_seq = oriented_seq[first_g:last_g + 1]
        right_trim = len(oriented_seq) - last_g - 1
        if strand == "+":
            start = union_start0 + first_g + 1
            end = union_end0 - right_trim
        else:
            start = union_start0 + right_trim + 1
            end = union_end0 - first_g
        seq_score = g4hunter_sequence_score(reported_seq)
        if strand == "-":
            seq_score = -abs(seq_score)
        rows.append(
            {
                "seqnames": seqname,
                "start": start,
                "end": end,
                "width": len(reported_seq),
                "strand": strand,
                "score": "{:.2f}".format(seq_score),
                "max_score": "{:.2f}".format(max_score),
                "threshold": "{:.2f}".format(threshold),
                "window": "{:.2f}".format(float(window_size)),
                "sequence": reported_seq,
            }
        )
    return rows


def write_website_format(output, rows):
    columns = [
        "seqnames",
        "start",
        "end",
        "width",
        "strand",
        "score",
        "max_score",
        "threshold",
        "window",
        "sequence",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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
    parser = argparse.ArgumentParser(description="Run G4Hunter for each FASTA file.")
    parser.add_argument("--input", type=Path, help="One FASTA file to process. Default: all files in data/.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--window", type=int, default=25)
    parser.add_argument("--score", type=float, default=1.2)
    args = parser.parse_args()

    inputs = [args.input] if args.input else fasta_files(args.data_dir)
    if not inputs:
        raise SystemExit(f"No FASTA files found in {args.data_dir}")

    for fasta in inputs:
        fasta = fasta.resolve()
        out_dir = (args.results_dir / fasta.stem).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        seqname, full_sequence = read_fasta_sequence(fasta)

        if platform.system() == "Windows":
            python2 = wsl_path(PYTHON2)
            script = wsl_path(G4HUNTER)
            fasta_arg = wsl_path(fasta)
            out_arg = wsl_path(out_dir)
            cwd = wsl_path(TOOL_DIR)
        else:
            python2 = str(PYTHON2)
            script = str(G4HUNTER)
            fasta_arg = str(fasta)
            out_arg = str(out_dir)
            cwd = str(TOOL_DIR)

        parts = [
            "cd", shlex.quote(cwd), "&&",
            shlex.quote(python2), shlex.quote(script),
            "-i", shlex.quote(fasta_arg),
            "-o", shlex.quote(out_arg),
            "-w", str(args.window),
            "-s", str(args.score),
        ]
        run_command(" ".join(parts))
        window_file = out_dir / "Results" / f"{fasta.stem}-{args.window}-nts"
        website_file = out_dir / "Results" / f"{fasta.stem}_g4hunter_website_format.tsv"
        if window_file.exists():
            windows = parse_window_file(window_file)
            rows = make_website_rows(seqname, full_sequence, windows, args.score, args.window)
            write_website_format(website_file, rows)
            print(f"Wrote {website_file}")
        print(f"Wrote {out_dir / 'Results'}")


if __name__ == "__main__":
    main()
