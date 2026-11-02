#!/usr/bin/env python3
"""Align common aptamers to all gene records found in the last Andes alignment."""

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from compare_common_sequences_to_best_genes import (
    clean_sequence,
    latest_common_sequences,
    read_tsv,
    safe_id,
    write_fasta,
    write_tsv,
    write_xlsx,
    rows_to_matrix,
    as_float,
    as_int,
    BLAST_COLUMNS,
    COMMON_DIR,
)


ROOT = Path(__file__).resolve().parent
ALIGNMENTS_DIR = ROOT / "alignments"
OUTPUT_DIR = ALIGNMENTS_DIR / "aptamers_vs_last_alignment_genes"


def read_fasta(path):
    title = None
    chunks = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if title is not None:
                    yield title, clean_sequence("".join(chunks))
                title = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if title is not None:
        yield title, clean_sequence("".join(chunks))


def load_db_records(path):
    records = {}
    for title, sequence in read_fasta(path):
        record_id = title.split()[0]
        records[record_id] = (title, sequence)
    return records


def wsl_path(path):
    path = Path(path).resolve()
    drive = path.drive.rstrip(":").lower()
    rest = path.as_posix().split(":", 1)[1].lstrip("/")
    return f"/mnt/{drive}/{rest}"


def run_blast(query_fasta, db_fasta, output):
    db_prefix = db_fasta.with_suffix("")
    command = (
        f"makeblastdb -in {wsl_path(db_fasta)!r} -dbtype nucl -parse_seqids -out {wsl_path(db_prefix)!r} >/dev/null && "
        f"blastn -query {wsl_path(query_fasta)!r} -db {wsl_path(db_prefix)!r} "
        "-task blastn-short "
        "-word_size 7 "
        "-dust no "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' "
        f"-out {wsl_path(output)!r}"
    )
    subprocess.run(["wsl", "bash", "-lc", command], check=True)


def read_blast(path):
    rows = []
    with path.open(newline="") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == len(BLAST_COLUMNS):
                rows.append(dict(zip(BLAST_COLUMNS, parts)))
    return rows


def species_from_title(title):
    title = re.sub(r"^NC_\S+\s+", "", title or "")
    title = re.split(r"\s+(segment|RNA|complete|strain|isolate|glycoprotein|nucleocapsid)\b", title, maxsplit=1)[0]
    return title.strip() or "Unknown"


def make_aptamer_queries(common_file):
    rows = []
    records = []
    for index, row in enumerate(read_tsv(common_file), 1):
        tools = row.get("tools", "")
        if "QGRS" not in tools or "PQS_Finder" not in tools:
            continue
        sequence = clean_sequence(row.get("sequence", ""))
        if not sequence:
            continue
        segment = row.get("input", "").strip()
        start = row.get("starts", "").split(";")[0]
        aptamer_id = f"aptamer_{index}|segment_{safe_id(segment)}|start_{safe_id(start)}|len_{len(sequence)}"
        metadata = {
            "aptamer_id": aptamer_id,
            "segment": segment,
            "andes_start": row.get("starts", ""),
            "andes_end": row.get("ends", ""),
            "sequence": sequence,
            "length_nt": len(sequence),
            "tools": tools,
        }
        rows.append(metadata)
        records.append((aptamer_id, f"common QGRS/PQS_Finder aptamer {index}", sequence))
    return rows, records


def make_gene_db_records(best_hits, alignments_dir):
    wanted_ids = defaultdict(set)
    hit_meta = {}
    for row in best_hits:
        segment = row["segment"]
        subject_id = row["subject_id"]
        wanted_ids[segment].add(subject_id)
        hit_meta[subject_id] = row

    fasta_records = []
    gene_meta = {}
    for segment, ids in sorted(wanted_ids.items()):
        db_fasta = alignments_dir / segment / f"{segment}_gene_records_db.fasta"
        records = load_db_records(db_fasta)
        for subject_id in sorted(ids):
            if subject_id not in records:
                continue
            title, sequence = records[subject_id]
            hit = hit_meta.get(subject_id, {})
            fasta_records.append((subject_id, title, sequence))
            gene_meta[subject_id] = {
                "gene_segment": segment,
                "gene_id": subject_id,
                "species": hit.get("species") or species_from_title(title),
                "gene_title": hit.get("gene_title") or title,
                "andes_identity_percent": hit.get("percent_identity", ""),
                "andes_alignment_length": hit.get("alignment_length", ""),
                "andes_query_coverage_percent": hit.get("query_coverage_percent", ""),
                "andes_bitscore": hit.get("bitscore", ""),
            }
    return fasta_records, gene_meta


def best_hit_per_aptamer_gene(blast_rows):
    best = {}
    for row in blast_rows:
        key = (row["qseqid"], row["sseqid"])
        score = (as_float(row["bitscore"]), as_int(row["length"]), as_float(row["pident"]))
        current = best.get(key)
        if current is None:
            best[key] = row
            continue
        current_score = (as_float(current["bitscore"]), as_int(current["length"]), as_float(current["pident"]))
        if score > current_score:
            best[key] = row
    return list(best.values())


def main():
    parser = argparse.ArgumentParser(description="Align common aptamers to all genes from last Andes alignment.")
    parser.add_argument("--common", type=Path, default=None)
    parser.add_argument("--alignments-dir", type=Path, default=ALIGNMENTS_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    common_file = args.common or latest_common_sequences(COMMON_DIR)
    best_hits_file = args.alignments_dir / "Andes_LMS_best_hits.tsv"
    best_hits = read_tsv(best_hits_file)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    aptamer_rows, aptamer_records = make_aptamer_queries(common_file)
    gene_records, gene_meta = make_gene_db_records(best_hits, args.alignments_dir)

    query_fasta = args.output_dir / "common_aptamers_qgrs_pqsfinder.fasta"
    db_fasta = args.output_dir / "genes_from_last_alignment.xlsx_source.fasta"
    raw_blast = args.output_dir / "aptamers_vs_last_alignment_genes_blast_raw.tsv"
    write_fasta(aptamer_records, query_fasta)
    write_fasta(gene_records, db_fasta)
    run_blast(query_fasta, db_fasta, raw_blast)

    aptamer_meta = {row["aptamer_id"]: row for row in aptamer_rows}
    best_rows = []
    for hit in best_hit_per_aptamer_gene(read_blast(raw_blast)):
        aptamer = aptamer_meta.get(hit["qseqid"], {})
        gene = gene_meta.get(hit["sseqid"], {})
        length = as_int(hit["length"])
        aptamer_len = int(aptamer.get("length_nt") or 0)
        best_rows.append({
            "aptamer_id": hit["qseqid"],
            "aptamer_segment": aptamer.get("segment", ""),
            "aptamer_sequence": aptamer.get("sequence", ""),
            "aptamer_length_nt": aptamer_len,
            "aptamer_andes_start": aptamer.get("andes_start", ""),
            "aptamer_andes_end": aptamer.get("andes_end", ""),
            "gene_segment": gene.get("gene_segment", ""),
            "gene_id": hit["sseqid"],
            "species": gene.get("species", ""),
            "gene_title": gene.get("gene_title", ""),
            "andes_vs_gene_identity_percent": gene.get("andes_identity_percent", ""),
            "aptamer_vs_gene_identity_percent": as_float(hit["pident"]),
            "aptamer_query_coverage_percent": round(length / aptamer_len * 100, 2) if aptamer_len else 0.0,
            "alignment_length": length,
            "aptamer_query_start": hit["qstart"],
            "aptamer_query_end": hit["qend"],
            "gene_start": hit["sstart"],
            "gene_end": hit["send"],
            "evalue": hit["evalue"],
            "bitscore": as_float(hit["bitscore"]),
            "aligned_aptamer": hit["qseq"],
            "aligned_gene": hit["sseq"],
        })

    fieldnames = [
        "aptamer_id", "aptamer_segment", "aptamer_sequence", "aptamer_length_nt",
        "aptamer_andes_start", "aptamer_andes_end", "gene_segment", "gene_id",
        "species", "gene_title", "andes_vs_gene_identity_percent",
        "aptamer_vs_gene_identity_percent", "aptamer_query_coverage_percent",
        "alignment_length", "aptamer_query_start", "aptamer_query_end", "gene_start",
        "gene_end", "evalue", "bitscore", "aligned_aptamer", "aligned_gene",
    ]
    selected_gene_fields = [
        "gene_segment", "gene_id", "species", "gene_title",
        "andes_identity_percent", "andes_alignment_length",
        "andes_query_coverage_percent", "andes_bitscore",
    ]
    selected_gene_rows = [gene_meta[key] for key in sorted(gene_meta)]

    write_tsv(args.output_dir / "aptamers_vs_all_last_alignment_genes_best_hits.tsv", best_rows, fieldnames)
    write_tsv(args.output_dir / "genes_used_from_last_alignment.tsv", selected_gene_rows, selected_gene_fields)

    by_aptamer = defaultdict(list)
    for row in best_rows:
        by_aptamer[row["aptamer_id"]].append(row)

    per_dir = args.output_dir / "per_aptamer"
    per_dir.mkdir(parents=True, exist_ok=True)
    for aptamer_id, rows in by_aptamer.items():
        label = safe_id(aptamer_id.split("|")[0] + "_" + aptamer_id.split("|")[1])
        write_xlsx(
            per_dir / f"{label}.xlsx",
            [("All gene alignments", rows_to_matrix(rows, fieldnames), [30, 10, 45, 12, 18, 18, 10, 24, 24, 70, 22, 24, 24, 16, 16, 16, 14, 14, 14, 12, 45, 45])],
        )

    write_xlsx(
        args.output_dir / "aptamers_vs_all_last_alignment_genes.xlsx",
        [
            ("Aptamer gene alignments", rows_to_matrix(best_rows, fieldnames), [30, 10, 45, 12, 18, 18, 10, 24, 24, 70, 22, 24, 24, 16, 16, 16, 14, 14, 14, 12, 45, 45]),
            ("Genes used", rows_to_matrix(selected_gene_rows, selected_gene_fields), [10, 24, 24, 70, 22, 22, 24, 14]),
        ],
    )

    print(f"Read {len(aptamer_rows)} aptamers from {common_file}")
    print(f"Used {len(gene_records)} genes from {best_hits_file}")
    print(f"Wrote {len(best_rows)} aptamer-gene best alignments to {args.output_dir}")


if __name__ == "__main__":
    main()
