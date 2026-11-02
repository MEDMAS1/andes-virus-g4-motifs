#!/usr/bin/env python3
"""BLAST Andes virus L/M/S gene records against all gene records and export Excel reports."""

import argparse
import csv
import html
import re
import subprocess
import zipfile
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parent
GENE_DIR = ROOT / "gene"
OUTPUT_DIR = ROOT / "alignments"

BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qseq",
    "sseq",
]


def clean_sequence(sequence):
    return re.sub(r"[^A-Za-z]", "", sequence or "").upper().replace("U", "T")


def safe_id(value):
    value = re.sub(r"[^A-Za-z0-9_.:-]+", "_", value.strip())
    return value.strip("_") or "NA"


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


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


def segment_from_path(path):
    name = path.stem.upper()
    if name.startswith("L"):
        return "L"
    if name.startswith("M"):
        return "M"
    if name.startswith("S"):
        return "S"
    return safe_id(path.stem)


def species_from_title(title):
    title = re.sub(r"^NC_\S+\s+", "", title or "")
    title = re.split(r"\s+(segment|RNA|complete|strain|isolate|glycoprotein|nucleocapsid)\b", title, maxsplit=1)[0]
    return title.strip() or "Unknown"


def wsl_path(path):
    path = Path(path).resolve()
    drive = path.drive.rstrip(":").lower()
    rest = path.as_posix().split(":", 1)[1].lstrip("/")
    return f"/mnt/{drive}/{rest}"


def write_fasta(records, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="\n") as handle:
        for record_id, title, sequence in records:
            handle.write(f">{record_id} {title}\n{sequence}\n")


def run_blast(query_fasta, db_fasta, output):
    db_prefix = db_fasta.with_suffix("")
    command = (
        f"makeblastdb -in {wsl_path(db_fasta)!r} -dbtype nucl -parse_seqids -out {wsl_path(db_prefix)!r} >/dev/null && "
        f"blastn -query {wsl_path(query_fasta)!r} -db {wsl_path(db_prefix)!r} "
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


def write_tsv(path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def col_letter(index):
    letters = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


def cell_xml(row_index, col_index, value):
    ref = f"{col_letter(col_index)}{row_index}"
    if value is None:
        value = ""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f'<c r="{ref}"><v>{value}</v></c>'
    text = html.escape(str(value), quote=False)
    return f'<c r="{ref}" t="inlineStr"><is><t>{text}</t></is></c>'


def worksheet_xml(rows, widths=None):
    widths = widths or []
    cols = ""
    if widths:
        cols = "<cols>" + "".join(
            f'<col min="{idx}" max="{idx}" width="{width}" customWidth="1"/>'
            for idx, width in enumerate(widths, 1)
        ) + "</cols>"
    sheet_rows = []
    for row_index, row in enumerate(rows, 1):
        cells = "".join(cell_xml(row_index, col_index, value) for col_index, value in enumerate(row, 1))
        sheet_rows.append(f'<row r="{row_index}">{cells}</row>')
    max_cols = max((len(row) for row in rows), default=1)
    max_rows = max(len(rows), 1)
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
        f'{cols}<sheetData>{"".join(sheet_rows)}</sheetData>'
        f'<autoFilter ref="A1:{col_letter(max_cols)}{max_rows}"/>'
        "</worksheet>"
    )


def write_xlsx(path, sheets):
    path.parent.mkdir(parents=True, exist_ok=True)
    sheet_names = [re.sub(r"[\[\]\:*?/\\]", "_", name)[:31] or "Sheet1" for name, _rows, _widths in sheets]
    sheet_count = len(sheets)
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        + "".join(
            f'<Override PartName="/xl/worksheets/sheet{idx}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            for idx in range(1, sheet_count + 1)
        )
        + "</Types>"
    )
    root_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
        + "".join(
            f'<sheet name="{html.escape(name, quote=True)}" sheetId="{idx}" r:id="rId{idx}"/>'
            for idx, name in enumerate(sheet_names, 1)
        )
        + "</sheets></workbook>"
    )
    workbook_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(
            f'<Relationship Id="rId{idx}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{idx}.xml"/>'
            for idx in range(1, sheet_count + 1)
        )
        + f'<Relationship Id="rId{sheet_count + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
        + "</Relationships>"
    )
    styles = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
        '<borders count="1"><border/></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        "</styleSheet>"
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", root_rels)
        zf.writestr("xl/workbook.xml", workbook)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        zf.writestr("xl/styles.xml", styles)
        for idx, (_name, rows, widths) in enumerate(sheets, 1):
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", worksheet_xml(rows, widths))
    verify_xlsx(path)


def verify_xlsx(path):
    with zipfile.ZipFile(path) as zf:
        ET.fromstring(zf.read("xl/workbook.xml"))
        for name in zf.namelist():
            if name.startswith("xl/worksheets/sheet"):
                ET.fromstring(zf.read(name))


def enrich_hits(segment, rows, id_to_title, query_length):
    enriched = []
    for row in rows:
        title = id_to_title.get(row["sseqid"], "")
        q_start = as_int(row["qstart"])
        q_end = as_int(row["qend"])
        if q_start > q_end:
            q_start, q_end = q_end, q_start
        coverage = as_float(row["length"]) / query_length * 100 if query_length else 0.0
        enriched.append({
            "segment": segment,
            "query_id": row["qseqid"],
            "subject_id": row["sseqid"],
            "species": species_from_title(title),
            "gene_title": title,
            "percent_identity": as_float(row["pident"]),
            "alignment_length": as_int(row["length"]),
            "query_coverage_percent": round(coverage, 2),
            "query_start": q_start,
            "query_end": q_end,
            "subject_start": row["sstart"],
            "subject_end": row["send"],
            "evalue": row["evalue"],
            "bitscore": as_float(row["bitscore"]),
            "aligned_query": row["qseq"],
            "aligned_subject": row["sseq"],
        })
    return enriched


def best_hit_per_subject(rows):
    grouped = {}
    for row in rows:
        current = grouped.get(row["subject_id"])
        score = (row["bitscore"], row["alignment_length"], row["percent_identity"])
        if current is None or score > (current["bitscore"], current["alignment_length"], current["percent_identity"]):
            grouped[row["subject_id"]] = row
    return list(grouped.values())


def conserved_regions(segment, rows, subject_ids, query_length, min_identity, min_length):
    subject_ids = set(subject_ids)
    if not subject_ids:
        return []
    coverage = [set() for _ in range(query_length + 1)]
    for row in rows:
        if row["subject_id"] not in subject_ids:
            continue
        if row["percent_identity"] < min_identity or row["alignment_length"] < min_length:
            continue
        for pos in range(row["query_start"], row["query_end"] + 1):
            if 1 <= pos <= query_length:
                coverage[pos].add(row["subject_id"])

    regions = []
    start = None
    for pos in range(1, query_length + 2):
        complete = pos <= query_length and len(coverage[pos]) == len(subject_ids)
        if complete and start is None:
            start = pos
        if not complete and start is not None:
            end = pos - 1
            regions.append({
                "segment": segment,
                "query_start": start,
                "query_end": end,
                "length": end - start + 1,
                "gene_records_covering_region": len(subject_ids),
                "threshold": f">={min_identity}% identity and >={min_length} aligned bases",
            })
            start = None
    return regions


def rows_to_matrix(rows, fieldnames):
    return [fieldnames] + [[row.get(field, "") for field in fieldnames] for row in rows]


def main():
    parser = argparse.ArgumentParser(description="Align Andes virus L/M/S genes against gene FASTA records.")
    parser.add_argument("--gene-dir", type=Path, default=GENE_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-length", type=int, default=20)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    combined_all = []
    combined_best = []
    combined_regions = []
    summary = []

    for fasta in sorted(args.gene_dir.glob("*.fasta")):
        segment = segment_from_path(fasta)
        records = list(read_fasta(fasta))
        if not records:
            continue
        andes = [(title, seq) for title, seq in records if "andes virus" in title.lower()]
        if not andes:
            print(f"Skipping {fasta}: no Andes virus query record found.")
            continue
        query_title, query_seq = andes[0]
        segment_dir = args.output_dir / segment
        segment_dir.mkdir(parents=True, exist_ok=True)

        query_id = f"Andes_{segment}|{safe_id(query_title.split()[0])}"
        query_fasta = segment_dir / f"Andes_{segment}_query.fasta"
        db_fasta = segment_dir / f"{segment}_gene_records_db.fasta"
        raw_blast = segment_dir / f"Andes_{segment}_vs_{segment}_gene_records_blast.tsv"

        db_records = []
        id_to_title = {}
        for index, (title, seq) in enumerate(records, 1):
            record_id = f"{segment}|{safe_id(title.split()[0])}|record_{index}"
            id_to_title[record_id] = title
            db_records.append((record_id, title, seq))

        write_fasta([(query_id, query_title, query_seq)], query_fasta)
        write_fasta(db_records, db_fasta)
        run_blast(query_fasta, db_fasta, raw_blast)

        enriched = enrich_hits(segment, read_blast(raw_blast), id_to_title, len(query_seq))
        best = best_hit_per_subject(enriched)
        non_query_subjects = [
            record_id for record_id, title in id_to_title.items()
            if "andes virus" not in title.lower()
        ]
        regions = conserved_regions(segment, enriched, non_query_subjects, len(query_seq), args.min_identity, args.min_length)

        fieldnames = [
            "segment", "query_id", "subject_id", "species", "gene_title", "percent_identity",
            "alignment_length", "query_coverage_percent", "query_start", "query_end",
            "subject_start", "subject_end", "evalue", "bitscore", "aligned_query", "aligned_subject",
        ]
        region_fields = ["segment", "query_start", "query_end", "length", "gene_records_covering_region", "threshold"]
        write_tsv(segment_dir / f"Andes_{segment}_all_hits.tsv", enriched, fieldnames)
        write_tsv(segment_dir / f"Andes_{segment}_best_hits.tsv", best, fieldnames)
        write_tsv(segment_dir / f"Andes_{segment}_similar_in_all_gene_records.tsv", regions, region_fields)
        write_xlsx(
            segment_dir / f"Andes_{segment}_blast_alignment.xlsx",
            [
                ("Best hits", rows_to_matrix(best, fieldnames), [9, 22, 24, 24, 70, 16, 16, 22, 12, 12, 14, 14, 14, 12, 50, 50]),
                ("All hits", rows_to_matrix(enriched, fieldnames), [9, 22, 24, 24, 70, 16, 16, 22, 12, 12, 14, 14, 14, 12, 50, 50]),
                ("Similar in all genes", rows_to_matrix(regions, region_fields), [9, 14, 14, 12, 28, 42]),
            ],
        )

        combined_all.extend(enriched)
        combined_best.extend(best)
        combined_regions.extend(regions)
        summary.append({
            "segment": segment,
            "andes_query": query_title,
            "query_length": len(query_seq),
            "gene_records_in_db": len(records),
            "non_andes_gene_records": len(non_query_subjects),
            "blast_hits": len(enriched),
            "best_hit_records": len(best),
            "regions_similar_in_all_non_andes_records": len(regions),
        })

    fieldnames = [
        "segment", "query_id", "subject_id", "species", "gene_title", "percent_identity",
        "alignment_length", "query_coverage_percent", "query_start", "query_end",
        "subject_start", "subject_end", "evalue", "bitscore", "aligned_query", "aligned_subject",
    ]
    region_fields = ["segment", "query_start", "query_end", "length", "gene_records_covering_region", "threshold"]
    summary_fields = [
        "segment", "andes_query", "query_length", "gene_records_in_db", "non_andes_gene_records",
        "blast_hits", "best_hit_records", "regions_similar_in_all_non_andes_records",
    ]

    write_tsv(args.output_dir / "Andes_LMS_all_blast_hits.tsv", combined_all, fieldnames)
    write_tsv(args.output_dir / "Andes_LMS_best_hits.tsv", combined_best, fieldnames)
    write_tsv(args.output_dir / "Andes_LMS_regions_similar_in_all_gene_records.tsv", combined_regions, region_fields)
    write_xlsx(
        args.output_dir / "Andes_LMS_blast_alignment_summary.xlsx",
        [
            ("Summary", rows_to_matrix(summary, summary_fields), [9, 80, 14, 18, 22, 12, 16, 34]),
            ("Best hits", rows_to_matrix(combined_best, fieldnames), [9, 22, 24, 24, 70, 16, 16, 22, 12, 12, 14, 14, 14, 12, 50, 50]),
            ("Similar in all genes", rows_to_matrix(combined_regions, region_fields), [9, 14, 14, 12, 28, 42]),
        ],
    )

    print(f"Wrote alignment results to {args.output_dir}")


if __name__ == "__main__":
    main()
