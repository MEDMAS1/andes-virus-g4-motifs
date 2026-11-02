#!/usr/bin/env python3
"""Compare QGRS/PQS_Finder common sequences to best non-Andes genes."""

import argparse
import csv
import html
import re
import subprocess
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parent
ALIGNMENTS_DIR = ROOT / "alignments"
COMMON_DIR = ROOT / "results" / "common"
OUTPUT_DIR = ALIGNMENTS_DIR / "common_sequences_vs_best_genes"

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
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", str(value).strip()).strip("_") or "NA"


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


def latest_common_sequences(common_dir):
    files = sorted(common_dir.glob("common_sequences*.tsv"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not files:
        raise SystemExit(f"No common_sequences*.tsv found in {common_dir}")
    return files[0]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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


def get_sequence_by_id(fasta_path, wanted_id):
    for title, sequence in read_fasta(fasta_path):
        record_id = title.split()[0]
        if record_id == wanted_id:
            return title, sequence
    raise SystemExit(f"Could not find {wanted_id} in {fasta_path}")


def choose_best_non_andes(best_hits):
    best_by_segment = {}
    for row in best_hits:
        title = row.get("gene_title", "")
        if "andes virus" in title.lower():
            continue
        segment = row["segment"]
        current = best_by_segment.get(segment)
        score = (as_float(row["bitscore"]), as_float(row["alignment_length"]), as_float(row["percent_identity"]))
        if current is None:
            best_by_segment[segment] = row
            continue
        current_score = (
            as_float(current["bitscore"]),
            as_float(current["alignment_length"]),
            as_float(current["percent_identity"]),
        )
        if score > current_score:
            best_by_segment[segment] = row
    return best_by_segment


def write_fasta(records, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="\n") as handle:
        for record_id, title, sequence in records:
            handle.write(f">{record_id} {title}\n{sequence}\n")


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


def best_hit_per_query(blast_rows):
    best = {}
    for row in blast_rows:
        current = best.get(row["qseqid"])
        score = (as_float(row["bitscore"]), as_float(row["length"]), as_float(row["pident"]))
        if current is None:
            best[row["qseqid"]] = row
            continue
        current_score = (as_float(current["bitscore"]), as_float(current["length"]), as_float(current["pident"]))
        if score > current_score:
            best[row["qseqid"]] = row
    return best


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
    with zipfile.ZipFile(path) as zf:
        ET.fromstring(zf.read("xl/workbook.xml"))
        for name in zf.namelist():
            if name.startswith("xl/worksheets/sheet"):
                ET.fromstring(zf.read(name))


def rows_to_matrix(rows, fieldnames):
    return [fieldnames] + [[row.get(field, "") for field in fieldnames] for row in rows]


def main():
    parser = argparse.ArgumentParser(description="Compare common QGRS/PQS_Finder sequences with best Andes-aligned genes.")
    parser.add_argument("--common", type=Path, default=None)
    parser.add_argument("--alignments-dir", type=Path, default=ALIGNMENTS_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    common_file = args.common or latest_common_sequences(COMMON_DIR)
    best_hits_file = args.alignments_dir / "Andes_LMS_best_hits.tsv"
    best_by_segment = choose_best_non_andes(read_tsv(best_hits_file))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected_records = []
    selected_summary = []
    selected_sequence_by_segment = {}
    for segment, row in sorted(best_by_segment.items()):
        db_fasta = args.alignments_dir / segment / f"{segment}_gene_records_db.fasta"
        title, sequence = get_sequence_by_id(db_fasta, row["subject_id"])
        selected_sequence_by_segment[segment] = sequence
        selected_records.append((row["subject_id"], title, sequence))
        selected_summary.append({
            "segment": segment,
            "selected_best_gene": row["subject_id"],
            "species": row["species"],
            "gene_title": row["gene_title"],
            "andes_alignment_identity_percent": row["percent_identity"],
            "andes_alignment_length": row["alignment_length"],
            "andes_query_coverage_percent": row["query_coverage_percent"],
            "bitscore": row["bitscore"],
        })

    common_rows = []
    query_records = []
    for index, row in enumerate(read_tsv(common_file), 1):
        tools = row.get("tools", "")
        if "QGRS" not in tools or "PQS_Finder" not in tools:
            continue
        segment = row.get("input", "").strip()
        sequence = clean_sequence(row.get("sequence", ""))
        if segment not in best_by_segment or not sequence:
            continue
        common_id = f"common_{index}|segment_{safe_id(segment)}|start_{safe_id(row.get('starts', '').split(';')[0])}|len_{len(sequence)}"
        common_rows.append({
            "common_id": common_id,
            "segment": segment,
            "sequence": sequence,
            "length_nt": len(sequence),
            "andes_start": row.get("starts", ""),
            "andes_end": row.get("ends", ""),
            "tools": tools,
            "selected_best_gene": best_by_segment[segment]["subject_id"],
            "selected_best_gene_title": best_by_segment[segment]["gene_title"],
            "selected_best_gene_species": best_by_segment[segment]["species"],
        })
        query_records.append((common_id, f"{segment} common QGRS/PQS_Finder sequence", sequence))

    query_fasta = args.output_dir / "common_qgrs_pqsfinder_sequences.fasta"
    db_fasta = args.output_dir / "best_andes_aligned_genes.fasta"
    raw_blast = args.output_dir / "common_sequences_vs_best_genes_blast_raw.tsv"
    write_fasta(query_records, query_fasta)
    write_fasta(selected_records, db_fasta)
    run_blast(query_fasta, db_fasta, raw_blast)

    blast_best = best_hit_per_query(read_blast(raw_blast))
    comparison_rows = []
    for row in common_rows:
        hit = blast_best.get(row["common_id"])
        selected_sequence = selected_sequence_by_segment[row["segment"]]
        exact_pos = selected_sequence.find(row["sequence"])
        if hit:
            alignment_length = as_int(hit["length"])
            identity = as_float(hit["pident"])
            coverage = round(alignment_length / row["length_nt"] * 100, 2) if row["length_nt"] else 0.0
            subject_start = hit["sstart"]
            subject_end = hit["send"]
            evalue = hit["evalue"]
            bitscore = as_float(hit["bitscore"])
            aligned_query = hit["qseq"]
            aligned_subject = hit["sseq"]
        else:
            alignment_length = 0
            identity = 0.0
            coverage = 0.0
            subject_start = ""
            subject_end = ""
            evalue = ""
            bitscore = 0.0
            aligned_query = ""
            aligned_subject = ""

        comparison_rows.append({
            **row,
            "best_gene_exact_match": "yes" if exact_pos >= 0 else "no",
            "best_gene_exact_start": exact_pos + 1 if exact_pos >= 0 else "",
            "best_gene_exact_end": exact_pos + row["length_nt"] if exact_pos >= 0 else "",
            "blast_identity_percent": identity,
            "blast_query_coverage_percent": coverage,
            "blast_alignment_length": alignment_length,
            "best_gene_blast_start": subject_start,
            "best_gene_blast_end": subject_end,
            "evalue": evalue,
            "bitscore": bitscore,
            "aligned_common_sequence": aligned_query,
            "aligned_best_gene_sequence": aligned_subject,
        })

    selected_fields = [
        "segment", "selected_best_gene", "species", "gene_title",
        "andes_alignment_identity_percent", "andes_alignment_length",
        "andes_query_coverage_percent", "bitscore",
    ]
    comparison_fields = [
        "common_id", "segment", "sequence", "length_nt", "andes_start", "andes_end", "tools",
        "selected_best_gene", "selected_best_gene_species", "selected_best_gene_title",
        "best_gene_exact_match", "best_gene_exact_start", "best_gene_exact_end",
        "blast_identity_percent", "blast_query_coverage_percent", "blast_alignment_length",
        "best_gene_blast_start", "best_gene_blast_end", "evalue", "bitscore",
        "aligned_common_sequence", "aligned_best_gene_sequence",
    ]
    write_tsv(args.output_dir / "selected_best_andes_aligned_genes.tsv", selected_summary, selected_fields)
    write_tsv(args.output_dir / "common_qgrs_pqsfinder_vs_best_genes.tsv", comparison_rows, comparison_fields)
    write_xlsx(
        args.output_dir / "common_qgrs_pqsfinder_vs_best_genes.xlsx",
        [
            ("Comparison", rows_to_matrix(comparison_rows, comparison_fields), [30, 9, 48, 10, 16, 16, 18, 24, 24, 70, 16, 18, 18, 18, 22, 18, 18, 18, 14, 12, 45, 45]),
            ("Selected genes", rows_to_matrix(selected_summary, selected_fields), [9, 24, 24, 80, 24, 22, 24, 12]),
        ],
    )

    print(f"Read common QGRS/PQS_Finder sequences from {common_file}")
    print(f"Selected {len(selected_summary)} best non-Andes genes")
    print(f"Compared {len(comparison_rows)} common sequences")
    print(f"Wrote reports to {args.output_dir}")


if __name__ == "__main__":
    main()
