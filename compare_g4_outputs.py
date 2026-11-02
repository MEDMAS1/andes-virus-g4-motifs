#!/usr/bin/env python3
"""Compare G4 predictions from all tool outputs.

The script normalizes available result files into calls with:
input, tool, start, end, sequence. It then writes reports for common
sequences, common start positions, overlaps, and sequences contained in
longer predictions from another tool.
"""

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
OUTPUT_DIR = RESULTS_DIR / "common"


@dataclass(frozen=True)
class Call:
    input_name: str
    tool: str
    start: int | None
    end: int | None
    sequence: str
    source_file: Path

    @property
    def norm_sequence(self):
        return normalize_sequence(self.sequence)


def normalize_sequence(sequence):
    return re.sub(r"[^A-Za-z]", "", sequence or "").upper().replace("U", "T")


def fasta_sequences(data_dir):
    sequences = {}
    for pattern in ("*.fasta", "*.fa", "*.fas"):
        for path in data_dir.glob(pattern):
            chunks = []
            with path.open() as handle:
                for line in handle:
                    line = line.strip()
                    if line and not line.startswith(">"):
                        chunks.append(line)
            sequences[path.stem] = normalize_sequence("".join(chunks))
    return sequences


def read_dict_rows(path, delimiter="\t"):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            yield {key.strip(): value.strip() for key, value in row.items() if key is not None}


def to_int(value):
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def call_from_fields(input_name, tool, source_file, start, end, sequence):
    sequence = normalize_sequence(sequence)
    if not sequence:
        return None
    if start is None:
        return None
    if end is None:
        end = start + len(sequence) - 1
    return Call(input_name, tool, start, end, sequence, source_file)


def reverse_complement(sequence):
    table = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return normalize_sequence(sequence).translate(table)[::-1]


def load_g4hunter(results_dir):
    calls = []
    for path in (results_dir / "G4Hunter").glob("*/*/*_g4hunter_website_format.tsv"):
        input_name = path.parents[1].name
        for row in read_dict_rows(path):
            calls.append(
                call_from_fields(
                    input_name,
                    "G4Hunter",
                    path,
                    to_int(row.get("start")),
                    to_int(row.get("end")),
                    row.get("sequence", ""),
                )
            )
    return [call for call in calls if call]


def load_g4catchall(results_dir):
    calls = []
    for path in (results_dir / "G4Catchall").glob("*/*_g4catchall.tsv"):
        input_name = path.parent.name
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 7:
                    continue
                start0 = to_int(parts[1])
                end0 = to_int(parts[2])
                if start0 is None or end0 is None:
                    continue
                sequence = parts[6] or parts[5]
                calls.append(call_from_fields(input_name, "G4Catchall", path, start0 + 1, end0, sequence))
    return [call for call in calls if call]


def load_g4imgrinder(results_dir):
    calls = []
    for path in (results_dir / "G4iMGrinder").glob("*/*_g4imgrinder_PQSM*.tsv"):
        input_name = path.parent.name
        for row in read_dict_rows(path):
            calls.append(
                call_from_fields(
                    input_name,
                    "G4iMGrinder",
                    path,
                    to_int(row.get("Start")),
                    to_int(row.get("Finish")),
                    row.get("Sequence", ""),
                )
            )
    return [call for call in calls if call]


def load_g4promfinder(results_dir, input_sequences):
    calls = []
    for path in (results_dir / "G4PromFinder").glob("*/*_g4promfinder_coordinates.txt"):
        input_name = path.parent.name
        full_sequence = input_sequences.get(input_name, "")
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("Region"):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                strand_label = parts[1].lower()
                start = to_int(parts[2])
                end = to_int(parts[3])
                if start is None or end is None or not full_sequence:
                    continue
                sequence = full_sequence[start - 1:end]
                if strand_label.startswith("neg"):
                    sequence = reverse_complement(sequence)
                calls.append(call_from_fields(input_name, "G4PromFinder", path, start, end, sequence))
    return [call for call in calls if call]


def load_g4predict(results_dir, input_sequences):
    calls = []
    for path in (results_dir / "G4Predict").glob("*/*_g4predict*.bed"):
        input_name = path.parent.name
        full_sequence = input_sequences.get(input_name, "")
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                start0 = to_int(parts[1])
                end0 = to_int(parts[2])
                if start0 is None or end0 is None or end0 <= start0:
                    continue
                sequence = full_sequence[start0:end0] if full_sequence else ""
                calls.append(call_from_fields(input_name, "G4Predict", path, start0 + 1, end0, sequence))
    return [call for call in calls if call]


def load_pqs_finder(results_dir):
    calls = []
    for path in (results_dir / "PQS_Finder").glob("*/*_pqs_finder.tsv"):
        input_name = path.parent.name
        for row in read_dict_rows(path):
            calls.append(
                call_from_fields(
                    input_name,
                    "PQS_Finder",
                    path,
                    to_int(row.get("start")),
                    to_int(row.get("end")),
                    row.get("sequence", ""),
                )
            )
    return [call for call in calls if call]


def load_qgrs(results_dir):
    calls = []
    for path in (results_dir / "QGRS").glob("*/*_qgrs_website_format.tsv"):
        if re.search(r"_website_format_\d{8}_\d{6}\.tsv$", path.name):
            continue
        input_name = path.parent.name
        for row in read_dict_rows(path):
            start = to_int(row.get("Position"))
            sequence = row.get("QGRS", "")
            length = to_int(row.get("Length"))
            end = start + length - 1 if start is not None and length is not None else None
            calls.append(call_from_fields(input_name, "QGRS", path, start, end, sequence))
    return [call for call in calls if call]


def load_qparse(results_dir):
    calls = []
    for path in (results_dir / "QPARSE").glob("*/*_qparse.txt"):
        input_name = path.parent.name
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith(("#", ">")):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                sequence = parts[0]
                start = to_int(parts[2])
                end = to_int(parts[3])
                calls.append(call_from_fields(input_name, "QPARSE", path, start, end, sequence))
    return [call for call in calls if call]


def load_quadron(results_dir):
    calls = []
    for path in (results_dir / "Quadron").glob("*/*_quadron.txt"):
        input_name = path.parent.name
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line.startswith("DATA:"):
                    continue
                parts = line.split(maxsplit=5)
                if len(parts) < 6 or parts[1].upper() == "NA":
                    continue
                start = to_int(parts[1])
                length = to_int(parts[3])
                sequence = parts[5]
                end = start + length - 1 if start is not None and length is not None else None
                calls.append(call_from_fields(input_name, "Quadron", path, start, end, sequence))
    return [call for call in calls if call]


def load_g4rna_screener(results_dir, min_g4nn):
    calls = []
    for path in (results_dir / "G4RNA_Screener").glob("*/*_g4rna_screener.tsv"):
        input_name = path.parent.name
        for row in read_dict_rows(path):
            if min_g4nn is not None:
                try:
                    if float(row.get("G4NN", "nan")) < min_g4nn:
                        continue
                except ValueError:
                    continue
            start = to_int(row.get("start"))
            sequence = row.get("sequence", "")
            end = start + len(normalize_sequence(sequence)) - 1 if start is not None else None
            calls.append(call_from_fields(input_name, "G4RNA_Screener", path, start, end, sequence))
    return [call for call in calls if call]


def load_g4mismatch(results_dir, input_sequences):
    calls = []
    for path in (results_dir / "G4mismatch").glob("*/*_g4mismatch.bedGraph"):
        input_name = path.parent.name
        full_sequence = input_sequences.get(input_name, "")
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                start0 = to_int(parts[1])
                end0 = to_int(parts[2])
                if start0 is None or end0 is None or end0 <= start0:
                    continue
                start = start0 + 1
                end = end0
                sequence = full_sequence[start0:end0] if full_sequence else ""
                calls.append(call_from_fields(input_name, "G4mismatch", path, start, end, sequence))
    return [call for call in calls if call]


def load_rg4detector_detection(results_dir, input_sequences):
    calls = []
    for path in (results_dir / "rG4detector").glob("*/rG4detector_detection.csv"):
        input_name = path.parent.name
        full_sequence = input_sequences.get(input_name, "")
        for row in read_dict_rows(path, delimiter=","):
            lower = {key.lower(): value for key, value in row.items()}
            start = to_int(lower.get("start"))
            end = to_int(lower.get("end"))
            sequence = lower.get("sequence", "")
            if not sequence and start is not None and end is not None and full_sequence:
                sequence = full_sequence[start - 1:end]
            calls.append(call_from_fields(input_name, "rG4detector", path, start, end, sequence))
    return [call for call in calls if call]


def load_calls(results_dir, data_dir, min_g4nn):
    input_sequences = fasta_sequences(data_dir)
    calls = []
    calls.extend(load_g4catchall(results_dir))
    calls.extend(load_g4hunter(results_dir))
    calls.extend(load_g4imgrinder(results_dir))
    calls.extend(load_g4predict(results_dir, input_sequences))
    calls.extend(load_g4promfinder(results_dir, input_sequences))
    calls.extend(load_pqs_finder(results_dir))
    calls.extend(load_qgrs(results_dir))
    calls.extend(load_qparse(results_dir))
    calls.extend(load_quadron(results_dir))
    calls.extend(load_g4rna_screener(results_dir, min_g4nn))
    calls.extend(load_g4mismatch(results_dir, input_sequences))
    calls.extend(load_rg4detector_detection(results_dir, input_sequences))
    return sorted(set(calls), key=lambda c: (c.input_name, c.start or 0, c.end or 0, c.tool, c.sequence))


def tools(calls):
    return sorted({call.tool for call in calls})


def source_names(calls, root):
    return ";".join(sorted(str(call.source_file.relative_to(root)) for call in calls))


def writable_output(path):
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a"):
            pass
        return path
    except PermissionError:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return path.with_name(f"{path.stem}_{stamp}{path.suffix}")


def write_normalized_calls(calls, output, root):
    output = writable_output(output)
    with output.open("w", newline="") as handle:
        fieldnames = ["input", "tool", "start", "end", "length", "sequence", "source_file"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for call in calls:
            writer.writerow(
                {
                    "input": call.input_name,
                    "tool": call.tool,
                    "start": call.start,
                    "end": call.end,
                    "length": len(call.norm_sequence),
                    "sequence": call.norm_sequence,
                    "source_file": call.source_file.relative_to(root),
                }
            )
    return output


def write_common_sequences(calls, output, min_tools, root):
    groups = defaultdict(list)
    for call in calls:
        groups[(call.input_name, call.norm_sequence)].append(call)

    output = writable_output(output)
    with output.open("w", newline="") as handle:
        fieldnames = [
            "input",
            "sequence",
            "length",
            "tool_count",
            "tools",
            "tool_starts",
            "starts",
            "ends",
            "source_files",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for (input_name, sequence), grouped in sorted(groups.items()):
            grouped_tools = tools(grouped)
            if len(grouped_tools) < min_tools:
                continue
            tool_starts = ";".join(f"{call.tool}:{call.start}" for call in grouped)
            writer.writerow(
                {
                    "input": input_name,
                    "sequence": sequence,
                    "length": len(sequence),
                    "tool_count": len(grouped_tools),
                    "tools": ";".join(grouped_tools),
                    "tool_starts": tool_starts,
                    "starts": ";".join(str(call.start) for call in grouped),
                    "ends": ";".join(str(call.end) for call in grouped),
                    "source_files": source_names(grouped, root),
                }
            )
    return output


def write_common_starts(calls, output, min_tools, root):
    groups = defaultdict(list)
    for call in calls:
        groups[(call.input_name, call.start)].append(call)

    output = writable_output(output)
    with output.open("w", newline="") as handle:
        fieldnames = ["input", "start", "tool_count", "tools", "ends", "sequences", "source_files"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for (input_name, start), grouped in sorted(groups.items()):
            grouped_tools = tools(grouped)
            if len(grouped_tools) < min_tools:
                continue
            writer.writerow(
                {
                    "input": input_name,
                    "start": start,
                    "tool_count": len(grouped_tools),
                    "tools": ";".join(grouped_tools),
                    "ends": ";".join(str(call.end) for call in grouped),
                    "sequences": ";".join(call.norm_sequence for call in grouped),
                    "source_files": source_names(grouped, root),
                }
            )
    return output


def write_overlaps(calls, output, min_overlap):
    by_input = defaultdict(list)
    for call in calls:
        by_input[call.input_name].append(call)

    output = writable_output(output)
    with output.open("w", newline="") as handle:
        fieldnames = [
            "input",
            "tool_a",
            "start_a",
            "end_a",
            "sequence_a",
            "tool_b",
            "start_b",
            "end_b",
            "sequence_b",
            "overlap_start",
            "overlap_end",
            "overlap_length",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for input_name, grouped in sorted(by_input.items()):
            for index, left in enumerate(grouped):
                for right in grouped[index + 1:]:
                    if left.tool == right.tool:
                        continue
                    overlap_start = max(left.start, right.start)
                    overlap_end = min(left.end, right.end)
                    overlap_length = overlap_end - overlap_start + 1
                    if overlap_length < min_overlap:
                        continue
                    writer.writerow(
                        {
                            "input": input_name,
                            "tool_a": left.tool,
                            "start_a": left.start,
                            "end_a": left.end,
                            "sequence_a": left.norm_sequence,
                            "tool_b": right.tool,
                            "start_b": right.start,
                            "end_b": right.end,
                            "sequence_b": right.norm_sequence,
                            "overlap_start": overlap_start,
                            "overlap_end": overlap_end,
                            "overlap_length": overlap_length,
                        }
                    )
    return output


def write_contained_sequences(calls, output):
    by_input = defaultdict(list)
    for call in calls:
        by_input[call.input_name].append(call)

    output = writable_output(output)
    with output.open("w", newline="") as handle:
        fieldnames = [
            "input",
            "contained_tool",
            "contained_start",
            "contained_end",
            "contained_sequence",
            "container_tool",
            "container_start",
            "container_end",
            "container_sequence",
            "included_at_position",
            "included_offset",
            "match_type",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for input_name, grouped in sorted(by_input.items()):
            for small in grouped:
                for large in grouped:
                    if small == large or small.tool == large.tool:
                        continue
                    if len(small.norm_sequence) > len(large.norm_sequence):
                        continue
                    offset = large.norm_sequence.find(small.norm_sequence)
                    interval_contained = large.start <= small.start and small.end <= large.end
                    if offset < 0 and not interval_contained:
                        continue
                    included_at = large.start + offset if offset >= 0 else small.start
                    if small == large:
                        continue
                    writer.writerow(
                        {
                            "input": input_name,
                            "contained_tool": small.tool,
                            "contained_start": small.start,
                            "contained_end": small.end,
                            "contained_sequence": small.norm_sequence,
                            "container_tool": large.tool,
                            "container_start": large.start,
                            "container_end": large.end,
                            "container_sequence": large.norm_sequence,
                            "included_at_position": included_at,
                            "included_offset": offset + 1 if offset >= 0 else "",
                            "match_type": "sequence" if offset >= 0 else "position",
                        }
                    )
    return output


def main():
    parser = argparse.ArgumentParser(description="Find common G4 predictions across tool outputs.")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--min-tools", type=int, default=2, help="Minimum tools required for common reports.")
    parser.add_argument("--min-overlap", type=int, default=1, help="Minimum shared bases for overlap report.")
    parser.add_argument(
        "--min-g4nn",
        type=float,
        default=None,
        help="Optional minimum G4NN score for G4RNA Screener windows.",
    )
    args = parser.parse_args()

    calls = load_calls(args.results_dir, args.data_dir, args.min_g4nn)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    written = [
        write_normalized_calls(calls, args.output_dir / "all_predictions_normalized.tsv", ROOT),
        write_common_sequences(calls, args.output_dir / "common_sequences.tsv", args.min_tools, ROOT),
        write_common_starts(calls, args.output_dir / "common_start_positions.tsv", args.min_tools, ROOT),
        write_overlaps(calls, args.output_dir / "overlapping_predictions.tsv", args.min_overlap),
        write_contained_sequences(calls, args.output_dir / "contained_sequences.tsv"),
    ]

    seen_tools = ", ".join(sorted({call.tool for call in calls})) or "none"
    print(f"Loaded {len(calls)} predictions from: {seen_tools}")
    print(f"Wrote reports to {args.output_dir}")
    for path in written:
        if path.parent == args.output_dir and re.search(r"_\d{8}_\d{6}\.", path.name):
            print(f"Locked output detected; wrote {path.name} instead.")
    if not any((args.results_dir / "rG4detector").glob("*/rG4detector_detection.csv")):
        print("Note: rG4detector prediction CSVs do not include start/end/sequence calls, so they were not compared.")


if __name__ == "__main__":
    main()
