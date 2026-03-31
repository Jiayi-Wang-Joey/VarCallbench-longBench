#!/usr/bin/env bash
set -euo pipefail

BAM=""
THREADS="1"
OUTDIR=""

while [ $# -gt 0 ]; do
    case "$1" in
        --align.bam|--bam)
            BAM="$2"; shift 2 ;;
        --threads)
            THREADS="$2"; shift 2 ;;
        --output_dir)
            OUTDIR="$2"; shift 2 ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1 ;;
    esac
done

if [ -z "$BAM" ] || [ -z "$OUTDIR" ]; then
    echo "Missing --align.bam/--bam or --output_dir" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

samtools stats -@ "$THREADS" "$BAM" > "$OUTDIR/samtools.stats"

python - "$BAM" "$OUTDIR/samtools.stats" "$OUTDIR/alignment_qc.csv" <<'PY'
import csv
import os
import sys

bam = sys.argv[1]
stats_file = sys.argv[2]
out_csv = sys.argv[3]

metrics = {
    "dataset_id": os.path.basename(bam).replace(".aligned.sorted.bam", "").replace(".bam", ""),
    "bam": bam,
    "total_sequences": "",
    "mapped_reads": "",
    "mapping_rate_pct": "",
    "average_length": "",
    "average_quality": "",
    "error_rate": "",
    "insert_size_average": ""
}

wanted = {
    "raw total sequences:": "total_sequences",
    "reads mapped:": "mapped_reads",
    "reads mapped and paired:": None,
    "average length:": "average_length",
    "average quality:": "average_quality",
    "error rate:": "error_rate",
    "insert size average:": "insert_size_average"
}

with open(stats_file) as fh:
    for line in fh:
        if not line.startswith("SN\t"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        key = parts[1]
        val = parts[2]
        if key in wanted and wanted[key] is not None:
            metrics[wanted[key]] = val

total = metrics["total_sequences"]
mapped = metrics["mapped_reads"]
try:
    total_f = float(total)
    mapped_f = float(mapped)
    metrics["mapping_rate_pct"] = 100 * mapped_f / total_f if total_f > 0 else ""
except Exception:
    metrics["mapping_rate_pct"] = ""

fieldnames = [
    "dataset_id",
    "bam",
    "total_sequences",
    "mapped_reads",
    "mapping_rate_pct",
    "average_length",
    "average_quality",
    "error_rate",
    "insert_size_average"
]

with open(out_csv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow(metrics)
PY