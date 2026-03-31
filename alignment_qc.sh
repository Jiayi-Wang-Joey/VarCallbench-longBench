#!/usr/bin/env bash
set -euo pipefail

BAM=""
OUTDIR=""
THREADS="1"
TASK=""

echo "ARGS: $*" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --task)
            TASK="$2"
            shift 2
            ;;
        --align.bam|--align_bam|--align-bam|--bam)
            BAM="$2"
            shift 2
            ;;
        --output_dir|--output-dir|--output.dir)
            OUTDIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [ -z "$BAM" ]; then
    echo "Missing --align.bam/--align_bam/--align-bam/--bam" >&2
    exit 1
fi

if [ -z "$OUTDIR" ]; then
    echo "Missing --output_dir" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

STATS="$OUTDIR/samtools.stats"
CSV="$OUTDIR/alignment_qc.csv"

samtools stats -@ "$THREADS" "$BAM" > "$STATS"

awk -F '\t' -v bam="$BAM" 'BEGIN {
    OFS=","
    print "dataset_id,bam,total_sequences,mapped_reads,mapping_rate_pct,average_length,average_quality,error_rate"
}
$1 == "SN" {
    gsub(":$", "", $2)
    val[$2] = $3
}
END {
    total = val["raw total sequences"]
    mapped = val["reads mapped"]
    avg_len = val["average length"]
    avg_qual = val["average quality"]
    err = val["error rate"]

    rate = 0
    if (total > 0) {
        rate = 100 * mapped / total
    }

    dataset = bam
    sub(/^.*\//, "", dataset)
    sub(/\.aligned\.bam$/, "", dataset)
    sub(/\.bam$/, "", dataset)

    print dataset, bam, total, mapped, rate, avg_len, avg_qual, err
}' "$STATS" > "$CSV"

echo "Wrote: $CSV" >&2