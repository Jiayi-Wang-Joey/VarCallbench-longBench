#!/usr/bin/env bash
set -euo pipefail

BAM=""
THREADS="1"

echo "ARGS: $*" >&2
echo "PWD: $(pwd)" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --task)
            shift 2
            ;;
        --align.bam|--align_bam|--align-bam|--bam)
            BAM="$2"
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
    BAM="$(find "$(pwd)/../.." -type f -name "*.aligned.bam" | head -n 1 || true)"
fi

if [ -z "$BAM" ]; then
    echo "Could not find an aligned BAM automatically" >&2
    exit 1
fi

echo "Using BAM: $BAM" >&2

dataset=$(basename "$BAM")
dataset=${dataset%.aligned.bam}
dataset=${dataset%.bam}

STATS="$(pwd)/${dataset}.samtools.stats"
CSV="$(pwd)/${dataset}.alignment_qc.csv"

samtools stats -@ "$THREADS" "$BAM" > "$STATS"

awk -F '\t' -v bam="$BAM" -v dataset="$dataset" 'BEGIN {
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

    print dataset, bam, total, mapped, rate, avg_len, avg_qual, err
}' "$STATS" > "$CSV"

echo "Wrote: $CSV" >&2
ls -lh "$(pwd)" >&2