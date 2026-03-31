#!/usr/bin/env bash
set -euo pipefail

BAM=""
THREADS="1"

echo "ARGS: $*" >&2

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

STDERR_PATH="$(readlink -f /proc/self/fd/2)"
OUTDIR="$(dirname "$STDERR_PATH")"
JOBDIR="$(pwd)"

echo "STDERR_PATH: $STDERR_PATH" >&2
echo "OUTDIR: $OUTDIR" >&2
echo "JOBDIR: $JOBDIR" >&2

# infer dataset from current path: out/rawdata/<dataset>/...
if [ -z "$BAM" ]; then
    DATASET="$(printf '%s\n' "$JOBDIR" | sed -n 's#.*out/rawdata/\([^/]*\)/.*#\1#p')"
    if [ -z "${DATASET:-}" ]; then
        echo "Could not infer dataset from working directory" >&2
        exit 1
    fi

    BAM_CANDIDATE="$(readlink -f "$JOBDIR/../../../${DATASET}.aligned.bam" || true)"
    if [ -n "${BAM_CANDIDATE:-}" ] && [ -f "$BAM_CANDIDATE" ]; then
        BAM="$BAM_CANDIDATE"
    fi
fi

if [ -z "$BAM" ]; then
    echo "Could not determine aligned BAM" >&2
    exit 1
fi

echo "Using BAM: $BAM" >&2

DATASET="$(basename "$BAM")"
DATASET="${DATASET%.aligned.bam}"
DATASET="${DATASET%.bam}"

STATS="$OUTDIR/${DATASET}.samtools.stats"
CSV="$OUTDIR/${DATASET}.alignment_qc.csv"

samtools stats -@ "$THREADS" "$BAM" > "$STATS"

awk -F '\t' -v bam="$BAM" -v dataset="$DATASET" 'BEGIN {
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
ls -lh "$OUTDIR" >&2