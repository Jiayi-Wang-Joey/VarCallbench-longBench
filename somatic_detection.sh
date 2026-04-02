#!/usr/bin/env bash
set -euo pipefail

BAM=""
THREADS="1"

echo "ARGS: $*" >&2

OUTDIR=""

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
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --output_dir)
            OUTDIR="$2"
            shift 2
            ;;
        --name|--dataset_id)
            DATASET="$2"
            shift 2
            ;;
        --variant.vcf|--variant_vcf|--variant-vcf)
            VCF="$2"
            shift 2
            ;;
        --somatic_vcf)
            SOMATIC_VCF="$2"
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

echo "STDERR_PATH: $STDERR_PATH" >&2
echo "OUTDIR: $OUTDIR" >&2

if [ -z "$BAM" ]; then
    DATASET="$(printf '%s\n' "$STDERR_PATH" | sed -n 's#.*out/rawdata/\([^/]*\)/.*#\1#p')"
    if [ -z "${DATASET:-}" ]; then
        echo "Could not infer dataset from stderr path" >&2
        exit 1
    fi

    ALIGN_DIR="$(printf '%s\n' "$STDERR_PATH" | sed 's#/alignment_qc/.*##')"
    BAM="${ALIGN_DIR}/${DATASET}.aligned.bam"
fi

if [ ! -f "$BAM" ]; then
    echo "Could not find BAM: $BAM" >&2
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