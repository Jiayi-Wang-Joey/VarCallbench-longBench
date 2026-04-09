#!/usr/bin/env bash
set -euo pipefail

TASK=""
OUTDIR=""
DATASET=""
VCF=""
SOMATIC_VCF=""
THREADS="1"

echo "ARGS: $*" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --task)
            TASK="$2"
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
        --somatic_dir)
            SOMATIC_DIR="$2"
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

if [ -z "$OUTDIR" ]; then
    echo "Missing required argument: --output_dir" >&2
    exit 1
fi

if [ -z "$DATASET" ]; then
    echo "Missing required argument: --name/--dataset_id" >&2
    exit 1
fi

if [ -z "$VCF" ]; then
    echo "Missing required argument: --variant.vcf" >&2
    exit 1
fi

if [ -z "$SOMATIC_VCF" ]; then
    [ -n "$SOMATIC_DIR" ] || { echo "Missing --somatic_dir or --somatic_vcf" >&2; exit 1; }

    BASE_DATASET="${DATASET%%_*}"
    SOMATIC_VCF="${SOMATIC_DIR}/${BASE_DATASET}.vcf.gz"
fi

if [ ! -f "$VCF" ]; then
    echo "Variant VCF not found: $VCF" >&2
    exit 1
fi

if [ ! -f "$SOMATIC_VCF" ]; then
    echo "Somatic truth VCF not found: $SOMATIC_VCF" >&2
    exit 1
fi

command -v bcftools >/dev/null 2>&1 || {
    echo "bcftools not found in PATH" >&2
    exit 1
}

mkdir -p "$OUTDIR"

echo "TASK: ${TASK:-somatic_detection}" >&2
echo "DATASET: $DATASET" >&2
echo "OUTDIR: $OUTDIR" >&2
echo "THREADS: $THREADS" >&2
echo "VCF: $VCF" >&2
echo "SOMATIC_VCF: $SOMATIC_VCF" >&2

TMPDIR="$OUTDIR/isec_tmp"
rm -rf "$TMPDIR"
mkdir -p "$TMPDIR"
tabix "$VCF"
bcftools isec \
    -p "$TMPDIR" \
    "$VCF" \
    "$SOMATIC_VCF"

DETECTED=0
MISSED=0
TOTAL_TRUTH=0
UNMATCHED_CALLS=0

if [ -f "$TMPDIR/0002.vcf" ]; then
    DETECTED=$(bcftools view -H "$TMPDIR/0002.vcf" | wc -l | tr -d ' ')
fi

if [ -f "$TMPDIR/0001.vcf" ]; then
    MISSED=$(bcftools view -H "$TMPDIR/0001.vcf" | wc -l | tr -d ' ')
fi

if [ -f "$TMPDIR/0000.vcf" ]; then
    UNMATCHED_CALLS=$(bcftools view -H "$TMPDIR/0000.vcf" | wc -l | tr -d ' ')
fi

TOTAL_TRUTH=$((DETECTED + MISSED))

if [ "$TOTAL_TRUTH" -gt 0 ]; then
    DETECTION_RATE=$(awk -v d="$DETECTED" -v t="$TOTAL_TRUTH" 'BEGIN{printf "%.6f", d/t}')
else
    DETECTION_RATE="NA"
fi

OUTFILE="$OUTDIR/${DATASET}.somatic_detection.csv"

{
    echo "dataset,detected,total_truth,missed,detection_rate,unmatched_calls"
    echo "${DATASET},${DETECTED},${TOTAL_TRUTH},${MISSED},${DETECTION_RATE},${UNMATCHED_CALLS}"
} > "$OUTFILE"

echo "DETECTED: $DETECTED" >&2
echo "MISSED: $MISSED" >&2
echo "TOTAL_TRUTH: $TOTAL_TRUTH" >&2
echo "DETECTION_RATE: $DETECTION_RATE" >&2
echo "UNMATCHED_CALLS: $UNMATCHED_CALLS" >&2
echo "Wrote: $OUTFILE" >&2