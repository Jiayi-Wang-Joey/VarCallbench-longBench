#!/usr/bin/env bash
set -euo pipefail

DATASET=""
VCF=""
OUTDIR=""
ANNOTATION_BED=""

echo "ARGS: $*" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --task)
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
        --output_dir)
            OUTDIR="$2"
            shift 2
            ;;
        --annotation_bed)
            ANNOTATION_BED="$2"
            shift 2
            ;;
        --min_dp)
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [ -z "$DATASET" ]; then
    echo "Missing --name/--dataset_id" >&2
    exit 1
fi

if [ -z "$VCF" ]; then
    echo "Missing --variant.vcf" >&2
    exit 1
fi

if [ -z "$OUTDIR" ]; then
    echo "Missing --output_dir" >&2
    exit 1
fi

if [ -z "$ANNOTATION_BED" ]; then
    echo "Missing --annotation_bed" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

OUTVCF="${OUTDIR}/${DATASET}.vcf.gz"

echo "Dataset: $DATASET" >&2
echo "Input VCF: $VCF" >&2
echo "Annotation BED: $ANNOTATION_BED" >&2
echo "Output VCF: $OUTVCF" >&2

bcftools view \
    -f PASS \
    -R "$ANNOTATION_BED" \
    -Oz \
    -o "$OUTVCF" \
    "$VCF"

bcftools index -t "$OUTVCF"

echo "Done: $OUTVCF" >&2