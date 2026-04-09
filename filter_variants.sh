#!/usr/bin/env bash
set -euo pipefail

DATASET=""
VCF=""
OUTDIR=""
ANNOTATION_BED=""
MIN_DP="5"

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
            MIN_DP="$2"
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

OUTVCF="${OUTDIR}/${DATASET}.filtered.vcf.gz"

echo "Input VCF: $VCF" >&2
echo "Annotation BED: $ANNOTATION_BED" >&2
echo "Minimum DP: $MIN_DP" >&2
echo "Output VCF: $OUTVCF" >&2

# First try INFO/DP; if DP is not defined there, fall back to FORMAT/DP.
if bcftools view -h "$VCF" | grep -q 'ID=DP,Number=1,Type=Integer'; then
    DP_EXPR="INFO/DP>=${MIN_DP}"
else
    DP_EXPR="FORMAT/DP>=${MIN_DP}"
fi

bcftools view \
    -f PASS \
    -i "$DP_EXPR" \
    -R "$ANNOTATION_BED" \
    -Oz \
    -o "$OUTVCF" \
    "$VCF"

bcftools index -t "$OUTVCF"

echo "Done: $OUTVCF" >&2