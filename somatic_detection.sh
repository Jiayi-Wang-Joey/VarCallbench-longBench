#!/usr/bin/env bash
set -euo pipefail

TASK=""
OUTDIR=""
DATASET=""
BAM=""
VCF=""
SOMATIC_VCF=""
REF=""
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
        --align.bam|--align_bam|--align-bam|--bam)
            BAM="$2"
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
        --reference_genome)
            REF="$2"
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

if [ -z "$TASK" ]; then
    echo "Missing required argument: --task" >&2
    exit 1
fi

if [ -z "$OUTDIR" ]; then
    echo "Missing required argument: --output_dir" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

echo "TASK: $TASK" >&2
echo "DATASET: ${DATASET:-NA}" >&2
echo "OUTDIR: $OUTDIR" >&2
echo "THREADS: $THREADS" >&2

case "$TASK" in
    align)
        if [ -z "$BAM" ]; then
            echo "Missing required BAM for task '$TASK': --align.bam" >&2
            exit 1
        fi

        echo "BAM: $BAM" >&2
        echo "REF: ${REF:-NA}" >&2

        exec bash "$(dirname "$0")/align.sh" \
            --name "$DATASET" \
            --bam "$BAM" \
            --reference_genome "$REF" \
            --threads "$THREADS" \
            --output_dir "$OUTDIR"
        ;;

    clair3_rna)
        if [ -z "$BAM" ]; then
            echo "Missing required BAM for task '$TASK': --align.bam" >&2
            exit 1
        fi

        echo "BAM: $BAM" >&2
        echo "REF: ${REF:-NA}" >&2

        exec bash "$(dirname "$0")/clair3_rna.sh" \
            --name "$DATASET" \
            --bam "$BAM" \
            --reference_genome "$REF" \
            --threads "$THREADS" \
            --output_dir "$OUTDIR"
        ;;

    deep_variant|deepvariant)
        if [ -z "$BAM" ]; then
            echo "Missing required BAM for task '$TASK': --align.bam" >&2
            exit 1
        fi

        echo "BAM: $BAM" >&2
        echo "REF: ${REF:-NA}" >&2

        exec bash "$(dirname "$0")/deep_variant.sh" \
            --name "$DATASET" \
            --bam "$BAM" \
            --reference_genome "$REF" \
            --threads "$THREADS" \
            --output_dir "$OUTDIR"
        ;;

    longcallR|longcallr)
        if [ -z "$BAM" ]; then
            echo "Missing required BAM for task '$TASK': --align.bam" >&2
            exit 1
        fi

        echo "BAM: $BAM" >&2
        echo "REF: ${REF:-NA}" >&2

        exec bash "$(dirname "$0")/longcallR.sh" \
            --name "$DATASET" \
            --bam "$BAM" \
            --reference_genome "$REF" \
            --threads "$THREADS" \
            --output_dir "$OUTDIR"
        ;;

    somatic_detection)
        if [ -z "$VCF" ]; then
            echo "Missing required VCF for task '$TASK': --variant.vcf" >&2
            exit 1
        fi

        if [ -z "$SOMATIC_VCF" ]; then
            echo "Missing required somatic truth VCF for task '$TASK': --somatic_vcf" >&2
            exit 1
        fi

        echo "VCF: $VCF" >&2
        echo "SOMATIC_VCF: $SOMATIC_VCF" >&2

        exec Rscript "$(dirname "$0")/somatic_detection.R" \
            --output_dir "$OUTDIR" \
            --variant.vcf "$VCF" \
            --somatic_vcf "$SOMATIC_VCF" \
            --name "$DATASET"
        ;;

    *)
        echo "ERROR: unknown task: $TASK" >&2
        exit 1
        ;;
esac