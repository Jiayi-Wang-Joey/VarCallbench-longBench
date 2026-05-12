#!/usr/bin/env bash
set -euo pipefail

DATASET=""
BAM=""
REF=""
GTF=""
TX_DB=""
THREADS="64"
OUTDIR=""

echo "ARGS: $@" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --name|--dataset_id)
            DATASET="$2"; shift 2 ;;
        --isolaser_annotate.bam|--bam)
            BAM="$2"; shift 2 ;;
        --reference_genome)
            REF="$2"; shift 2 ;;
        --gtf)
            GTF="$2"; shift 2 ;;
        --transcriptome_db)
            TX_DB="$2"; shift 2 ;;
        --threads)
            THREADS="$2"; shift 2 ;;
        --output_dir)
            OUTDIR="$2"; shift 2 ;;
        --task)
            shift 2 ;;
        *)
            echo "Unknown arg: $1" >&2
            shift ;;
    esac
done

sample_lc="$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')"

if [[ "$sample_lc" == *"drna"* ]]; then
    PLATFORM="Nanopore"
elif [[ "$sample_lc" == *"ont"* ]]; then
    PLATFORM="Nanopore"
elif [[ "$sample_lc" == *"pb"* ]]; then
    PLATFORM="PacBio"
else
    echo "Unrecognized platform for sample '$DATASET'" >&2
    exit 1
fi

mkdir -p "$OUTDIR/calls"

isolaser \
    -b "$BAM" \
    -o "$OUTDIR/calls/${DATASET}" \
    -t "$TX_DB" \
    -f "$REF" \
    -n "$THREADS" \
    --DP=0 \
    --platform="$PLATFORM" \
    > "$OUTDIR/isolaser.log" 2>&1

# Convert gVCF to VCF: remove reference-only blocks, trim residual <NON_REF> alleles
bcftools view -e 'ALT="<NON_REF>"' "$OUTDIR/calls/${DATASET}.gvcf" \
    | bcftools norm --trim-alt-alleles -Ov \
    | bcftools sort \
        -T "$OUTDIR/tmp.sort" \
        -Oz \
        -o "$OUTDIR/${DATASET}.vcf.gz"

tabix -f -p vcf "$OUTDIR/${DATASET}.vcf.gz"
