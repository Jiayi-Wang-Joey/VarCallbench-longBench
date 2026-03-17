#!/usr/bin/env bash
set -euo pipefail

DATASET=""
BAM=""
REF=""
THREADS="10"
OUTDIR=""

echo "ARGS: $@" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --name|--dataset_id)
            DATASET="$2"; shift 2 ;;
        --align.bam|--bam)
            BAM="$2"; shift 2 ;;
        --reference_genome)
            REF="$2"; shift 2 ;;
        --threads)
            THREADS="$2"; shift 2 ;;
        --output_dir)
            OUTDIR="$2"; shift 2 ;;
        *)
            echo "Unknown arg: $1" >&2
            shift ;;
    esac
done

sample_lc="$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')"

if [[ "$sample_lc" == *"drna_ont"* ]]; then
    PLATFORM="ont-drna"
elif [[ "$sample_lc" == *"ont"* ]]; then
    PLATFORM="ont-cdna"
elif [[ "$sample_lc" == *"pb"* ]]; then
    PLATFORM="hifi-masseq"
else
    echo "Unrecognized platform for sample '$DATASET'" >&2
    exit 1
fi

if [ ! -f "${BAM}.bai" ] || [ "$BAM" -nt "${BAM}.bai" ]; then
    samtools index "$BAM"
fi

mkdir -p "$OUTDIR"

longcallR \
    --bam-path "$BAM" \
    --ref-path "$REF" \
    --output "$OUTDIR/output" \
    --preset "$PLATFORM" \
    -t "$THREADS" \
    --no-bam-output \
    --min-depth 0 \
    --max-depth 1000000 \
    > "$OUTDIR/longcallR.log" 2>&1

bgzip "$OUTDIR/output.vcf"
mv "$OUTDIR/output.vcf.gz" "$OUTDIR/output.vcf.gz.tmp"
mv "$OUTDIR/output.vcf.gz.tmp" "$OUTDIR/longcallR.vcf.gz"
tabix -p vcf "$OUTDIR/longcallR.vcf.gz"