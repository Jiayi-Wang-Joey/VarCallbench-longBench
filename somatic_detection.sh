#!/usr/bin/env bash
set -euo pipefail

TASK=""
OUTDIR=""
DATASET=""
VCF=""
SOMATIC_VCF=""
SOMATIC_DIR=""
EXON_BED=""
THREADS="1"

echo "ARGS: $*" >&2

########################################
# parse arguments
########################################
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
        --exon_bed)
            EXON_BED="$2"
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

########################################
# checks
########################################
[ -n "$OUTDIR" ] || { echo "Missing --output_dir" >&2; exit 1; }
[ -n "$DATASET" ] || { echo "Missing --name/--dataset_id" >&2; exit 1; }
[ -n "$VCF" ] || { echo "Missing --variant.vcf" >&2; exit 1; }
[ -n "$EXON_BED" ] || { echo "Missing --exon_bed" >&2; exit 1; }

[ -f "$VCF" ] || { echo "VCF not found: $VCF" >&2; exit 1; }
[ -f "$EXON_BED" ] || { echo "Exon BED not found: $EXON_BED" >&2; exit 1; }

if [ -z "$SOMATIC_VCF" ]; then
    [ -n "$SOMATIC_DIR" ] || { echo "Missing --somatic_dir" >&2; exit 1; }
    BASE_DATASET="${DATASET%%_*}"
    SOMATIC_VCF="${SOMATIC_DIR}/${BASE_DATASET}.vcf.gz"
fi

[ -f "$SOMATIC_VCF" ] || { echo "Somatic VCF not found: $SOMATIC_VCF" >&2; exit 1; }

command -v bcftools >/dev/null 2>&1 || { echo "bcftools not found" >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found" >&2; exit 1; }

########################################
# extract caller
########################################
CALLER=$(echo "$VCF" | sed -n 's#.*/variant_call/\([^/]*\)/.*#\1#p')

case "$CALLER" in
    clair3_rna) CALLER="Clair3-RNA" ;;
    deep_variant) CALLER="DeepVariant" ;;
    longcallR) CALLER="longcallR" ;;
    longcallR_nn) CALLER="longcallR-nn" ;;
    *) CALLER="unknown" ;;
esac

########################################
# prepare
########################################
mkdir -p "$OUTDIR"

echo "TASK: ${TASK:-somatic_detection}" >&2
echo "DATASET: $DATASET" >&2
echo "CALLER: $CALLER" >&2
echo "VCF: $VCF" >&2
echo "SOMATIC_VCF: $SOMATIC_VCF" >&2
echo "EXON_BED: $EXON_BED" >&2

########################################
# somatic detection
########################################
TMPDIR="$OUTDIR/isec_tmp"
rm -rf "$TMPDIR"
mkdir -p "$TMPDIR"

tabix -f "$VCF"

bcftools isec -p "$TMPDIR" "$VCF" "$SOMATIC_VCF"

DETECTED=0
MISSED=0
UNMATCHED_CALLS=0

[ -f "$TMPDIR/0002.vcf" ] && DETECTED=$(bcftools view -H "$TMPDIR/0002.vcf" | wc -l | tr -d ' ')
[ -f "$TMPDIR/0001.vcf" ] && MISSED=$(bcftools view -H "$TMPDIR/0001.vcf" | wc -l | tr -d ' ')
[ -f "$TMPDIR/0000.vcf" ] && UNMATCHED_CALLS=$(bcftools view -H "$TMPDIR/0000.vcf" | wc -l | tr -d ' ')

TOTAL_TRUTH=$((DETECTED + MISSED))

if [ "$TOTAL_TRUTH" -gt 0 ]; then
    DETECTION_RATE=$(awk -v d="$DETECTED" -v t="$TOTAL_TRUTH" 'BEGIN{printf "%.6f", d/t}')
else
    DETECTION_RATE="NA"
fi

########################################
# exon / intron classification
########################################
VAR_BED="$OUTDIR/variants.bed"

bcftools view -H "$VCF" | \
awk 'BEGIN{OFS="\t"} {
    start=$2-1;
    end=start+length($4);
    if (end <= start) end=start+1;
    print $1, start, end
}' > "$VAR_BED"

bedtools intersect -a "$VAR_BED" -b "$EXON_BED" -u > "$OUTDIR/exonic.bed"
bedtools intersect -a "$VAR_BED" -b "$EXON_BED" -v > "$OUTDIR/intronic.bed"

EXONIC_CALLS=$(wc -l < "$OUTDIR/exonic.bed" | tr -d ' ')
INTRONIC_CALLS=$(wc -l < "$OUTDIR/intronic.bed" | tr -d ' ')

########################################
# output
########################################
OUTFILE="$OUTDIR/${DATASET}.somatic_detection.csv"

{
    echo "dataset,caller,detected,total_truth,missed,detection_rate,unmatched_calls,exonic_calls,intronic_calls"
    echo "${DATASET},${CALLER},${DETECTED},${TOTAL_TRUTH},${MISSED},${DETECTION_RATE},${UNMATCHED_CALLS},${EXONIC_CALLS},${INTRONIC_CALLS}"
} > "$OUTFILE"

########################################
# logs
########################################
echo "DETECTED: $DETECTED" >&2
echo "MISSED: $MISSED" >&2
echo "TOTAL_TRUTH: $TOTAL_TRUTH" >&2
echo "DETECTION_RATE: $DETECTION_RATE" >&2
echo "UNMATCHED_CALLS: $UNMATCHED_CALLS" >&2
echo "EXONIC_CALLS: $EXONIC_CALLS" >&2
echo "INTRONIC_CALLS: $INTRONIC_CALLS" >&2
echo "Wrote: $OUTFILE" >&2