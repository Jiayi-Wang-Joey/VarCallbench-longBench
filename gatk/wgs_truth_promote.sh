#!/bin/sh
set -eu

# Filters the WGS hard-filter output down to PASS-only records and stages it at a
# stable, non-hashed path (data/truth/<sample>.vcf.gz) so downstream RNA-seq-side
# omnibenchmark stages -- which live under a completely separate dataset root
# (rawdata/align/variant_call, not wgs_download) and therefore cannot reference this
# stage's output via normal omnibenchmark input wiring -- can pick it up as a static
# param, the same way the (now-retired) gnomad/somatic truth comparisons did.

VCF=""
TRUTH_DIR=""
OUTPUT_DIR=""
NAME=""

echo "WGS_TRUTH_PROMOTE ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --wgs_ground_truth.vcf.gz|--wgs_ground_truth_vcf_gz)
            VCF="$2"
            shift 2
            ;;
        --truth_dir)
            TRUTH_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --name)
            NAME="$2"
            shift 2
            ;;
        --task)
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

[ -n "$VCF" ]        || { echo "ERROR: missing --wgs_ground_truth.vcf.gz" >&2; exit 2; }
[ -n "$TRUTH_DIR" ]  || { echo "ERROR: missing --truth_dir" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]       || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR" "$TRUTH_DIR"

# Dataset names are H211_WGS / H526_WGS -- strip the _WGS suffix to get the sample
# prefix shared with the RNA-seq dataset names (H211_bulk_PB, H526_dRNA_ONT, ...).
SAMPLE=$(echo "$NAME" | sed -E 's/_WGS$//')

OUT="$OUTPUT_DIR/$NAME.vcf.gz"
bcftools view -f PASS -Oz -o "$OUT" "$VCF"
bcftools index -t -f "$OUT"

cp "$OUT" "$TRUTH_DIR/$SAMPLE.vcf.gz"
cp "$OUT.tbi" "$TRUTH_DIR/$SAMPLE.vcf.gz.tbi"
