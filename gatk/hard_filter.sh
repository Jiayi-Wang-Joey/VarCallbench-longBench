#!/bin/sh
set -eu

VCF=""
REFERENCE=""
OUTPUT_DIR=""
NAME=""

echo "HARD_FILTER ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --wgs_genotypegvcfs.vcf.gz|--wgs_genotypegvcfs_vcf_gz)
            VCF="$2"
            shift 2
            ;;
        --reference_genome)
            REFERENCE="$2"
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

[ -n "$VCF" ]        || { echo "ERROR: missing --wgs_genotypegvcfs.vcf.gz" >&2; exit 2; }
[ -n "$REFERENCE" ]  || { echo "ERROR: missing --reference_genome" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]       || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

SNPS="$OUTPUT_DIR/$NAME.snps.vcf.gz"
INDELS="$OUTPUT_DIR/$NAME.indels.vcf.gz"
SNPS_FILTERED="$OUTPUT_DIR/$NAME.snps.filtered.vcf.gz"
INDELS_FILTERED="$OUTPUT_DIR/$NAME.indels.filtered.vcf.gz"

gatk --java-options "-Xmx8g" SelectVariants -R "$REFERENCE" -V "$VCF" --select-type-to-include SNP -O "$SNPS"
gatk --java-options "-Xmx8g" SelectVariants -R "$REFERENCE" -V "$VCF" --select-type-to-include INDEL -O "$INDELS"

# GATK's recommended germline hard-filter expressions (used in place of VQSR, which
# needs a large multi-sample training cohort we don't have here).
gatk --java-options "-Xmx8g" VariantFiltration -R "$REFERENCE" -V "$SNPS" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    -O "$SNPS_FILTERED"

gatk --java-options "-Xmx8g" VariantFiltration -R "$REFERENCE" -V "$INDELS" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10" \
    -O "$INDELS_FILTERED"

gatk --java-options "-Xmx8g" MergeVcfs \
    -I "$SNPS_FILTERED" \
    -I "$INDELS_FILTERED" \
    -O "$OUTPUT_DIR/$NAME.vcf.gz"
