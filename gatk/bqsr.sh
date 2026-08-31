#!/bin/sh
set -eu

BAM=""
REFERENCE=""
KNOWN_SITES_DBSNP=""
KNOWN_SITES_MILLS=""
OUTPUT_DIR=""
NAME=""

echo "BQSR ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --wgs_markdup.bam|--wgs_markdup_bam)
            BAM="$2"
            shift 2
            ;;
        --reference_genome)
            REFERENCE="$2"
            shift 2
            ;;
        --known_sites_dbsnp)
            KNOWN_SITES_DBSNP="$2"
            shift 2
            ;;
        --known_sites_mills)
            KNOWN_SITES_MILLS="$2"
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

[ -n "$BAM" ]               || { echo "ERROR: missing --wgs_markdup.bam" >&2; exit 2; }
[ -n "$REFERENCE" ]         || { echo "ERROR: missing --reference_genome" >&2; exit 2; }
[ -n "$KNOWN_SITES_DBSNP" ] || { echo "ERROR: missing --known_sites_dbsnp" >&2; exit 2; }
[ -n "$KNOWN_SITES_MILLS" ] || { echo "ERROR: missing --known_sites_mills" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ]        || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]              || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

RECAL_TABLE="$OUTPUT_DIR/$NAME.recal.table"
RECAL_BAM="$OUTPUT_DIR/$NAME.recal.bam"

gatk --java-options "-Xmx28g" BaseRecalibrator \
    -I "$BAM" \
    -R "$REFERENCE" \
    --known-sites "$KNOWN_SITES_DBSNP" \
    --known-sites "$KNOWN_SITES_MILLS" \
    -O "$RECAL_TABLE"

gatk --java-options "-Xmx28g" ApplyBQSR \
    -I "$BAM" \
    -R "$REFERENCE" \
    --bqsr-recal-file "$RECAL_TABLE" \
    -O "$RECAL_BAM"

samtools index "$RECAL_BAM"
