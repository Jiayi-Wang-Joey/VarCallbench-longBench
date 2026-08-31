#!/bin/sh
set -eu

BAM=""
REFERENCE=""
OUTPUT_DIR=""
NAME=""

echo "HAPLOTYPECALLER ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --wgs_bqsr.bam|--wgs_bqsr_bam)
            BAM="$2"
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

[ -n "$BAM" ]        || { echo "ERROR: missing --wgs_bqsr.bam" >&2; exit 2; }
[ -n "$REFERENCE" ]  || { echo "ERROR: missing --reference_genome" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]       || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

gatk --java-options "-Xmx28g" HaplotypeCaller \
    -R "$REFERENCE" \
    -I "$BAM" \
    -O "$OUTPUT_DIR/$NAME.g.vcf.gz" \
    -ERC GVCF
