#!/bin/sh
set -eu

BAM=""
THREADS=""
OUTPUT_DIR=""
NAME=""

echo "MOSDEPTH_COVERAGE ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --generate_bed_variant_call_bam_index.bam|--generate_bed_variant_call_bam_index_bam|--generate_bed_isolaser_bam_index.bam|--generate_bed_isolaser_bam_index_bam)
            BAM="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
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

[ -n "$BAM" ]        || { echo "ERROR: missing --*_bam_index.bam" >&2; exit 2; }
[ -n "$THREADS" ]    || { echo "ERROR: missing --threads" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]       || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

mosdepth --threads "$THREADS" "$OUTPUT_DIR/$NAME" "$BAM"
