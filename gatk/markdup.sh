#!/bin/sh
set -eu

SAM=""
THREADS=""
OUTPUT_DIR=""
NAME=""

echo "MARKDUP ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --wgs_align.sam|--wgs_align_sam)
            SAM="$2"
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

[ -n "$SAM" ]        || { echo "ERROR: missing --wgs_align.sam" >&2; exit 2; }
[ -n "$THREADS" ]    || { echo "ERROR: missing --threads" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]       || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

SORTED_BAM="$OUTPUT_DIR/$NAME.sorted.bam"
DEDUP_BAM="$OUTPUT_DIR/$NAME.dedup.bam"

# Primary alignments only (-F 0x900 excludes secondary + supplementary), then
# coordinate-sort. bwa-mem2's raw SAM has no sorting/filtering applied yet.
samtools view -@ "$THREADS" -F 0x900 -bS "$SAM" | \
    samtools sort -@ "$THREADS" -o "$SORTED_BAM"
samtools index -@ "$THREADS" "$SORTED_BAM"

gatk --java-options "-Xmx28g" MarkDuplicates \
    -I "$SORTED_BAM" \
    -O "$DEDUP_BAM" \
    -M "$OUTPUT_DIR/$NAME.dedup_metrics.txt"

samtools index -@ "$THREADS" "$DEDUP_BAM"
