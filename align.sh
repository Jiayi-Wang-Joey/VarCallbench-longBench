#!/bin/sh
set -eu

READS=""
REFERENCE_GENOME=""
TRANSCRIPTOME_BED=""
ALIGN_MAP_BAM_THREADS=""
ALIGN_SORT_BAM_THREADS=""
ALIGN_SORT_BAM_MEMORY_GB=""
OUTPUT_DIR=""
NAME=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --rawdata_fastq)
            READS="$2"
            shift 2
            ;;
        --reference_genome)
            REFERENCE_GENOME="$2"
            shift 2
            ;;
        --transcriptome_bed)
            TRANSCRIPTOME_BED="$2"
            shift 2
            ;;
        --align_map_bam_threads)
            ALIGN_MAP_BAM_THREADS="$2"
            shift 2
            ;;
        --align_sort_bam_threads)
            ALIGN_SORT_BAM_THREADS="$2"
            shift 2
            ;;
        --align_sort_bam_memory_gb)
            ALIGN_SORT_BAM_MEMORY_GB="$2"
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
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [ -z "$READS" ] || [ -z "$REFERENCE_GENOME" ] || [ -z "$TRANSCRIPTOME_BED" ] || \
   [ -z "$ALIGN_MAP_BAM_THREADS" ] || [ -z "$ALIGN_SORT_BAM_THREADS" ] || \
   [ -z "$ALIGN_SORT_BAM_MEMORY_GB" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$NAME" ]; then
    echo "Missing required arguments" >&2
    exit 2
fi

case "$NAME" in
    *IsoSeq*|*MasSeq*)
        PRESET="-ax splice:hq -uf"
        ;;
    *dRNA*)
        PRESET="-ax splice -uf -k14"
        ;;
    *cDNA*)
        PRESET="-ax splice"
        ;;
    *)
        echo "Cannot determine minimap2 preset from dataset name: $NAME" >&2
        exit 2
        ;;
esac

mkdir -p "$OUTPUT_DIR"

minimap2 $PRESET --junc-bed "$TRANSCRIPTOME_BED" \
    -t "$ALIGN_MAP_BAM_THREADS" \
    "$REFERENCE_GENOME" "$READS" | \
samtools sort \
    -@ "$ALIGN_SORT_BAM_THREADS" \
    -m "${ALIGN_SORT_BAM_MEMORY_GB}g" \
    -o "$OUTPUT_DIR/$NAME.aligned.bam"