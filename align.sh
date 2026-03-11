#!/bin/sh
set -eu

READS=""
REFERENCE=""
JUNC_BED=""
MAP_THREADS=""
SORT_THREADS=""
SORT_MEM=""
OUTPUT_DIR=""
NAME=""

echo "ALIGN ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --rawdata.fastq|--rawdata_fastq)
            READS="$2"
            shift 2
            ;;
        --reference_genome)
            REFERENCE="$2"
            shift 2
            ;;
        --transcriptome_bed)
            JUNC_BED="$2"
            shift 2
            ;;
        --align_map_bam_threads)
            MAP_THREADS="$2"
            shift 2
            ;;
        --align_sort_bam_threads)
            SORT_THREADS="$2"
            shift 2
            ;;
        --align_sort_bam_memory_gb)
            SORT_MEM="$2"
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

[ -n "$READS" ] || { echo "ERROR: missing reads input" >&2; exit 2; }
[ -n "$REFERENCE" ] || { echo "ERROR: missing --reference_genome" >&2; exit 2; }
[ -n "$JUNC_BED" ] || { echo "ERROR: missing --transcriptome_bed" >&2; exit 2; }
[ -n "$MAP_THREADS" ] || { echo "ERROR: missing --align_map_bam_threads" >&2; exit 2; }
[ -n "$SORT_THREADS" ] || { echo "ERROR: missing --align_sort_bam_threads" >&2; exit 2; }
[ -n "$SORT_MEM" ] || { echo "ERROR: missing --align_sort_bam_memory_gb" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ] || { echo "ERROR: missing --name" >&2; exit 2; }

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
        echo "ERROR: Cannot determine minimap2 preset from dataset name: $NAME" >&2
        exit 2
        ;;
esac

mkdir -p "$OUTPUT_DIR"

minimap2 $PRESET --junc-bed "$JUNC_BED" \
    -t "$MAP_THREADS" \
    "$REFERENCE" "$READS" | \
samtools sort \
    -@ "$SORT_THREADS" \
    -m "${SORT_MEM}g" \
    -o "$OUTPUT_DIR/$NAME.aligned.bam"