#!/bin/sh
set -eu

READ1=""
READ2=""
REFERENCE=""
THREADS=""
OUTPUT_DIR=""
NAME=""

echo "WGS ALIGN ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --wgs_R1.fastq|--wgs_R1_fastq)
            READ1="$2"
            shift 2
            ;;
        --wgs_R2.fastq|--wgs_R2_fastq)
            READ2="$2"
            shift 2
            ;;
        --reference_genome)
            REFERENCE="$2"
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

[ -n "$REFERENCE" ]  || { echo "ERROR: missing --reference_genome" >&2; exit 2; }
[ -n "$THREADS" ]    || { echo "ERROR: missing --threads" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]       || { echo "ERROR: missing --name" >&2; exit 2; }
[ -n "$READ1" ]      || { echo "ERROR: missing reads input (--wgs_R1.fastq)" >&2; exit 2; }
[ -n "$READ2" ]      || { echo "ERROR: missing reads input (--wgs_R2.fastq)" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

# bwa-mem2 container has no samtools, so this stage only produces raw SAM. Sorting,
# primary-alignment filtering (-F 0x900) and indexing happen at the start of the
# markdup stage (gatk_env container, which bundles samtools).
bwa-mem2 mem \
    -t "$THREADS" \
    -R "@RG\tID:$NAME\tSM:$NAME\tPL:ILLUMINA\tLB:$NAME" \
    "$REFERENCE" "$READ1" "$READ2" \
    > "$OUTPUT_DIR/$NAME.aligned.sam"
