#!/usr/bin/env bash
set -euo pipefail

DATASET=""
BAM=""
TRANSCRIPTOME_FA=""
GTF=""
THREADS="10"
OUTDIR=""

echo "ARGS: $@" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --name|--dataset_id)
            DATASET="$2"; shift 2 ;;
        --align.bam|--bam)
            BAM="$2"; shift 2 ;;
        --transcriptome_fa)
            TRANSCRIPTOME_FA="$2"; shift 2 ;;
        --gtf)
            GTF="$2"; shift 2 ;;
        --threads)
            THREADS="$2"; shift 2 ;;
        --output_dir)
            OUTDIR="$2"; shift 2 ;;
        --task)
            shift 2 ;;
        *)
            echo "Unknown arg: $1" >&2
            shift ;;
    esac
done

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/tmp"

if [ ! -f "${BAM}.bai" ] || [ "$BAM" -nt "${BAM}.bai" ]; then
    samtools index "$BAM"
fi

SAM="$OUTDIR/tmp/${DATASET}.aligned.sam"
UNSORTED="$OUTDIR/tmp/${DATASET}.anno.unsorted.bam"
OUT_BAM="$OUTDIR/${DATASET}.anno.bam"

samtools bam2fq "$BAM" \
    | minimap2 -t "$THREADS" -ax splice:hq -uf --MD "$TRANSCRIPTOME_FA" - \
    > "$SAM" 2>> "$OUTDIR/isolaser_annotate.log"

isolaser_annotate \
    -b "$BAM" \
    -t "$SAM" \
    -g "$GTF" \
    -o "$UNSORTED" 2>> "$OUTDIR/isolaser_annotate.log"

samtools sort -@ "$THREADS" -m 4G \
    -T "$OUTDIR/tmp/${DATASET}_sort" \
    -o "$OUT_BAM" "$UNSORTED" 2>> "$OUTDIR/isolaser_annotate.log"

samtools index "$OUT_BAM"

rm -f "$SAM" "$UNSORTED"
