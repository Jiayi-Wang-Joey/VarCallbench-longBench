#!/bin/sh
set -eu

BAM=""
OUTPUT_DIR=""
NAME=""

echo "BAM_INDEX ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --align.bam|--align_bam)
            BAM="$2"
            shift 2
            ;;
        # Chained after variant_call/isolaser_call purely to give this stage's
        # own output path a per-tool ancestor (omnibenchmark resolves a stage's
        # inputs via string-prefix ancestry on node ids, so two stages that both
        # just depend on align.bam in parallel can't both feed a later stage) --
        # the actual VCF content isn't used here.
        --variant.vcf|--variant_vcf|--isolaser.vcf|--isolaser_vcf|--isolaser_filtered.vcf|--isolaser_filtered_vcf)
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

[ -n "$BAM" ]        || { echo "ERROR: missing --align.bam" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]       || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

# mosdepth (in its own, minimal container downstream) needs the BAM and its
# index sitting side by side; symlink the (large, already-computed) BAM into
# this stage's output dir rather than copying it, then index the symlink so
# the .bai lands next to it.
ln -sf "$BAM" "$OUTPUT_DIR/$NAME.bam"
samtools index "$OUTPUT_DIR/$NAME.bam"
