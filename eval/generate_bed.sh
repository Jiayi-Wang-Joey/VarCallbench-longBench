#!/bin/sh
set -eu

COVERAGE_GZ=""
GTF=""
MIN_COVERAGE=""
OUTPUT_DIR=""
NAME=""

echo "GENERATE_BED ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --generate_bed_variant_call_mosdepth.coverage.bed.gz|--generate_bed_variant_call_mosdepth_coverage_bed_gz|--generate_bed_isolaser_mosdepth.coverage.bed.gz|--generate_bed_isolaser_mosdepth_coverage_bed_gz)
            COVERAGE_GZ="$2"
            shift 2
            ;;
        --gtf)
            GTF="$2"
            shift 2
            ;;
        --min_coverage)
            MIN_COVERAGE="$2"
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

[ -n "$COVERAGE_GZ" ]   || { echo "ERROR: missing --*_mosdepth.coverage.bed.gz" >&2; exit 2; }
[ -n "$GTF" ]           || { echo "ERROR: missing --gtf" >&2; exit 2; }
[ -n "$MIN_COVERAGE" ]  || { echo "ERROR: missing --min_coverage" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ]    || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]          || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

COVERAGE_BED="$OUTPUT_DIR/$NAME.coverage.bed"

zcat "$COVERAGE_GZ" \
    | awk -v min_cov="$MIN_COVERAGE" '$4 >= min_cov' \
    | bedtools merge -d 1 -c 4 -o mean -i - \
    > "$COVERAGE_BED"

# No confidence/high-confidence region is available for this project (unlike GIAB) --
# the benchmarking region is simply coverage-filtered RNA-seq alignment intersected
# with GENCODE exons.
grep -v "^#" "$GTF" | awk '$3 == "exon"' \
    | bedtools intersect -a "$COVERAGE_BED" -b - \
    > "$OUTPUT_DIR/$NAME.bed"

rm -f "$COVERAGE_BED"

echo "$MIN_COVERAGE" > "$OUTPUT_DIR/$NAME.min_coverage.txt"
