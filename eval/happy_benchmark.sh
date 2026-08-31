#!/bin/sh
set -eu

CALLER_VCF=""
CALLER_VCF_SOURCE=""
BED=""
TRUTH_DIR=""
REFERENCE=""
THREADS=""
OUTPUT_DIR=""
NAME=""

echo "HAPPY_BENCHMARK ARGS: $*" >&2

while [ "$#" -gt 0 ]; do
    case "$1" in
        --variant.vcf|--variant_vcf)
            CALLER_VCF="$2"
            CALLER_VCF_SOURCE="variant_call"
            shift 2
            ;;
        --isolaser.vcf|--isolaser_vcf|--isolaser_filtered.vcf|--isolaser_filtered_vcf)
            CALLER_VCF="$2"
            CALLER_VCF_SOURCE="isolaser_call"
            shift 2
            ;;
        --generate_bed_variant_call.bed|--generate_bed_variant_call_bed|--generate_bed_isolaser.bed|--generate_bed_isolaser_bed)
            BED="$2"
            shift 2
            ;;
        --truth_dir)
            TRUTH_DIR="$2"
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

[ -n "$CALLER_VCF" ] || { echo "ERROR: missing --variant.vcf or --isolaser.vcf" >&2; exit 2; }
[ -n "$BED" ]         || { echo "ERROR: missing --generate_bed.bed" >&2; exit 2; }
[ -n "$TRUTH_DIR" ]   || { echo "ERROR: missing --truth_dir" >&2; exit 2; }
[ -n "$REFERENCE" ]   || { echo "ERROR: missing --reference_genome" >&2; exit 2; }
[ -n "$THREADS" ]     || { echo "ERROR: missing --threads" >&2; exit 2; }
[ -n "$OUTPUT_DIR" ]  || { echo "ERROR: missing --output_dir" >&2; exit 2; }
[ -n "$NAME" ]        || { echo "ERROR: missing --name" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

# Tool name: isoLASER is its own single-module stage; the other 4 callers are all
# modules of the variant_call stage, distinguished by path (same convention already
# used by the retired somatic_detection.sh / gnomad_detection.sh scripts).
if [ "$CALLER_VCF_SOURCE" = "isolaser_call" ]; then
    TOOL="isoLASER"
else
    TOOL=$(echo "$CALLER_VCF" | sed -En 's#.*/variant_call/([^/]+)/.*#\1#p')
    [ -n "$TOOL" ] || TOOL="unknown"
fi

# Sample prefix (H211/H526) shared with the WGS ground-truth dataset names, derived
# from the RNA-seq dataset name (H211_bulk_PB, H526_dRNA_ONT, ...).
SAMPLE=$(echo "$NAME" | sed -E 's/_(bulk|dRNA)_(PB|ONT)$//')
TRUTH_VCF="$TRUTH_DIR/$SAMPLE.vcf.gz"
[ -f "$TRUTH_VCF" ] || { echo "ERROR: truth VCF not found: $TRUTH_VCF" >&2; exit 2; }

# generate_bed is invoked once per --min_coverage value (a benchmark.yaml parameter
# list); the value itself isn't visible in its (hash-based) output path, so read it
# back from the sidecar file generate_bed.sh drops next to the BED.
COVERAGE_FILE="$(dirname "$BED")/$NAME.min_coverage.txt"
[ -f "$COVERAGE_FILE" ] || { echo "ERROR: coverage sidecar file not found: $COVERAGE_FILE" >&2; exit 2; }
COVERAGE=$(cat "$COVERAGE_FILE")

SRC="$OUTPUT_DIR/$NAME.src.vcf.gz"
cp "$CALLER_VCF" "$SRC"

if ! bcftools index -t -f "$SRC" >/dev/null 2>&1; then
    echo "Caller VCF is not sorted/indexed -- sorting..." >&2
    mv "$SRC" "$SRC.unsorted"
    bcftools sort -Oz -o "$SRC" "$SRC.unsorted"
    bcftools index -t -f "$SRC"
    rm -f "$SRC.unsorted"
fi

QUERY="$OUTPUT_DIR/$NAME.query.vcf.gz"
bcftools view -f PASS -Oz -o "$QUERY" "$SRC"
bcftools index -t -f "$QUERY"
rm -f "$SRC" "$SRC.tbi"

HAPPY_PREFIX="$OUTPUT_DIR/$NAME"

# No RTG Tools/vcfeval available in this container -- use hap.py's built-in xcmp
# engine. No -f/confidence-region flag: this project has no confident-region BED
# (unlike GIAB); -T alone restricts the whole comparison to the benchmarking region.
HGREF="$REFERENCE" hap.py \
    "$TRUTH_VCF" \
    "$QUERY" \
    -r "$REFERENCE" \
    -T "$BED" \
    -o "$HAPPY_PREFIX" \
    --threads "$THREADS"

rm -f "$QUERY" "$QUERY.tbi"

SUMMARY="$HAPPY_PREFIX.summary.csv"
FINAL="$OUTPUT_DIR/$NAME.happy_summary.csv"

awk -v tool="$TOOL" \
    -v aligner="minimap2" \
    -v sample="$NAME" \
    -v coverage="$COVERAGE" \
    'BEGIN {FS=OFS=","}
     NR==1 {print $0, "tool", "aligner", "sample", "coverage"}
     NR>1  {print $0, tool, aligner, sample, coverage}' \
    "$SUMMARY" > "$FINAL"
