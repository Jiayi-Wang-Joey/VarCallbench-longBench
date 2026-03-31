#!/usr/bin/env bash
set -euo pipefail

OUT="gnomad_detection.csv"
THREADS=4
GTF=""
GNOMAD_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.vcf.bgz"

VCFS=()

while [ $# -gt 0 ]; do
    case "$1" in
        --variant.vcf|--variant_vcf|--variant-vcf)
            shift
            while [ $# -gt 0 ] && [[ "$1" != --* ]]; do
                VCFS+=("$1")
                shift
            done
            ;;
        --gtf)
            GTF="$2"; shift 2 ;;
        --output)
            OUT="$2"; shift 2 ;;
        --threads)
            THREADS="$2"; shift 2 ;;
        --gnomad_url)
            GNOMAD_URL="$2"; shift 2 ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [ ${#VCFS[@]} -eq 0 ]; then
    echo "ERROR: no VCF files provided" >&2
    exit 1
fi

if [ -z "$GTF" ]; then
    echo "ERROR: missing --gtf" >&2
    exit 1
fi

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

EXON_BED="$TMPDIR/exons.merged.bed"

# GTF -> merged exon BED (0-based, half-open BED format)
awk '
BEGIN { OFS="\t" }
$0 ~ /^#/ { next }
$3 == "exon" {
    start = $4 - 1
    end = $5
    if (start < 0) start = 0
    print $1, start, end
}
' "$GTF" \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> "$EXON_BED"

echo "dataset,caller,total_pass_exonic_variants,gnomad_overlap,detection_rate" > "$OUT"

for VCF in "${VCFS[@]}"; do
    echo "Processing $VCF" >&2

    base=$(basename "$VCF")
    dataset="${base%.vcf.gz}"
    dataset="${dataset%.vcf}"

    caller="unknown"
    case "$VCF" in
        *"/variant_call/clair3_rna/"*)   caller="Clair3-RNA" ;;
        *"/variant_call/deep_variant/"*) caller="DeepVariant" ;;
        *"/variant_call/longcallR/"*)    caller="longcallR" ;;
        *"/variant_call/longcallR_nn/"*) caller="longcallR-nn" ;;
    esac

    PASS_VCF="$TMPDIR/${dataset}.pass.vcf.gz"
    EXONIC_VCF="$TMPDIR/${dataset}.pass.exonic.vcf.gz"
    POS_BED="$TMPDIR/${dataset}.sites.bed"
    GNOMAD_SUB="$TMPDIR/${dataset}.gnomad.subset.vcf.gz"

    # 1) keep PASS only
    bcftools view -f PASS -Oz -o "$PASS_VCF" "$VCF"
    tabix -f "$PASS_VCF"

    # 2) keep only variants overlapping exon BED from GTF
    bcftools view -R "$EXON_BED" -Oz -o "$EXONIC_VCF" "$PASS_VCF"
    tabix -f "$EXONIC_VCF"

    # 3) total number of PASS exonic variants
    TOTAL=$(bcftools view -H "$EXONIC_VCF" | wc -l | awk '{print $1}')

    if [ "$TOTAL" -eq 0 ]; then
        echo "$dataset,$caller,0,0,0" >> "$OUT"
        continue
    fi

    # 4) get unique variant positions from filtered VCF
    bcftools query -f '%CHROM\t%POS0\t%END\n' "$EXONIC_VCF" \
        | sort -u > "$POS_BED"

    # 5) query gnomAD only at these positions
    bcftools view -R "$POS_BED" "$GNOMAD_URL" -Oz -o "$GNOMAD_SUB"
    tabix -f "$GNOMAD_SUB"

    # 6) overlap count (exact CHROM,POS,REF,ALT match)
    OVERLAP=$(bcftools isec -n=2 -w1 "$EXONIC_VCF" "$GNOMAD_SUB" \
        | bcftools view -H \
        | wc -l \
        | awk '{print $1}')

    RATE=$(awk -v o="$OVERLAP" -v t="$TOTAL" 'BEGIN { printf "%.6f", (t > 0 ? o/t : 0) }')

    echo "$dataset,$caller,$TOTAL,$OVERLAP,$RATE" >> "$OUT"
done

echo "Wrote: $OUT"