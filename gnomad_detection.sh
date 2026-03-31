#!/usr/bin/env bash
set -euo pipefail

GTF=""
GNOMAD_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.vcf.bgz"
THREADS="1"

echo "ARGS: $*" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --task)
            shift 2
            ;;
        --gtf)
            GTF="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --gnomad_url|--gnomad-url)
            GNOMAD_URL="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

STDERR_PATH="$(readlink -f /proc/self/fd/2)"
OUTDIR="$(dirname "$STDERR_PATH")"
PARAM_JSON="$OUTDIR/parameters.json"

echo "STDERR_PATH: $STDERR_PATH" >&2
echo "OUTDIR: $OUTDIR" >&2
echo "PARAM_JSON: $PARAM_JSON" >&2

if [ -z "$GTF" ] && [ -f "$PARAM_JSON" ]; then
    GTF="$(python - <<'PY' "$PARAM_JSON"
import json, sys
p = json.load(open(sys.argv[1]))
params = p.get("parameters", {})
gtf = params.get("gtf", "")
print(gtf)
PY
)"
fi

if [ -z "$GTF" ]; then
    echo "Missing gtf parameter" >&2
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "Could not find GTF: $GTF" >&2
    exit 1
fi

DATASET="$(printf '%s\n' "$STDERR_PATH" | sed -n 's#.*out/rawdata/\([^/]*\)/.*#\1#p')"
if [ -z "${DATASET:-}" ]; then
    echo "Could not infer dataset from stderr path" >&2
    exit 1
fi

VC_DIR="$(printf '%s\n' "$STDERR_PATH" | sed 's#/gnomad_detection/.*##')"
VCF="${VC_DIR}/${DATASET}.vcf.gz"

if [ ! -f "$VCF" ]; then
    echo "Could not find VCF: $VCF" >&2
    exit 1
fi

echo "Using VCF: $VCF" >&2
echo "Using GTF: $GTF" >&2

CSV="$OUTDIR/${DATASET}.gnomad_detection.csv"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

EXON_BED="$TMPDIR/exons.merged.bed"
PASS_VCF="$TMPDIR/${DATASET}.pass.vcf.gz"
EXONIC_VCF="$TMPDIR/${DATASET}.pass.exonic.vcf.gz"
POS_BED="$TMPDIR/${DATASET}.sites.bed"
GNOMAD_SUB="$TMPDIR/${DATASET}.gnomad.subset.vcf.gz"

awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {print $1, $4-1, $5}' "$GTF" \
    | sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge -i - > "$EXON_BED"

bcftools view -f PASS -Oz -o "$PASS_VCF" "$VCF"
tabix -f "$PASS_VCF"

bcftools view -R "$EXON_BED" -Oz -o "$EXONIC_VCF" "$PASS_VCF"
tabix -f "$EXONIC_VCF"

TOTAL=$(bcftools view -H "$EXONIC_VCF" | wc -l | awk '{print $1}')

if [ "$TOTAL" -eq 0 ]; then
    {
        echo "dataset_id,total_pass_exonic_variants,gnomad_overlap,detection_rate"
        echo "${DATASET},0,0,0"
    } > "$CSV"
    echo "Wrote: $CSV" >&2
    exit 0
fi

bcftools query -f '%CHROM\t%POS0\t%END\n' "$EXONIC_VCF" | sort -u > "$POS_BED"

bcftools view -R "$POS_BED" "$GNOMAD_URL" -Oz -o "$GNOMAD_SUB"
tabix -f "$GNOMAD_SUB"

OVERLAP=$(bcftools isec -n=2 -w1 "$EXONIC_VCF" "$GNOMAD_SUB" | bcftools view -H | wc -l | awk '{print $1}')
RATE=$(awk -v o="$OVERLAP" -v t="$TOTAL" 'BEGIN{printf "%.6f", (t>0 ? o/t : 0)}')

{
    echo "dataset_id,total_pass_exonic_variants,gnomad_overlap,detection_rate"
    echo "${DATASET},${TOTAL},${OVERLAP},${RATE}"
} > "$CSV"

echo "Wrote: $CSV" >&2
ls -lh "$OUTDIR" >&2