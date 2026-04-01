#!/usr/bin/env bash
set -euo pipefail

SOMATIC_VCF=""
THREADS="1"

echo "ARGS: $*" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --task)
            shift 2
            ;;
        --somatic_vcf|--somatic-vcf)
            SOMATIC_VCF="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
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

if [ -z "$SOMATIC_VCF" ] && [ -f "$PARAM_JSON" ]; then
    SOMATIC_VCF="$(python - <<'PY' "$PARAM_JSON"
import json, sys

def find_key(obj, key):
    if isinstance(obj, dict):
        if key in obj and isinstance(obj[key], str):
            return obj[key]
        for v in obj.values():
            out = find_key(v, key)
            if out:
                return out
    elif isinstance(obj, list):
        for v in obj:
            out = find_key(v, key)
            if out:
                return out
    return ""

with open(sys.argv[1]) as fh:
    data = json.load(fh)

print(find_key(data, "somatic_vcf"))
PY
)"
fi

echo "SOMATIC_VCF: $SOMATIC_VCF" >&2

if [ -z "$SOMATIC_VCF" ]; then
    echo "Missing somatic_vcf parameter" >&2
    exit 1
fi

if [ ! -f "$SOMATIC_VCF" ]; then
    echo "Could not find somatic VCF: $SOMATIC_VCF" >&2
    exit 1
fi

if [ ! -f "${SOMATIC_VCF}.tbi" ] && [ ! -f "${SOMATIC_VCF}.csi" ]; then
    echo "Could not find index for somatic VCF: ${SOMATIC_VCF}.tbi or .csi" >&2
    exit 1
fi

DATASET="$(printf '%s\n' "$STDERR_PATH" | sed -n 's#.*out/rawdata/\([^/]*\)/.*#\1#p')"
if [ -z "${DATASET:-}" ]; then
    echo "Could not infer dataset from stderr path" >&2
    exit 1
fi

VC_DIR="$(printf '%s\n' "$STDERR_PATH" | sed 's#/somatic_detection/.*##')"
VCF="${VC_DIR}/${DATASET}.vcf.gz"

if [ ! -f "$VCF" ]; then
    echo "Could not find VCF: $VCF" >&2
    exit 1
fi

echo "Using VCF: $VCF" >&2
echo "Using somatic truth VCF: $SOMATIC_VCF" >&2

CSV_OUT="$OUTDIR/${DATASET}.somatic_detection.csv"

WORKDIR="$OUTDIR/tmp"
mkdir -p "$WORKDIR"
export TMPDIR="$WORKDIR"

echo "WORKDIR: $WORKDIR" >&2
echo "TMPDIR: $TMPDIR" >&2

TRUTH_TSV="$WORKDIR/${DATASET}.somatic_truth.tsv"
CALLS_TSV="$WORKDIR/${DATASET}.calls.tsv"
TRUTH_KEY="$WORKDIR/${DATASET}.somatic_truth.key"
CALLS_KEY="$WORKDIR/${DATASET}.calls.key"
MATCH_KEY="$WORKDIR/${DATASET}.matches.key"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$SOMATIC_VCF" \
    | sort -T "$WORKDIR" -u > "$TRUTH_TSV"

bcftools view -f PASS "$VCF" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    | sort -T "$WORKDIR" -u > "$CALLS_TSV"

awk '{print $1":"$2":"$3":"$4}' "$TRUTH_TSV" \
    | sort -T "$WORKDIR" -u > "$TRUTH_KEY"

awk '{print $1":"$2":"$3":"$4}' "$CALLS_TSV" \
    | sort -T "$WORKDIR" -u > "$CALLS_KEY"

comm -12 "$TRUTH_KEY" "$CALLS_KEY" > "$MATCH_KEY"

TOTAL=$(wc -l < "$TRUTH_KEY" | awk '{print $1}')
DETECTED=$(wc -l < "$MATCH_KEY" | awk '{print $1}')
RATE=$(awk -v d="$DETECTED" -v t="$TOTAL" 'BEGIN{printf "%.6f", (t>0 ? d/t : 0)}')

{
    echo "dataset_id,total_somatic_variants,detected_somatic_variants,detection_rate"
    echo "${DATASET},${TOTAL},${DETECTED},${RATE}"
} > "$CSV_OUT"

echo "Wrote: $CSV_OUT" >&2
ls -lh "$OUTDIR" >&2