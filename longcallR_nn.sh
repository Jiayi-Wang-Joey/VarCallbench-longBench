#!/usr/bin/env bash
set -euo pipefail

DATASET=""
BAM=""
REF=""
THREADS=""
OUTDIR=""

echo "ARGS: $@" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --name) DATASET="$2"; shift 2 ;;
        --dataset_id) DATASET="$2"; shift 2 ;;
        --align.bam) BAM="$2"; shift 2 ;;
        --bam) BAM="$2"; shift 2 ;;
        --reference_genome) REF="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --output_dir) OUTDIR="$2"; shift 2 ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

echo "DATASET=${DATASET}" >&2
echo "BAM=${BAM}" >&2
echo "REF=${REF}" >&2
echo "THREADS=${THREADS}" >&2
echo "OUTDIR=${OUTDIR}" >&2

[ -n "${DATASET}" ] || { echo "Missing --name/--dataset_id" >&2; exit 2; }
[ -n "${BAM}" ] || { echo "Missing --align.bam/--bam" >&2; exit 2; }
[ -n "${REF}" ] || { echo "Missing --reference_genome" >&2; exit 2; }
[ -n "${THREADS}" ] || { echo "Missing --threads" >&2; exit 2; }
[ -n "${OUTDIR}" ] || { echo "Missing --output_dir" >&2; exit 2; }

mkdir -p "${OUTDIR}"
TMPDIR="${OUTDIR}/tmp"
mkdir -p "${TMPDIR}"

# choose model by dataset name
dataset_lc="$(printf '%s' "${DATASET}" | tr '[:upper:]' '[:lower:]')"

if echo "${dataset_lc}" | grep -q 'pb'; then
    MODEL_CONFIG="/models/pb_masseq_config.yaml"
    MODEL_CKPT="/models/pb_masseq_model.chkpt"
elif echo "${dataset_lc}" | grep -q 'drna'; then
    MODEL_CONFIG="/models/ont_drna_config.yaml"
    MODEL_CKPT="/models/ont_drna_model.chkpt"
else
    MODEL_CONFIG="/models/ont_cdna_config.yaml"
    MODEL_CKPT="/models/ont_cdna_model.chkpt"
fi

echo "MODEL_CONFIG=${MODEL_CONFIG}" >&2
echo "MODEL_CKPT=${MODEL_CKPT}" >&2

[ -f "${MODEL_CONFIG}" ] || { echo "Missing model config: ${MODEL_CONFIG}" >&2; exit 2; }
[ -f "${MODEL_CKPT}" ] || { echo "Missing model checkpoint: ${MODEL_CKPT}" >&2; exit 2; }

RAW_VCF="${TMPDIR}/${DATASET}.raw.vcf"
FINAL_VCF_GZ="${OUTDIR}/${DATASET}.vcf.gz"

# example command - adapt to actual CLI of longcallR_nn
longcallR_nn \
    --bam "${BAM}" \
    --ref "${REF}" \
    --thread "${THREADS}" \
    --model_config "${MODEL_CONFIG}" \
    --model "${MODEL_CKPT}" \
    --output "${RAW_VCF}"

[ -f "${RAW_VCF}" ] || { echo "Expected output not found: ${RAW_VCF}" >&2; exit 2; }

bgzip -c "${RAW_VCF}" > "${FINAL_VCF_GZ}"
tabix -f -p vcf "${FINAL_VCF_GZ}"

echo "Wrote ${FINAL_VCF_GZ}" >&2