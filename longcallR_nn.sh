#!/usr/bin/env bash
set -euo pipefail

DATASET=""
BAM=""
REF=""
THREADS=""
OUTDIR=""
TASK=""

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
        --task) TASK="$2"; shift 2 ;;
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
echo "TASK=${TASK}" >&2

[ -n "${DATASET}" ] || { echo "Missing --name/--dataset_id" >&2; exit 2; }
[ -n "${BAM}" ] || { echo "Missing --align.bam/--bam" >&2; exit 2; }
[ -n "${REF}" ] || { echo "Missing --reference_genome" >&2; exit 2; }
[ -n "${THREADS}" ] || { echo "Missing --threads" >&2; exit 2; }
[ -n "${OUTDIR}" ] || { echo "Missing --output_dir" >&2; exit 2; }

[ -f "${BAM}" ] || { echo "BAM not found: ${BAM}" >&2; exit 2; }
[ -f "${REF}" ] || { echo "Reference not found: ${REF}" >&2; exit 2; }

mkdir -p "${OUTDIR}"
TMPDIR="${OUTDIR}/tmp"
DATADIR="${OUTDIR}/data"
VCFDIR="${OUTDIR}/vcf"
LOGDIR="${OUTDIR}/logs"
mkdir -p "${TMPDIR}" "${DATADIR}" "${VCFDIR}" "${LOGDIR}"

dataset_lc="$(printf '%s' "${DATASET}" | tr '[:upper:]' '[:lower:]')"

if echo "${dataset_lc}" | grep -q 'pb'; then
    MODEL_CONFIG="/models/pb_masseq_config.yaml"
    MODEL_CKPT="/models/pb_masseq_model.chkpt"
    MIN_BASEQ=""
elif echo "${dataset_lc}" | grep -q 'drna'; then
    MODEL_CONFIG="/models/ont_drna_config.yaml"
    MODEL_CKPT="/models/ont_drna_model.chkpt"
    MIN_BASEQ="--min-baseq 10"
else
    MODEL_CONFIG="/models/ont_cdna_config.yaml"
    MODEL_CKPT="/models/ont_cdna_model.chkpt"
    MIN_BASEQ="--min-baseq 10"
fi

echo "MODEL_CONFIG=${MODEL_CONFIG}" >&2
echo "MODEL_CKPT=${MODEL_CKPT}" >&2
echo "MIN_BASEQ=${MIN_BASEQ}" >&2

[ -f "${MODEL_CONFIG}" ] || { echo "Missing model config: ${MODEL_CONFIG}" >&2; exit 2; }
[ -f "${MODEL_CKPT}" ] || { echo "Missing model checkpoint: ${MODEL_CKPT}" >&2; exit 2; }

FINAL_VCF_GZ="${OUTDIR}/${DATASET}.vcf.gz"

# collect contigs from BAM header, similar to per-chromosome workflow
mapfile -t CHRS < <(
    samtools idxstats "${BAM}" \
    | awk '$3 > 0 {print $1}' \
    | grep -E '^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|([1-9]|1[0-9]|2[0-2]|X|Y|MT))$'
)

[ "${#CHRS[@]}" -gt 0 ] || { echo "No contigs found in BAM header" >&2; exit 2; }

echo "Detected contigs: ${#CHRS[@]}" >&2

VCF_LIST="${TMPDIR}/vcf.list"
: > "${VCF_LIST}"

for CHR in "${CHRS[@]}"; do
    echo "=== Processing ${CHR} ===" >&2

    CHR_DATA_DIR="${DATADIR}/${CHR}"
    CHR_VCF="${VCFDIR}/${CHR}.vcf"
    CHR_LOG="${LOGDIR}/${CHR}.log"

    mkdir -p "${CHR_DATA_DIR}"

    echo "[${CHR}] longcallR_dp predict" >&2
    if [ -n "${MIN_BASEQ}" ]; then
        /usr/local/bin/longcallR_dp \
            --mode predict \
            --bam-path "${BAM}" \
            --ref-path "${REF}" \
            --threads "${THREADS}" \
            --contigs "${CHR}" \
            --output "${CHR_DATA_DIR}" \
            --min-baseq 10 \
            >> "${CHR_LOG}" 2>&1
    else
        /usr/local/bin/longcallR_dp \
            --mode predict \
            --bam-path "${BAM}" \
            --ref-path "${REF}" \
            --threads "${THREADS}" \
            --contigs "${CHR}" \
            --output "${CHR_DATA_DIR}" \
            >> "${CHR_LOG}" 2>&1
    fi

    echo "[${CHR}] longcallR_nn call" >&2
    /usr/local/bin/longcallR_nn call \
        -config "${MODEL_CONFIG}" \
        -model "${MODEL_CKPT}" \
        -data "${CHR_DATA_DIR}" \
        -ref "${REF}" \
        -output "${CHR_VCF}" \
        --no_cuda \
        -max_depth 200 \
        -batch_size 256 \
        >> "${CHR_LOG}" 2>&1

    [ -f "${CHR_VCF}" ] || { echo "Expected per-contig VCF not found: ${CHR_VCF}" >&2; exit 2; }

    printf '%s\n' "${CHR_VCF}" >> "${VCF_LIST}"
done

[ -s "${VCF_LIST}" ] || { echo "No per-contig VCFs produced" >&2; exit 2; }

echo "=== Merging VCFs ===" >&2
bcftools concat $(cat "${VCF_LIST}") 2>> "${LOGDIR}/merge.log" \
    | bcftools sort -Oz -o "${FINAL_VCF_GZ}" 2>> "${LOGDIR}/merge.log"

tabix -f -p vcf "${FINAL_VCF_GZ}" 2>> "${LOGDIR}/merge.log"

[ -f "${FINAL_VCF_GZ}" ] || { echo "Final VCF missing: ${FINAL_VCF_GZ}" >&2; exit 2; }
[ -f "${FINAL_VCF_GZ}.tbi" ] || { echo "Final VCF index missing: ${FINAL_VCF_GZ}.tbi" >&2; exit 2; }

echo "Wrote ${FINAL_VCF_GZ}" >&2