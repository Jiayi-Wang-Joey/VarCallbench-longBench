#!/usr/bin/env bash
set -euo pipefail

DATASET=""
BAM=""
REF=""
THREADS=""
OUTDIR=""
TASK=""

echo "ARGS: $*" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --name) DATASET="$2"; shift 2 ;;
        --dataset_id) DATASET="$2"; shift 2 ;;
        --align.bam|--bam) BAM="$2"; shift 2 ;;
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

command -v samtools >/dev/null 2>&1 || { echo "samtools not found in PATH" >&2; exit 2; }
command -v bcftools >/dev/null 2>&1 || { echo "bcftools not found in PATH" >&2; exit 2; }
command -v tabix >/dev/null 2>&1 || { echo "tabix not found in PATH" >&2; exit 2; }

BAM="$(readlink -f "${BAM}")"
REF="$(readlink -f "${REF}")"
mkdir -p "${OUTDIR}"
OUTDIR="$(readlink -f "${OUTDIR}")"

echo "ABS_BAM=${BAM}" >&2
echo "ABS_REF=${REF}" >&2
echo "ABS_OUTDIR=${OUTDIR}" >&2

TMPDIR="${OUTDIR}/tmp"
DATADIR="${OUTDIR}/data"
VCFDIR="${OUTDIR}/vcf"
LOGDIR="${OUTDIR}/logs"
mkdir -p "${TMPDIR}" "${DATADIR}" "${VCFDIR}" "${LOGDIR}"

dataset_lc="$(printf '%s' "${DATASET}" | tr '[:upper:]' '[:lower:]')"

if echo "${dataset_lc}" | grep -q 'pb'; then
    MODEL_CONFIG="/models/pb_masseq_config.yaml"
    MODEL_CKPT="/models/pb_masseq_model.chkpt"
    MIN_BASEQ_FLAG=0
elif echo "${dataset_lc}" | grep -q 'drna'; then
    MODEL_CONFIG="/models/ont_drna_config.yaml"
    MODEL_CKPT="/models/ont_drna_model.chkpt"
    MIN_BASEQ_FLAG=1
else
    MODEL_CONFIG="/models/ont_cdna_config.yaml"
    MODEL_CKPT="/models/ont_cdna_model.chkpt"
    MIN_BASEQ_FLAG=1
fi

echo "MODEL_CONFIG=${MODEL_CONFIG}" >&2
echo "MODEL_CKPT=${MODEL_CKPT}" >&2
echo "MIN_BASEQ_FLAG=${MIN_BASEQ_FLAG}" >&2

[ -f "${MODEL_CONFIG}" ] || { echo "Missing model config: ${MODEL_CONFIG}" >&2; exit 2; }
[ -f "${MODEL_CKPT}" ] || { echo "Missing model checkpoint: ${MODEL_CKPT}" >&2; exit 2; }

[ -x /usr/local/bin/longcallR_dp ] || { echo "Missing executable: /usr/local/bin/longcallR_dp" >&2; exit 2; }
[ -x /usr/local/bin/longcallR_nn ] || { echo "Missing executable: /usr/local/bin/longcallR_nn" >&2; exit 2; }

FINAL_VCF_GZ="${OUTDIR}/${DATASET}.vcf.gz"

if ! samtools idxstats "${BAM}" > "${LOGDIR}/samtools_idxstats.initial.txt" 2> "${LOGDIR}/samtools_idxstats.initial.log"; then
    echo "BAM index missing or unusable for: ${BAM}" >&2
    echo "Attempting to create BAM index with samtools index..." >&2

    samtools index -@ "${THREADS}" "${BAM}" > "${LOGDIR}/samtools_index.log" 2>&1 || {
        rc=$?
        echo "samtools index failed with exit code ${rc}" >&2
        tail -n 50 "${LOGDIR}/samtools_index.log" >&2 || true
        exit "${rc}"
    }

    samtools idxstats "${BAM}" > "${LOGDIR}/samtools_idxstats.after_index.txt" 2> "${LOGDIR}/samtools_idxstats.after_index.log" || {
        rc=$?
        echo "samtools idxstats still failed after indexing for BAM: ${BAM}" >&2
        tail -n 50 "${LOGDIR}/samtools_idxstats.after_index.log" >&2 || true
        exit "${rc}"
    }
fi

mapfile -t CHRS < <(
    samtools idxstats "${BAM}" 2> "${LOGDIR}/samtools_idxstats.contigs.log" \
    | awk '$3 > 0 {print $1}' \
    | grep -E '^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|([1-9]|1[0-9]|2[0-2]|X|Y|MT))$'
)

[ "${#CHRS[@]}" -gt 0 ] || {
    echo "No primary contigs with mapped reads found in BAM: ${BAM}" >&2
    exit 2
}

if [ -n "${DEBUG_CHR:-}" ]; then
    echo "DEBUG_CHR set: restricting to ${DEBUG_CHR}" >&2
    found=0
    for c in "${CHRS[@]}"; do
        if [ "${c}" = "${DEBUG_CHR}" ]; then
            found=1
            CHRS=("${DEBUG_CHR}")
            break
        fi
    done
    [ "${found}" -eq 1 ] || {
        echo "DEBUG_CHR ${DEBUG_CHR} not found among selected contigs" >&2
        printf 'Selected contigs: %s\n' "${CHRS[*]}" >&2
        exit 2
    }
fi

echo "Detected primary contigs with mapped reads: ${#CHRS[@]}" >&2
printf '%s\n' "${CHRS[@]}" >&2

VCF_LIST="${TMPDIR}/vcf.list"
: > "${VCF_LIST}"

process_chr() {
    local CHR="$1"
    local CHR_DATA_DIR="${DATADIR}/${CHR}"
    local CHR_VCF="${VCFDIR}/${CHR}.vcf"
    local DP_LOG="${LOGDIR}/${CHR}.longcallR_dp.log"
    local NN_LOG="${LOGDIR}/${CHR}.longcallR_nn.log"

    echo "=== Processing ${CHR} ===" >&2

    mkdir -p "${CHR_DATA_DIR}"
    rm -f "${CHR_VCF}"

    if [ -f "${BAM}.bai" ]; then
        :
    elif [ -f "${BAM%.bam}.bai" ]; then
        :
    elif [ -f "${BAM}.csi" ]; then
        :
    else
        echo "[${CHR}] No BAM index found for ${BAM}" >&2
        exit 2
    fi

    echo "[${CHR}] longcallR_dp predict" >&2
    if [ "${MIN_BASEQ_FLAG}" -eq 1 ]; then
        /usr/local/bin/longcallR_dp \
            --mode predict \
            --bam-path "${BAM}" \
            --ref-path "${REF}" \
            --threads 1 \
            --contigs "${CHR}" \
            --output "${CHR_DATA_DIR}" \
            --min-baseq 10 \
            > "${DP_LOG}" 2>&1
    else
        /usr/local/bin/longcallR_dp \
            --mode predict \
            --bam-path "${BAM}" \
            --ref-path "${REF}" \
            --threads 1 \
            --contigs "${CHR}" \
            --output "${CHR_DATA_DIR}" \
            > "${DP_LOG}" 2>&1
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
        > "${NN_LOG}" 2>&1

    [ -f "${CHR_VCF}" ] || {
        echo "[${CHR}] Expected per-contig VCF not found: ${CHR_VCF}" >&2
        exit 2
    }

    [ -s "${CHR_VCF}" ] || {
        echo "[${CHR}] Per-contig VCF is empty: ${CHR_VCF}" >&2
        exit 2
    }

    echo "${CHR_VCF}"
}

export BAM REF DATADIR VCFDIR LOGDIR MODEL_CONFIG MODEL_CKPT MIN_BASEQ_FLAG
export -f process_chr

JOBLOG="${LOGDIR}/parallel_jobs.log"
: > "${JOBLOG}"

pids=()
tmp_vcf_lists=()

for CHR in "${CHRS[@]}"; do
    tmp_vcf_list_chr="${TMPDIR}/${CHR}.vcf.path.txt"
    tmp_vcf_lists+=("${tmp_vcf_list_chr}")

    (
        process_chr "${CHR}" > "${tmp_vcf_list_chr}"
    ) >> "${JOBLOG}" 2>&1 &

    pids+=("$!")

    while [ "$(jobs -rp | wc -l)" -ge "${THREADS}" ]; do
        sleep 1
    done
done

fail=0
for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        fail=1
    fi
done

[ "${fail}" -eq 0 ] || {
    echo "At least one chromosome job failed. Check logs in ${LOGDIR}" >&2
    exit 1
}

: > "${VCF_LIST}"
for f in "${tmp_vcf_lists[@]}"; do
    [ -s "${f}" ] || {
        echo "Missing chromosome VCF path file: ${f}" >&2
        exit 2
    }
    cat "${f}" >> "${VCF_LIST}"
done

[ -s "${VCF_LIST}" ] || {
    echo "No per-contig VCFs produced" >&2
    exit 2
}

echo "=== Merging VCFs ===" >&2
bcftools concat -f "${VCF_LIST}" 2>> "${LOGDIR}/merge.log" \
    | bcftools sort -Oz -o "${FINAL_VCF_GZ}" 2>> "${LOGDIR}/merge.log"

tabix -f -p vcf "${FINAL_VCF_GZ}" 2>> "${LOGDIR}/merge.log"

[ -f "${FINAL_VCF_GZ}" ] || { echo "Final VCF missing: ${FINAL_VCF_GZ}" >&2; exit 2; }
[ -f "${FINAL_VCF_GZ}.tbi" ] || { echo "Final VCF index missing: ${FINAL_VCF_GZ}.tbi" >&2; exit 2; }

echo "Wrote ${FINAL_VCF_GZ}" >&2