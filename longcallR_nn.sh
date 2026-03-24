#!/usr/bin/env bash
set -euo pipefail

DATASET=""
BAM=""
REF=""
THREADS=""
OUTDIR=""
MODEL_DIR="/models"
NO_CUDA="true"

while [ $# -gt 0 ]; do
    case "$1" in
        --task)
            shift 2
            ;;
        --name|--dataset_id)
            DATASET="$2"
            shift 2
            ;;
        --align.bam|--bam)
            BAM="$2"
            shift 2
            ;;
        --reference_genome)
            REF="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --output_dir)
            OUTDIR="$2"
            shift 2
            ;;
        --model_dir)
            MODEL_DIR="$2"
            shift 2
            ;;
        --no_cuda)
            NO_CUDA="$2"
            shift 2
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [ -z "$DATASET" ] || [ -z "$BAM" ] || [ -z "$REF" ] || [ -z "$THREADS" ] || [ -z "$OUTDIR" ]; then
    echo "ERROR: missing required arguments" >&2
    echo "DATASET=$DATASET" >&2
    echo "BAM=$BAM" >&2
    echo "REF=$REF" >&2
    echo "THREADS=$THREADS" >&2
    echo "OUTDIR=$OUTDIR" >&2
    exit 2
fi

pick_longcallr_nn_model() {
    local sample="$1"
    local model_dir="$2"
    local s
    s=$(printf '%s' "$sample" | tr '[:upper:]' '[:lower:]')

    if [[ "$s" == *"drna"* ]]; then
        echo "${model_dir}/ont_drna_config.yaml ${model_dir}/ont_drna_model.chkpt"
    elif [[ "$s" == *"pb"* ]]; then
        echo "${model_dir}/pb_masseq_config.yaml ${model_dir}/pb_masseq_model.chkpt"
    elif [[ "$s" == *"ont"* ]]; then
        echo "${model_dir}/ont_cdna_config.yaml ${model_dir}/ont_cdna_model.chkpt"
    else
        echo "ERROR: sample '$sample' does not contain pb / ont / drna" >&2
        exit 1
    fi
}

read -r CFG CKPT < <(pick_longcallr_nn_model "$DATASET" "$MODEL_DIR")

MIN_BASEQ_ARGS=()
DATASET_LC=$(printf '%s' "$DATASET" | tr '[:upper:]' '[:lower:]')
if [[ "$DATASET_LC" != *"pb"* ]]; then
    MIN_BASEQ_ARGS=(--min-baseq 10)
fi

mkdir -p "$OUTDIR"
WORKDIR="${OUTDIR}/work"
DATADIR="${WORKDIR}/data"
VCFDIR="${WORKDIR}/vcf"
LOGDIR="${WORKDIR}/logs"
mkdir -p "$DATADIR" "$VCFDIR" "$LOGDIR"

if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
    samtools index "$BAM"
fi

CONTIGS=$(samtools idxstats "$BAM" | awk '$1 != "*" && $3 > 0 {print $1}')

if [ -z "$CONTIGS" ]; then
    echo "ERROR: no contigs with reads found in $BAM" >&2
    exit 1
fi

for CHR in $CONTIGS; do
    mkdir -p "${DATADIR}/${CHR}"

    /opt/longcallR-nn/longcallR_dp/target/release/longcallR_dp --mode predict \
        --bam-path "$BAM" \
        --ref-path "$REF" \
        --threads "$THREADS" \
        --contigs "$CHR" \
        --output "${DATADIR}/${CHR}" \
        "${MIN_BASEQ_ARGS[@]}" \
        > "${LOGDIR}/${CHR}.predict.log" 2>&1

    CALL_ARGS=()
    if [ "$NO_CUDA" = "true" ]; then
        CALL_ARGS+=(--no_cuda)
    fi

    /opt/conda/bin/longcallR_nn call \
        -config "$CFG" \
        -model "$CKPT" \
        -data "${DATADIR}/${CHR}" \
        -ref "$REF" \
        -output "${VCFDIR}/${CHR}.vcf" \
        "${CALL_ARGS[@]}" \
        -max_depth 200 \
        -batch_size 256 \
        > "${LOGDIR}/${CHR}.call.log" 2>&1
done

find "$VCFDIR" -name "*.vcf" | sort > "${WORKDIR}/vcf.list"

if [ ! -s "${WORKDIR}/vcf.list" ]; then
    echo "ERROR: no per-contig VCF files were generated" >&2
    exit 1
fi

bcftools concat -f "${WORKDIR}/vcf.list" \
    | bcftools sort -Oz -o "${OUTDIR}/${DATASET}.vcf.gz"

tabix -p vcf "${OUTDIR}/${DATASET}.vcf.gz"