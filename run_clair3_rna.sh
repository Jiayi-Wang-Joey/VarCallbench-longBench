#!/usr/bin/env bash
set -euo pipefail

DATASET=""
BAM=""
REF=""
THREADS=""
OUTVCF=""

while [ $# -gt 0 ]; do
    case "$1" in
        --dataset_id) DATASET="$2"; shift 2 ;;
        --bam) BAM="$2"; shift 2 ;;
        --reference_genome) REF="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --output_vcf) OUTVCF="$2"; shift 2 ;;
        *) shift ;;
    esac
done

dataset_lc=$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')

if [[ "$dataset_lc" == *pb* ]]; then
    PLATFORM="hifi_mas_minimap2"
elif [[ "$dataset_lc" == *ont* && "$dataset_lc" == *drna* ]]; then
    PLATFORM="ont_dorado_drna004"
elif [[ "$dataset_lc" == *ont* ]]; then
    PLATFORM="ont_r10_dorado_cdna"
else
    echo "ERROR: cannot infer Clair3-RNA platform from dataset: $DATASET" >&2
    exit 1
fi

echo "DATASET=$DATASET"
echo "PLATFORM=$PLATFORM"

TOOL_IMAGE="/home/jiayiwang/tools/clair3-rna-latest.simg"
OUTDIR=$(dirname "$OUTVCF")/clair3_tmp_${DATASET}
mkdir -p "$OUTDIR"

if [ ! -f "${BAM}.bai" ]; then
    samtools index "$BAM"
fi

singularity exec \
    -B "$(dirname "$BAM")":"$(dirname "$BAM")" \
    -B "$(dirname "$REF")":"$(dirname "$REF")" \
    -B "$(dirname "$OUTDIR")":"$(dirname "$OUTDIR")" \
    "$TOOL_IMAGE" \
    /bin/bash -c "
        source /opt/conda/bin/activate /opt/conda/envs/clair3_rna && \
        /opt/bin/run_clair3_rna \
            --bam_fn '$BAM' \
            --ref_fn '$REF' \
            --threads '$THREADS' \
            --platform '$PLATFORM' \
            --tag_variant_using_readiportal \
            --remove_intermediate_dir \
            --output_dir '$OUTDIR' \
            --conda_prefix /opt/conda/envs/clair3_rna
    "

mv "$OUTDIR/output.vcf.gz" "$OUTVCF"