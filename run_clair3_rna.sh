#!/usr/bin/env bash
set -euo pipefail

DATASET="$1"
BAM="$2"
REF="$3"
THREADS="$4"
OUTVCF="$5"

dataset_lc=$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')

if [[ "$dataset_lc" == *pb* ]]; then
    PLATFORM="hifi_mas_minimap2"

elif [[ "$dataset_lc" == *ont* && "$dataset_lc" == *drna* ]]; then
    PLATFORM="ont_dorado_drna004"

elif [[ "$dataset_lc" == *ont* ]]; then
    PLATFORM="ont_r10_dorado_cdna"

else
    echo "ERROR: cannot infer Clair3-RNA platform from dataset: $DATASET"
    exit 1
fi

echo "Dataset: $DATASET"
echo "Detected platform: $PLATFORM"

OUTDIR=$(dirname "$OUTVCF")/clair3_tmp_${DATASET}
mkdir -p "$OUTDIR"

if [ ! -f "${BAM}.bai" ]; then
    samtools index "$BAM"
fi

source /opt/conda/bin/activate /opt/conda/envs/clair3_rna

/opt/bin/run_clair3_rna \
    --bam_fn "$BAM" \
    --ref_fn "$REF" \
    --threads "$THREADS" \
    --platform "$PLATFORM" \
    --tag_variant_using_readiportal \
    --remove_intermediate_dir \
    --output_dir "$OUTDIR" \
    --conda_prefix /opt/conda/envs/clair3_rna

mv "$OUTDIR/output.vcf.gz" "$OUTVCF"