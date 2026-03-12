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
*) shift ;;
esac
done
if [ -z "$BAM" ]; then echo "ERROR: bam was not provided" >&2; exit 1; fi
if [ -z "$REF" ]; then echo "ERROR: reference_genome was not provided" >&2; exit 1; fi
if [ -z "$THREADS" ]; then echo "ERROR: threads was not provided" >&2; exit 1; fi
if [ -z "$OUTDIR" ]; then echo "ERROR: output_dir was not provided" >&2; exit 1; fi
if [ -z "$DATASET" ] || [ "$DATASET" = "{dataset}" ]; then
    DATASET=$(basename "$BAM"); DATASET=${DATASET%.aligned.bam}; DATASET=${DATASET%.bam}
fi
echo "DATASET=$DATASET" >&2
echo "BAM=$BAM" >&2
echo "REF=$REF" >&2
echo "THREADS=$THREADS" >&2
echo "OUTDIR=$OUTDIR" >&2
dataset_lc=$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')
if [[ "$dataset_lc" == *pb* ]]; then
    PLATFORM="hifi_mas_minimap2"
elif [[ "$dataset_lc" == *ont* && "$dataset_lc" == *drna* ]]; then
    PLATFORM="ont_dorado_drna004"
elif [[ "$dataset_lc" == *ont* ]]; then
    PLATFORM="ont_r10_dorado_cdna"
else
    echo "ERROR: cannot infer Clair3-RNA platform from dataset: $DATASET" >&2; exit 1
fi
echo "PLATFORM=$PLATFORM" >&2

TMPDIR="${OUTDIR}/clair3_tmp_${DATASET}"
FINALVCF="${OUTDIR}/${DATASET}.vcf.gz"
mkdir -p "$TMPDIR" "$OUTDIR"

if [ ! -f "${BAM}.bai" ]; then
    samtools index "$BAM"
fi

TOOL_IMAGE="/home/jiayiwang/VarCallbench/envs/clair3-rna-latest.simg"

set -x
singularity exec \
    -B "$(dirname "$BAM")":"$(dirname "$BAM")" \
    -B "$(dirname "$REF")":"$(dirname "$REF")" \
    -B "$TMPDIR":"$TMPDIR" \
    -B "$OUTDIR":"$OUTDIR" \
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
            --output_dir '$TMPDIR' \
            --conda_prefix /opt/conda/envs/clair3_rna
    "
set +x

if [ ! -f "${TMPDIR}/output.vcf.gz" ]; then
    echo "ERROR: expected output not found: ${TMPDIR}/output.vcf.gz" >&2; exit 1
fi
mv "${TMPDIR}/output.vcf.gz" "$FINALVCF"