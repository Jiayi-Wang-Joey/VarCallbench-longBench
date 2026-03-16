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
    DATASET=$(basename "$BAM")
    DATASET=${DATASET%.aligned.bam}
    DATASET=${DATASET%.bam}
fi

echo "DATASET=$DATASET" >&2
echo "BAM=$BAM" >&2
echo "REF=$REF" >&2
echo "THREADS=$THREADS" >&2
echo "OUTDIR=$OUTDIR" >&2

dataset_lc=$(echo "$DATASET" | tr '[:upper:]' '[:lower:]')

if [[ "$dataset_lc" == *pb*]]; then
    MODEL="MASSEQ"
    EXTRA_ARGS=()
elif [[ "$dataset_lc" == *ont* ]]; then
    MODEL="ONT_R104"
    EXTRA_ARGS=(--make_examples_extra_args "split_skip_reads=true")
else
    echo "ERROR: cannot infer DeepVariant model from dataset: $DATASET" >&2
    exit 1
fi

echo "MODEL=$MODEL" >&2

TMPDIR="${OUTDIR}/deepvariant_tmp_${DATASET}"
FINALVCF="${OUTDIR}/${DATASET}.vcf.gz"
mkdir -p "$TMPDIR" "$OUTDIR"

if [ ! -f "${BAM}.bai" ]; then
    samtools index "$BAM"
fi

set -x

/opt/deepvariant/bin/run_deepvariant \
    --model_type "$MODEL" \
    --ref "$REF" \
    --reads "$BAM" \
    --output_vcf "$FINALVCF" \
    --num_shards "$THREADS" \
    --intermediate_results_dir "$TMPDIR" \
    "${EXTRA_ARGS[@]}"

set +x

if [ ! -f "$FINALVCF" ]; then
    echo "ERROR: expected output not found: $FINALVCF" >&2
    exit 1
fi