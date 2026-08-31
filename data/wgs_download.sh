#!/usr/bin/env bash
set -euo pipefail

DATASET=""
SRR=""
OUTDIR=""
THREADS="4"

echo "ARGS: $@" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --name|--dataset_id)
            DATASET="$2"; shift 2 ;;
        --srr)
            SRR="$2"; shift 2 ;;
        --output_dir)
            OUTDIR="$2"; shift 2 ;;
        --threads)
            THREADS="$2"; shift 2 ;;
        --task)
            shift 2 ;;
        *)
            echo "Unknown arg: $1" >&2
            shift ;;
    esac
done

[ -n "$DATASET" ] || { echo "Missing --name" >&2; exit 1; }
[ -n "$SRR" ]     || { echo "Missing --srr (provide verified SRR accession from NCBI)" >&2; exit 1; }
[ -n "$OUTDIR" ]  || { echo "Missing --output_dir" >&2; exit 1; }

mkdir -p "$OUTDIR"
echo "Downloading $SRR for $DATASET from ENA..." >&2

# ENA mirrors SRA runs as plain FASTQ.gz over HTTPS -- no prefetch/fasterq-dump toolchain
# needed (that path repeatedly failed late in extraction on these large WGS runs due to
# SRA's reference-compression + on-demand reference reconstruction). Path layout depends
# on accession length: 9 chars -> no extra subdir; 10/11/12 chars -> an extra subdir of
# the last 1/2/3 digits, zero-padded to 3.
PREFIX="${SRR:0:6}"
LEN=${#SRR}
case "$LEN" in
    9)  ENA_DIR="${PREFIX}/${SRR}" ;;
    10) ENA_DIR="${PREFIX}/00${SRR: -1}/${SRR}" ;;
    11) ENA_DIR="${PREFIX}/0${SRR: -2}/${SRR}" ;;
    12) ENA_DIR="${PREFIX}/${SRR: -3}/${SRR}" ;;
    *)  echo "ERROR: unsupported accession length for $SRR ($LEN chars)" >&2; exit 1 ;;
esac

BASE_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/${ENA_DIR}"
R1="${OUTDIR}/${SRR}_1.fastq.gz"
R2="${OUTDIR}/${SRR}_2.fastq.gz"

wget -c --tries=5 --timeout=60 -O "$R1" "${BASE_URL}/${SRR}_1.fastq.gz"
wget -c --tries=5 --timeout=60 -O "$R2" "${BASE_URL}/${SRR}_2.fastq.gz"

# Verify both files are actually complete, valid gzip streams before declaring success --
# a truncated/interrupted download is still a file on disk, gzip -t catches that.
gzip -t "$R1"
gzip -t "$R2"

mv "$R1" "${OUTDIR}/${DATASET}_R1.fastq.gz"
mv "$R2" "${OUTDIR}/${DATASET}_R2.fastq.gz"

echo "Done: ${OUTDIR}/${DATASET}_R1.fastq.gz" >&2
