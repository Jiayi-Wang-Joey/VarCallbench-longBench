#!/bin/sh
set -eu

S3_URL=""
OUTPUT_DIR=""
NAME=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --s3_url)
            S3_URL="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --name)
            NAME="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

if [ -z "$S3_URL" ]; then
    echo "ERROR: --s3_url is required" >&2
    exit 2
fi

if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: --output_dir is required" >&2
    exit 2
fi

if [ -z "$NAME" ]; then
    echo "ERROR: --name is required" >&2
    exit 2
fi

mkdir -p "$OUTPUT_DIR"
aws s3 cp "$S3_URL" "$OUTPUT_DIR/$NAME.fastq.gz" --no-sign-request