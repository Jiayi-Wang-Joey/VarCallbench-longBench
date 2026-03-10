#!/bin/bash
set -euo pipefail

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --s3_url) S3_URL="$2"; shift ;;
    esac
    shift
done

aws s3 cp "$S3_URL" "$OMNI_OUTPUT_rawdata_fastq" --no-sign-request
