#!/bin/bash
set -euo pipefail

S3_URL=$1
OUTPUT_PATH=$2

mkdir -p $(dirname $OUTPUT_PATH)

aws s3 cp $S3_URL $OUTPUT_PATH
