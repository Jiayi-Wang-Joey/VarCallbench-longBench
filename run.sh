#!/bin/sh
set -eu

# detect stage from arguments
for arg in "$@"; do
    case "$arg" in
        --s3_url)
            exec ./download.sh "$@"
            ;;
        --reference_genome)
            exec ./align.sh "$@"
            ;;
    esac
done

echo "ERROR: cannot determine module type from arguments" >&2
exit 2