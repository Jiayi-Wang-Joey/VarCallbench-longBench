#!/bin/sh
set -eu

DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

for arg in "$@"; do
    case "$arg" in
        --s3_url)
            exec "$DIR/download.sh" "$@"
            ;;
        --reference_genome)
            exec "$DIR/align.sh" "$@"
            ;;
    esac
done

echo "ERROR: cannot determine module type from arguments" >&2
exit 2