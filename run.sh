#!/bin/sh
set -eu

TASK=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --task)
            TASK="$2"
            shift 2
            ;;
        *)
            break
            ;;
    esac
done

if [ -z "$TASK" ]; then
    echo "ERROR: --task is required" >&2
    exit 2
fi

case "$TASK" in
    download)
        exec ./download.sh "$@"
        ;;
    align)
        exec ./align.sh "$@"
        ;;
    *)
        echo "ERROR: unknown task: $TASK" >&2
        exit 2
        ;;
esac