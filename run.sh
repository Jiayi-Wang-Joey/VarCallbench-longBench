#!/bin/sh
set -eu

DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

TASK=""

prev=""
for arg in "$@"; do
    if [ "$prev" = "--task" ]; then
        TASK="$arg"
        break
    fi
    prev="$arg"
done

case "$TASK" in
    clair3_rna)
        exec "$DIR/clair3_rna.sh" "$@"
        ;;
    deep_variant)
        exec "$DIR/deep_variant.sh" "$@"
        ;;
    longcallR)
        exec "$DIR/longcallR.sh" "$@"
        ;;
    longcallR_nn)
        exec "$DIR/longcallR_nn.sh" "$@"
        ;;
    plot_upset)
        exec Rscript "$DIR/plot_upset.R" "$@"
        ;;
    align)
        exec "$DIR/align.sh" "$@"
        ;;
    "")
        # backward-compatible fallback
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
        ;;
    *)
        echo "ERROR: unknown task: $TASK" >&2
        exit 2
        ;;
esac

echo "ERROR: cannot determine module type from arguments" >&2
exit 2