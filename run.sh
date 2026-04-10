#!/usr/bin/env sh
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
    deep_variant|deepvariant)
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
    alignment_qc)
        exec "$DIR/alignment_qc.sh" "$@"
        ;;
    gnomad_detection)
        exec "$DIR/gnomad_detection.sh" "$@"
        ;;
    somatic_detection)
        exec "$DIR/somatic_detection.sh" "$@"
        ;;
    filter_variants)
        exec "$DIR/filter_variants.sh" "$@"
        ;;
    align)
        exec "$DIR/align.sh" "$@"
        ;;
    "")
        # collector autodetection
        for arg in "$@"; do
            case "$arg" in
                --filtered.vcf|--filtered_vcf|--filtered-vcf)
                    exec Rscript "$DIR/plot_upset.R" "$@"
                    ;;
            esac
        done

        # old-style CLI autodetection
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

        # Omnibenchmark rawdata jobs may pass no args at all
        STDERR_PATH="$(readlink -f /proc/self/fd/2 2>/dev/null || true)"
        case "$STDERR_PATH" in
            *"/out/rawdata/"*)
                exec "$DIR/download.sh"
                ;;
        esac
        ;;
    *)
        echo "ERROR: unknown task: $TASK" >&2
        exit 2
        ;;
esac

echo "ERROR: cannot determine module type from arguments" >&2
echo "ARGS: $*" >&2
exit 2