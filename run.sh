#!/bin/sh
set -eu

DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

TASK=""

while [ $# -gt 0 ]; do
    case "$1" in
        --task)
            if [ $# -lt 2 ]; then
                echo "ERROR: --task requires a value" >&2
                exit 2
            fi
            TASK="$2"
            break
            ;;
    esac
    shift
done

dispatch_collector() {
    for arg in "$@"; do
        case "$arg" in
            --variant.vcf|--variant_vcf|--variant-vcf)
                exec Rscript "$DIR/plot_upset.R" "$@"
                ;;
        esac
    done
}

dispatch_legacy() {
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
}

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
    alignment_qc)
        exec "$DIR/alignment_qc.sh" "$@"
        ;;
    gnomad_detection)
        exec "$DIR/gnomad_detection.sh" "$@"
        ;;
    somatic_detection)
        exec "$DIR/somatic_detection.sh" "$@"
        ;;
    align)
        exec "$DIR/align.sh" "$@"
        ;;
    "")
        dispatch_collector "$@"
        dispatch_legacy "$@"
        ;;
    *)
        echo "ERROR: unknown task: $TASK" >&2
        exit 2
        ;;
esac

echo "ERROR: cannot determine module type from arguments" >&2
echo "ARGS: $*" >&2
exit 2