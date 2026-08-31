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
        exec "$DIR/variant_call/clair3_rna.sh" "$@"
        ;;
    deep_variant|deepvariant)
        exec "$DIR/variant_call/deep_variant.sh" "$@"
        ;;
    longcallR)
        exec "$DIR/variant_call/longcallR.sh" "$@"
        ;;
    longcallR_nn)
        exec "$DIR/variant_call/longcallR_nn.sh" "$@"
        ;;
    align)
        exec "$DIR/align/align.sh" "$@"
        ;;
    wgs_align)
        exec "$DIR/align/wgs_align.sh" "$@"
        ;;
    alignment_qc)
        exec "$DIR/align/alignment_qc.sh" "$@"
        ;;
    isolaser_annotate)
        exec "$DIR/variant_call/isolaser_annotate.sh" "$@"
        ;;
    isolaser_run)
        exec "$DIR/variant_call/isolaser_run.sh" "$@"
        ;;
    isolaser_filter)
        exec "$DIR/variant_call/isolaser_filter.sh" "$@"
        ;;
    wgs_markdup)
        exec "$DIR/gatk/markdup.sh" "$@"
        ;;
    wgs_bqsr)
        exec "$DIR/gatk/bqsr.sh" "$@"
        ;;
    wgs_haplotypecaller)
        exec "$DIR/gatk/haplotypecaller.sh" "$@"
        ;;
    wgs_genotypegvcfs)
        exec "$DIR/gatk/genotypegvcfs.sh" "$@"
        ;;
    wgs_hard_filter)
        exec "$DIR/gatk/hard_filter.sh" "$@"
        ;;
    wgs_truth_promote)
        exec "$DIR/gatk/wgs_truth_promote.sh" "$@"
        ;;
    bam_index)
        exec "$DIR/eval/bam_index.sh" "$@"
        ;;
    mosdepth_coverage)
        exec "$DIR/eval/mosdepth_coverage.sh" "$@"
        ;;
    generate_bed)
        exec "$DIR/eval/generate_bed.sh" "$@"
        ;;
    happy_benchmark)
        exec "$DIR/eval/happy_benchmark.sh" "$@"
        ;;
    happy_summary_collector)
        exec "$DIR/eval/happy_summary_collector.R" "$@"
        ;;
    "")
        # old-style CLI autodetection
        for arg in "$@"; do
            case "$arg" in
                --s3_url)
                    exec "$DIR/data/download.sh" "$@"
                    ;;
                --reference_genome)
                    exec "$DIR/align/align.sh" "$@"
                    ;;
                --srr)
                    exec "$DIR/data/wgs_download.sh" "$@"
                    ;;
            esac
        done

        # Omnibenchmark rawdata jobs may pass no args at all
        STDERR_PATH="$(readlink -f /proc/self/fd/2 2>/dev/null || true)"
        case "$STDERR_PATH" in
            *"/out/rawdata/"*)
                exec "$DIR/data/download.sh"
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