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

dispatch_from_parameters_json() {
    STDERR_PATH="$(readlink -f /proc/self/fd/2 2>/dev/null || true)"
    [ -n "$STDERR_PATH" ] || return 1

    OUTDIR="$(dirname "$STDERR_PATH")"
    PARAM_JSON="$OUTDIR/parameters.json"
    [ -f "$PARAM_JSON" ] || return 1

    # extract task if present
    JSON_TASK="$(python - <<'PY' "$PARAM_JSON"
import json, sys

def find_key(obj, key):
    if isinstance(obj, dict):
        if key in obj and isinstance(obj[key], str):
            return obj[key]
        for v in obj.values():
            out = find_key(v, key)
            if out:
                return out
    elif isinstance(obj, list):
        for v in obj:
            out = find_key(v, key)
            if out:
                return out
    return ""

with open(sys.argv[1]) as fh:
    data = json.load(fh)

print(find_key(data, "task"))
PY
)"

    if [ -n "$JSON_TASK" ]; then
        case "$JSON_TASK" in
            clair3_rna)
                exec "$DIR/clair3_rna.sh" --task "$JSON_TASK"
                ;;
            deep_variant)
                exec "$DIR/deep_variant.sh" --task "$JSON_TASK"
                ;;
            longcallR)
                exec "$DIR/longcallR.sh" --task "$JSON_TASK"
                ;;
            longcallR_nn)
                exec "$DIR/longcallR_nn.sh" --task "$JSON_TASK"
                ;;
            alignment_qc)
                exec "$DIR/alignment_qc.sh" --task "$JSON_TASK"
                ;;
            gnomad_detection)
                exec "$DIR/gnomad_detection.sh" --task "$JSON_TASK"
                ;;
            somatic_detection)
                exec "$DIR/somatic_detection.sh" --task "$JSON_TASK"
                ;;
            align)
                exec "$DIR/align.sh" --task "$JSON_TASK"
                ;;
        esac
    fi

    # rawdata fallback: detect s3_url in parameters.json
    JSON_S3_URL="$(python - <<'PY' "$PARAM_JSON"
import json, sys

def find_key(obj, key):
    if isinstance(obj, dict):
        if key in obj and isinstance(obj[key], str):
            return obj[key]
        for v in obj.values():
            out = find_key(v, key)
            if out:
                return out
    elif isinstance(obj, list):
        for v in obj:
            out = find_key(v, key)
            if out:
                return out
    return ""

with open(sys.argv[1]) as fh:
    data = json.load(fh)

print(find_key(data, "s3_url"))
PY
)"

    if [ -n "$JSON_S3_URL" ]; then
        exec "$DIR/download.sh" --s3_url "$JSON_S3_URL"
    fi

    # align fallback if reference genome present
    JSON_REF="$(python - <<'PY' "$PARAM_JSON"
import json, sys

def find_key(obj, key):
    if isinstance(obj, dict):
        if key in obj and isinstance(obj[key], str):
            return obj[key]
        for v in obj.values():
            out = find_key(v, key)
            if out:
                return out
    elif isinstance(obj, list):
        for v in obj:
            out = find_key(v, key)
            if out:
                return out
    return ""

with open(sys.argv[1]) as fh:
    data = json.load(fh)

print(find_key(data, "reference_genome"))
PY
)"

    if [ -n "$JSON_REF" ]; then
        exec "$DIR/align.sh" --reference_genome "$JSON_REF"
    fi

    return 1
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
        dispatch_from_parameters_json
        ;;
    *)
        echo "ERROR: unknown task: $TASK" >&2
        exit 2
        ;;
esac

echo "ERROR: cannot determine module type from arguments" >&2
echo "ARGS: $*" >&2
exit 2