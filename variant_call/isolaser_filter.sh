#!/usr/bin/env bash
set -euo pipefail

DATASET=""
VCF=""
OUTDIR=""

echo "ARGS: $@" >&2

while [ $# -gt 0 ]; do
    case "$1" in
        --name|--dataset_id)
            DATASET="$2"; shift 2 ;;
        --isolaser.vcf|--isolaser_vcf)
            VCF="$2"; shift 2 ;;
        --output_dir)
            OUTDIR="$2"; shift 2 ;;
        --task)
            shift 2 ;;
        *)
            echo "Unknown arg: $1" >&2
            shift ;;
    esac
done

mkdir -p "$OUTDIR"

# isoLASER-specific post-filtering, per the author's own evaluation criteria
# (allelic ratio > 25%, GQ > 25, and collapsing multiple candidate alleles at
# the same position -- e.g. multiple INDEL alleles, or the same allele reported
# once per overlapping gene annotation -- down to the 2 highest-QUAL records).
# Deliberately its own stage (not folded into isolaser_run.sh) so that changing
# these filter thresholds only reruns this cheap step and everything downstream,
# not the expensive isoLASER calling step itself.

UNSORTED="$OUTDIR/tmp.unsorted.vcf"

bcftools view -e 'GT="0/0" || FORMAT/AB<=0.25 || FORMAT/GQ<=25' "$VCF" \
    | awk '
        BEGIN { FS = OFS = "\t" }
        /^#/ { print; next }
        {
            key = $1 SUBSEP $2
            if (key != prevkey && prevkey != "") flush()
            n[key]++
            idx = n[key]
            q = ($6 == "." || $6 == "nan") ? 0 : $6 + 0
            qual[key, idx] = q
            line[key, idx] = $0
            prevkey = key
        }
        END { if (prevkey != "") flush() }
        function flush(   i, j, tq, tl, cnt, top) {
            cnt = n[prevkey]
            for (i = 1; i <= cnt; i++) {
                for (j = i + 1; j <= cnt; j++) {
                    if (qual[prevkey, j] > qual[prevkey, i]) {
                        tq = qual[prevkey, i]; qual[prevkey, i] = qual[prevkey, j]; qual[prevkey, j] = tq
                        tl = line[prevkey, i]; line[prevkey, i] = line[prevkey, j]; line[prevkey, j] = tl
                    }
                }
            }
            top = (cnt < 2) ? cnt : 2
            for (i = 1; i <= top; i++) print line[prevkey, i]
            for (i = 1; i <= cnt; i++) { delete qual[prevkey, i]; delete line[prevkey, i] }
            delete n[prevkey]
        }' \
    > "$UNSORTED"

bcftools sort -T "$OUTDIR/tmp.sort" -Oz -o "$OUTDIR/${DATASET}.vcf.gz" "$UNSORTED"
bcftools index -t -f "$OUTDIR/${DATASET}.vcf.gz"
rm -f "$UNSORTED"
