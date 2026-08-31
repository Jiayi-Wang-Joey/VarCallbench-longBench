# VarCallbench-longBench

Module repository for [VarCallbench](https://github.com/Jiayi-Wang-Joey/VarCallBench), an omnibenchmark for long-read variant calling and phasing.

Each stage in `benchmark.yaml` (in VarCallbench) dispatches into this repo via `run.sh --task <task>`, which routes to the corresponding script under `align/`, `variant_call/`, `gatk/`, or `eval/`.

## Layout

- `align/` — read alignment and alignment QC (minimap2, WGS alignment)
- `variant_call/` — variant callers (Clair3-RNA, DeepVariant, longcallR, longcallR-nn, isoLASER)
- `gatk/` — GATK-based WGS ground truth generation (MarkDuplicates, BQSR, HaplotypeCaller, etc.)
- `eval/` — evaluation and metric collection (hap.py, coverage BED generation)
- `data/` — raw data download scripts (WGS, rawdata)
- `arch/` — archived/unused scripts
