#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(data.table)
})

argv <- commandArgs(trailingOnly = TRUE)

parse_args <- function(argv) {
    res <- list(output_dir = NULL, vcf_files = character())

    i <- 1
    while (i <= length(argv)) {
        key <- argv[i]
        if (key == "--output_dir") {
            res$output_dir <- argv[i + 1]; i <- i + 2
        } else if (key %in% c("--task", "--name")) {
            i <- i + 2
        } else if (key %in% c(
            "--filtered.vcf", "--filtered_vcf", "--filtered-vcf",
            "--variant.vcf",  "--variant_vcf",  "--variant-vcf"
        )) {
            i <- i + 1
            vals <- character()
            while (i <= length(argv) && !startsWith(argv[i], "--")) {
                vals <- c(vals, argv[i]); i <- i + 1
            }
            res$vcf_files <- c(res$vcf_files, vals)
        } else {
            stop("Unknown argument: ", key)
        }
    }
    res
}

extract_meta <- function(path) {
    bn      <- basename(path)
    dataset <- sub("\\.vcf(\\.gz)?$", "", bn)

    caller <- sub(".*/variant_call/([^/]+)/.*", "\\1", path)
    if (identical(caller, path)) caller <- NA_character_

    tech <- fcase(
        grepl("bulk_ONT", dataset), "cDNA",
        grepl("dRNA_ONT", dataset), "dRNA",
        grepl("bulk_PB",  dataset), "Kinnex",
        default = NA_character_
    )

    cell_line <- sub("_.*$", "", dataset)

    data.table(file = path, dataset = dataset, caller = caller,
               cell_line = cell_line, tech = tech)
}

vcf_to_ids <- function(path) {
    vcf <- VariantAnnotation::readVcf(path)
    if (length(vcf) == 0L) return(character())
    rr  <- rowRanges(vcf)
    alt <- vapply(VariantAnnotation::alt(vcf),
                  function(a) paste(as.character(a), collapse = ","),
                  character(1))
    unique(paste(as.character(seqnames(rr)), start(rr),
                 as.character(VariantAnnotation::ref(vcf)), alt, sep = "_"))
}

summarize_group <- function(ids_list, cell_line, caller) {
    techs   <- names(ids_list)
    all_ids <- unique(unlist(ids_list, use.names = FALSE))

    if (length(all_ids) == 0L) return(NULL)

    presence <- as.data.table(
        setNames(
            lapply(ids_list, function(v) as.integer(all_ids %in% v)),
            techs
        )
    )
    presence[, variant := all_ids]

    expected <- c("cDNA", "Kinnex", "dRNA")
    for (nm in expected) if (!nm %in% names(presence)) presence[, (nm) := 0L]

    presence[, n_tech := (cDNA > 0L) + (Kinnex > 0L) + (dRNA > 0L)]

    presence[, share_class := fcase(
        n_tech == 3L, "shared_by_3",
        n_tech == 2L, "shared_by_2",
        default      = "shared_by_1"
    )]

    presence[, detail_class := fcase(
        cDNA > 0L & Kinnex > 0L & dRNA > 0L, "cDNA_Kinnex_dRNA",
        cDNA > 0L & Kinnex > 0L,              "cDNA_Kinnex",
        cDNA > 0L & dRNA > 0L,                "cDNA_dRNA",
        Kinnex > 0L & dRNA > 0L,              "Kinnex_dRNA",
        cDNA > 0L,                             "cDNA_only",
        Kinnex > 0L,                           "Kinnex_only",
        default                              = "dRNA_only"
    )]

    presence[, `:=`(cell_line = cell_line, caller = caller)]

    variant_table <- presence[, .(
        cell_line, caller, variant,
        cDNA   = as.integer(cDNA   > 0L),
        Kinnex = as.integer(Kinnex > 0L),
        dRNA   = as.integer(dRNA   > 0L),
        n_tech, share_class, detail_class
    )]

    summary_table <- variant_table[, .N, by = .(cell_line, caller, share_class)]
    summary_table[, proportion := N / sum(N), by = .(cell_line, caller)]

    detail_table <- variant_table[, .N, by = .(cell_line, caller, detail_class)]
    detail_table[, proportion := N / sum(N), by = .(cell_line, caller)]

    list(variant = variant_table, summary = summary_table, detail = detail_table)
}

main <- function() {
    opt <- parse_args(argv)

    if (is.null(opt$output_dir) || !length(opt$vcf_files))
        stop("Need --output_dir and at least one VCF input")

    dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

    meta <- rbindlist(lapply(opt$vcf_files, extract_meta), fill = TRUE)

    # drop empty VCFs and unknown tech/caller
    meta <- meta[!is.na(tech) & !is.na(caller)]

    groups <- unique(meta[, .(caller, cell_line)])

    variant_list <- summary_list <- detail_list <- vector("list", nrow(groups))
    n <- 0L

    for (i in seq_len(nrow(groups))) {
        cl  <- groups$caller[i]
        cll <- groups$cell_line[i]
        sub <- meta[caller == cl & cell_line == cll]

        # need all three technologies; skip if any missing
        if (!all(c("cDNA", "Kinnex", "dRNA") %in% sub$tech)) {
            warning(sprintf("Skipping caller=%s cell_line=%s: missing techs (%s)",
                            cl, cll, paste(sort(sub$tech), collapse = ",")))
            next
        }

        ids_list <- setNames(
            lapply(c("cDNA", "Kinnex", "dRNA"), function(t) {
                vcf_to_ids(sub[tech == t, file][1])
            }),
            c("cDNA", "Kinnex", "dRNA")
        )

        res <- summarize_group(ids_list, cll, cl)
        if (is.null(res)) next

        n <- n + 1L
        variant_list[[n]] <- res$variant
        summary_list[[n]] <- res$summary
        detail_list[[n]]  <- res$detail
    }

    if (n == 0L) stop("No valid groups produced results")

    fwrite(rbindlist(variant_list[seq_len(n)]),
           file.path(opt$output_dir, "cross_platform_variant_table.csv"))
    fwrite(rbindlist(summary_list[seq_len(n)]),
           file.path(opt$output_dir, "cross_platform_summary.csv"))
    fwrite(rbindlist(detail_list[seq_len(n)]),
           file.path(opt$output_dir, "cross_platform_detail_summary.csv"))

    message("Wrote cross_platform_*.csv to ", opt$output_dir)
}

main()
