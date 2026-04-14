#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

argv <- commandArgs(trailingOnly = TRUE)

parse_args <- function(argv) {
    res <- list(
        output_dir = NULL,
        vcf_files = character(),
        task = NULL,
        name = NULL
    )
    
    i <- 1
    while (i <= length(argv)) {
        key <- argv[i]
        
        if (key == "--output_dir") {
            if (i == length(argv)) stop("--output_dir requires a value")
            res$output_dir <- argv[i + 1]
            i <- i + 2
            
        } else if (key == "--task") {
            if (i == length(argv)) stop("--task requires a value")
            res$task <- argv[i + 1]
            i <- i + 2
            
        } else if (key == "--name") {
            if (i == length(argv)) stop("--name requires a value")
            res$name <- argv[i + 1]
            i <- i + 2
            
        } else if (key %in% c(
            "--filtered.vcf", "--filtered_vcf", "--filtered-vcf",
            "--variant.vcf", "--variant_vcf", "--variant-vcf"
        )) {
            i <- i + 1
            vals <- character()
            
            while (i <= length(argv) && !startsWith(argv[i], "--")) {
                vals <- c(vals, argv[i])
                i <- i + 1
            }
            
            res$vcf_files <- c(res$vcf_files, vals)
            
        } else {
            stop("Unknown argument: ", key)
        }
    }
    
    res
}

extract_meta <- function(path) {
    bn <- basename(path)
    dataset <- sub("\\.vcf(\\.gz)?$", "", bn)
    
    caller <- sub(".*/variant_call/([^/]+)/.*", "\\1", path)
    if (identical(caller, path)) {
        caller <- NA_character_
    }
    
    tech <- fcase(
        grepl("bulk_ONT$", dataset), "cDNAxR10",
        grepl("dRNA_ONT$", dataset), "dRNA004",
        grepl("bulk_PB$", dataset), "Kinnex",
        default = NA_character_
    )
    
    cell_line <- sub("_.*$", "", dataset)
    
    data.table(
        file = path,
        dataset = dataset,
        caller = caller,
        cell_line = cell_line,
        tech = tech
    )
}

read_vcf_ids <- function(vcf_file) {
    cmd <- sprintf(
        "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' %s",
        shQuote(vcf_file)
    )
    
    dt <- fread(
        cmd = cmd,
        sep = "\t",
        header = FALSE,
        col.names = c("CHROM", "POS", "REF", "ALT")
    )
    
    unique(paste(dt$CHROM, dt$POS, dt$REF, dt$ALT, sep = "_"))
}

calc_overlap <- function(v_cdna, v_drna, v_pb) {
    all_ids <- unique(c(v_cdna, v_drna, v_pb))
    
    dt <- data.table(id = all_ids)
    dt[, cdna := id %in% v_cdna]
    dt[, drna := id %in% v_drna]
    dt[, pb := id %in% v_pb]
    dt[, pattern := paste0(as.integer(cdna), as.integer(drna), as.integer(pb))]
    
    out <- data.table(
        pattern = c("111", "110", "101", "011", "100", "010", "001"),
        category = c(
            "shared_all",
            "shared_cdna_drna",
            "shared_cdna_pb",
            "shared_drna_pb",
            "unique_cdna",
            "unique_drna",
            "unique_pb"
        )
    )
    
    counts <- dt[, .N, by = pattern]
    out <- merge(out, counts, by = "pattern", all.x = TRUE)
    out[is.na(N), N := 0L]
    
    total_union <- sum(out$N)
    out[, pct_union := if (total_union > 0) N / total_union else NA_real_]
    
    out
}

main <- function() {
    opt <- parse_args(argv)
    
    if (is.null(opt$output_dir) || !length(opt$vcf_files)) {
        stop("Need --output_dir and at least one VCF input")
    }
    
    dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
    
    meta <- rbindlist(lapply(opt$vcf_files, extract_meta), fill = TRUE)
    
    if (any(is.na(meta$caller))) {
        stop("Failed to parse caller from some VCF paths")
    }
    if (any(is.na(meta$tech))) {
        stop("Failed to parse tech from some dataset names")
    }
    
    id_list <- lapply(opt$vcf_files, read_vcf_ids)
    names(id_list) <- opt$vcf_files
    
    res <- list()
    groups <- unique(meta[, .(caller, cell_line)])
    
    for (i in seq_len(nrow(groups))) {
        caller_i <- groups$caller[i]
        cell_i <- groups$cell_line[i]
        
        subm <- meta[caller == caller_i & cell_line == cell_i]
        expected <- c("Kinnex", "cDNAxR10", "dRNA004")
        
        if (!all(expected %in% subm$tech)) {
            warning(sprintf(
                "Skipping caller=%s, cell_line=%s; found techs: %s",
                caller_i, cell_i, paste(sort(unique(subm$tech)), collapse = ",")
            ))
            next
        }
        
        f_cdna <- subm[tech == "cDNAxR10", file][1]
        f_drna <- subm[tech == "dRNA004", file][1]
        f_pb <- subm[tech == "Kinnex", file][1]
        
        ov <- calc_overlap(
            v_cdna = id_list[[f_cdna]],
            v_drna = id_list[[f_drna]],
            v_pb = id_list[[f_pb]]
        )
        
        ov[, `:=`(
            caller = caller_i,
            cell_line = cell_i,
            n_cdna = length(id_list[[f_cdna]]),
            n_drna = length(id_list[[f_drna]]),
            n_pb = length(id_list[[f_pb]])
        )]
        
        res[[length(res) + 1]] <- ov
    }
    
    if (length(res) == 0) {
        stop("No valid caller/cell_line groups produced results")
    }
    
    final_dt <- rbindlist(res, fill = TRUE)
    
    out_file <- file.path(opt$output_dir, "cross_platform_overlap.csv")
    fwrite(final_dt, out_file)
    
    if (!file.exists(out_file)) {
        stop("Failed to create output file: ", out_file)
    }
    
    message("Wrote: ", out_file)
}

main()