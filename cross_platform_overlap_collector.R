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
        col.names = c("CHROM", "POS", "REF", "ALT"),
        showProgress = FALSE
    )
    
    if (nrow(dt) == 0L) {
        return(character())
    }
    
    unique(paste(dt$CHROM, dt$POS, dt$REF, dt$ALT, sep = "_"))
}

calc_overlap_fast <- function(v_cdna, v_drna, v_pb) {
    cdna_drna <- intersect(v_cdna, v_drna)
    cdna_pb   <- intersect(v_cdna, v_pb)
    drna_pb   <- intersect(v_drna, v_pb)
    shared_all <- Reduce(intersect, list(v_cdna, v_drna, v_pb))
    
    n_111 <- length(shared_all)
    n_110 <- length(cdna_drna) - n_111
    n_101 <- length(cdna_pb)   - n_111
    n_011 <- length(drna_pb)   - n_111
    
    n_100 <- length(setdiff(v_cdna, union(v_drna, v_pb)))
    n_010 <- length(setdiff(v_drna, union(v_cdna, v_pb)))
    n_001 <- length(setdiff(v_pb,   union(v_cdna, v_drna)))
    
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
        ),
        N = c(n_111, n_110, n_101, n_011, n_100, n_010, n_001)
    )
    
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
    
    groups <- unique(meta[, .(caller, cell_line)])
    res <- vector("list", nrow(groups))
    res_i <- 0L
    
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
        f_pb   <- subm[tech == "Kinnex", file][1]
        
        v_cdna <- read_vcf_ids(f_cdna)
        v_drna <- read_vcf_ids(f_drna)
        v_pb   <- read_vcf_ids(f_pb)
        
        ov <- calc_overlap_fast(v_cdna, v_drna, v_pb)
        ov[, `:=`(
            caller = caller_i,
            cell_line = cell_i,
            n_cdna = length(v_cdna),
            n_drna = length(v_drna),
            n_pb = length(v_pb)
        )]
        
        res_i <- res_i + 1L
        res[[res_i]] <- ov
    }
    
    if (res_i == 0L) {
        stop("No valid caller/cell_line groups produced results")
    }
    
    final_dt <- rbindlist(res[seq_len(res_i)], fill = TRUE)
    
    out_file <- file.path(opt$output_dir, "cross_platform_overlap.csv")
    fwrite(final_dt, out_file)
    
    if (!file.exists(out_file)) {
        stop("Failed to create output file: ", out_file)
    }
    
    message("Wrote: ", out_file)
}

main()