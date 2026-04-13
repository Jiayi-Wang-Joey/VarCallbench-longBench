#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    res <- list(
        output_dir = NULL,
        csv_files = character()
    )
    
    i <- 1
    while (i <= length(args)) {
        key <- args[i]
        
        if (key == "--output_dir") {
            res$output_dir <- args[i + 1]
            i <- i + 2
            
        } else if (key %in% c(
            "--alignment_qc.csv",
            "--alignment_qc_csv",
            "--alignment-qc-csv"
        )) {
            i <- i + 1
            vals <- character()
            
            while (i <= length(args) && !startsWith(args[i], "--")) {
                vals <- c(vals, args[i])
                i <- i + 1
            }
            
            res$csv_files <- c(res$csv_files, vals)
            
        } else {
            i <- i + 1
        }
    }
    
    res
}

args <- parse_args(args)

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

dt_list <- lapply(args$csv_files, fread)
dt <- rbindlist(dt_list, fill = TRUE)

fwrite(dt, file.path(args$output_dir, "alignment_qc_merged.csv"))