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
            
        } else if (key %in% c("--somatic_detection.csv",
                              "--somatic_detection_csv",
                              "--somatic-detection-csv")) {
            i <- i + 1
            vals <- character()
            
            while (i <= length(args) && !startsWith(args[i], "--")) {
                vals <- c(vals, args[i])
                i <- i + 1
            }
            
            res$csv_files <- c(res$csv_files, vals)
            
        } else if (key %in% c("--name", "--task")) {
            i <- i + 2
            
        } else {
            stop("Unknown argument: ", key)
        }
    }
    
    res
}

opt <- parse_args(args)

if (is.null(opt$output_dir) || !length(opt$csv_files)) {
    stop("Need --output_dir and at least one --somatic_detection.csv input")
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

dt_list <- lapply(opt$csv_files, fread)
dt <- rbindlist(dt_list, fill = TRUE)

fwrite(dt, file.path(opt$output_dir, "somatic_detection_merged.csv"))