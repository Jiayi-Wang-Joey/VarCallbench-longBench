#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(VariantAnnotation)
    library(ComplexUpset)
    library(ggplot2)
    library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    res <- list(
        output_dir = NULL,
        vcf_files = character()
    )
    
    i <- 1
    while (i <= length(args)) {
        key <- args[i]
        
        if (key == "--output_dir") {
            res$output_dir <- args[i + 1]
            i <- i + 2
            
        } else if (key %in% c("--variant.vcf", "--variant_vcf", "--variant-vcf")) {
            i <- i + 1
            vals <- character()
            
            while (i <= length(args) && !startsWith(args[i], "--")) {
                vals <- c(vals, args[i])
                i <- i + 1
            }
            
            res$vcf_files <- c(res$vcf_files, vals)
            
        } else {
            i <- i + 1
        }
    }
    
    res$vcf_files <- unique(res$vcf_files)
    res
}

extract_dataset <- function(vcf_path) {
    x <- basename(vcf_path)
    x <- sub("\\.vcf\\.gz$", "", x)
    x <- sub("\\.vcf$", "", x)
    x
}

extract_caller <- function(vcf_path) {
    p <- normalizePath(vcf_path, mustWork = FALSE)
    
    m <- regexpr("variant_call/([^/]+)/", p, perl = TRUE)
    if (m[1] != -1) {
        hit <- regmatches(p, m)
        caller <- sub("^variant_call/", "", hit)
        caller <- sub("/$", "", caller)
        return(caller)
    }
    
    parts <- strsplit(p, "/", fixed = TRUE)[[1]]
    
    known <- c("clair3_rna", "deep_variant", "longcallR", "longcallR_nn",
               "clair3-rna", "deepvariant", "longcallR-nn")
    hit <- parts[parts %in% known]
    if (length(hit) > 0) {
        return(hit[1])
    }
    
    if (length(parts) >= 2) {
        return(parts[length(parts) - 1])
    }
    
    "unknown"
}

clean_caller_name <- function(x) {
    map <- c(
        "clair3_rna"   = "Clair3-RNA",
        "deep_variant" = "DeepVariant",
        "longcallR"    = "longcallR",
        "longcallR_nn" = "longcallR-nn",
        "clair3-rna"   = "Clair3-RNA",
        "deepvariant"  = "DeepVariant",
        "longcallR-nn" = "longcallR-nn"
    )
    
    ifelse(x %in% names(map), unname(map[x]), x)
}

read_vcf_ids <- function(vcf_path) {
    vcf <- readVcf(vcf_path)
    
    rr <- rowRanges(vcf)
    
    ref <- as.character(ref(vcf))
    alt <- vapply(alt(vcf), function(a) paste(as.character(a), collapse = ","), character(1))
    
    data.table(
        seqnames = as.character(seqnames(rr)),
        pos = start(rr),
        ref = ref,
        alt = alt
    )[
        ,
        variant_id := paste(seqnames, pos, ref, alt, sep = ":")
    ][
        ,
        .(variant_id)
    ][
        !is.na(variant_id)
    ][
        unique(variant_id)
    ]
}

build_incidence <- function(meta_dt) {
    long_dt <- rbindlist(
        lapply(seq_len(nrow(meta_dt)), function(i) {
            ids <- read_vcf_ids(meta_dt$vcf[i])
            ids[, `:=`(
                dataset = meta_dt$dataset[i],
                caller = meta_dt$caller[i]
            )]
            ids
        }),
        use.names = TRUE,
        fill = TRUE
    )
    
    long_dt <- unique(long_dt, by = c("dataset", "caller", "variant_id"))
    
    dcast(
        long_dt,
        dataset + variant_id ~ caller,
        value.var = "caller",
        fun.aggregate = length
    )
}

make_upset_plot <- function(dt_sub, dataset_name) {
    caller_cols <- setdiff(colnames(dt_sub), c("dataset", "variant_id"))
    
    if (length(caller_cols) < 2) {
        return(
            ggplot() +
                theme_void() +
                annotate(
                    "text",
                    x = 0.5, y = 0.5,
                    label = paste0("Dataset: ", dataset_name, "\nNot enough callers to build UpSet plot")
                )
        )
    }
    
    plot_dt <- copy(dt_sub)[, c("dataset", "variant_id") := NULL]
    for (cc in caller_cols) {
        set(plot_dt, j = cc, value = as.integer(plot_dt[[cc]] > 0))
    }
    
    ComplexUpset::upset(
        plot_dt,
        intersect = caller_cols,
        name = "Variants",
        width_ratio = 0.18,
        min_size = 1
    ) +
        ggtitle(dataset_name) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
}

main <- function() {
    parsed <- parse_args(args)
    
    if (is.null(parsed$output_dir) || parsed$output_dir == "") {
        stop("Missing --output_dir")
    }
    
    if (length(parsed$vcf_files) == 0) {
        stop("No VCF files found from --variant.vcf")
    }
    
    dir.create(parsed$output_dir, recursive = TRUE, showWarnings = FALSE)
    
    meta_dt <- data.table(
        vcf = parsed$vcf_files,
        dataset = vapply(parsed$vcf_files, extract_dataset, character(1)),
        caller_raw = vapply(parsed$vcf_files, extract_caller, character(1))
    )[
        ,
        caller := clean_caller_name(caller_raw)
    ]
    
    incidence_dt <- build_incidence(meta_dt)
    
    datasets <- unique(incidence_dt$dataset)
    pdf_file <- file.path(parsed$output_dir, "upset.pdf")
    
    pdf(pdf_file, width = 10, height = 7)
    
    for (ds in datasets) {
        dt_sub <- incidence_dt[dataset == ds]
        p <- make_upset_plot(dt_sub, ds)
        print(p)
    }
    
    dev.off()
    
    message("Wrote: ", pdf_file)
}

main()