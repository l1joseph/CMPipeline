#!/usr/bin/env Rscript
# ==============================================================================
#  BATCH-CORRECTION & NORMALISATION PIPELINE
# ==============================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("pacman", quietly = TRUE))
      install.packages("pacman", repos = "https://cloud.r-project.org")
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
  pacman::p_load(tidyverse, vegan, ConQuR, argparse,
               zCompositions, devtools,
               doParallel, foreach, iterators) 
  ## DEICODE (for RCLR) --------------------------------------------------------
  if (!requireNamespace("DEICODE", quietly = TRUE)) {
    try(BiocManager::install("DEICODE", update = FALSE, ask = FALSE), silent = TRUE)
    if (!requireNamespace("DEICODE", quietly = TRUE))
      try(devtools::install_github("DEICODE-dev/DEICODE", quiet = TRUE), silent = TRUE)
  }
  if (requireNamespace("DEICODE", quietly = TRUE)) {
    library(DEICODE);  rclr_f <- DEICODE::rclr
  } else {
    message("WARN: DEICODE unavailable — using manual RCLR")
    rclr_f <- function(mat){
      t(apply(mat, 1, function(x){
        nz <- x>0; g <- exp(mean(log(x[nz]))); y <- rep(NA, length(x))
        y[nz] <- log(x[nz]/g); y
      }))
    }
  }
})

# ---------- small helpers -----------------------------------------------------
save_table <- function(df, path){
  if (!is.data.frame(df))             
      df <- as.data.frame(df)         
  rownames_to_column(df, "SampleID") %>%
    write_tsv(path, progress = FALSE)
}

load_table <- function(path)
  read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
  column_to_rownames("SampleID") %>% as.data.frame()

# ---------- PERMANOVA ---------------------------------------------------------
run_permanova <- function(mat, meta, rhs_formula, dist_meth){
  smp  <- intersect(rownames(mat), rownames(meta))
  if (anyNA(mat[smp, ])) mat[smp, ][is.na(mat[smp, ])] <- 0
  mat <- mat[ , colSums(mat) > 0, drop = FALSE]
  dist <- vegdist(mat[smp, ], dist_meth)
  res  <- adonis2(dist ~ ., data = meta[smp, all.vars(rhs_formula)],
                  permutations = 999, by = "margin")
  as_tibble(res, rownames = "Term") %>% dplyr::select(Term, Df, R2, `Pr(>F)`)
}

# ---------- ConQuR ------------------------------------------------------------
run_conqur <- function(tab, meta, batch_col, covars,
                       method = c("default", "lasso"), 
                       ref_batch = NULL,
                       lambda_quantile = "2p/n",
                       corr_dir = ".", 
                       verbose = TRUE) {
  method <- match.arg(method)
  if (!dir.exists(corr_dir)) {
    dir.create(corr_dir, recursive = TRUE)
  }
  
  cache_file <- file.path(corr_dir, paste0("ConQuR_", method, ".tsv"))
  if (file.exists(cache_file)) {
    if (verbose) message("SKIP  ▶  ", basename(cache_file), " (cached)")
    return(load_table(cache_file))
  }
                                                                   
  if (verbose) message("INFO: Running ConQuR with method = ", method)
  smp  <- intersect(rownames(tab), rownames(meta))
  if (length(smp) == 0) {
    stop("No common samples found between tab and meta")
  }
  tab  <- tab[smp, ];  meta <- meta[smp, ]
  if (is.null(ref_batch)){
    batch_counts <- table(meta[[batch_col]])
    ref_batch <- names(which.max(batch_counts))
    if (verbose) {
      message("INFO: Reference batch = ", ref_batch, 
              " (largest N = ", max(batch_counts), ")")
    }
  }
  if (!ref_batch %in% unique(meta[[batch_col]])) {
    stop("Reference batch '", ref_batch, "' not found in data")
  }
    
  ## ---- drop constant covariates ---------------------------------
  keep <- vapply(meta[, covars, drop = FALSE], function(v) {
    if (is.numeric(v)) var(v, na.rm = TRUE) > 0      
    else length(unique(v[!is.na(v)])) > 1            
  }, logical(1))
  covars <- covars[keep]   
  meta   <- droplevels(meta) 
  
  ## ---------------------------------------------------------------------
  if (method == "lasso") {
      res <- ConQuR(
        tax_tab        = tab,
        batchid        = factor(meta[[batch_col]]),
        covariates     = meta[, covars, drop = FALSE],
        batch_ref      = ref_batch,
        logistic_lasso = TRUE,
        quantile_type  = "lasso",
        interplt       = TRUE,  
        # lambda_quantile = lambda_quantile,
        num_core       = min(32, parallel::detectCores())
      )
    } else {
      res <- ConQuR(
        tax_tab        = tab,
        batchid        = factor(meta[[batch_col]]),
        covariates     = meta[, covars, drop = FALSE],
        batch_ref      = ref_batch,
        num_core       = min(32, parallel::detectCores())
      )
    }
    
  res <- as.data.frame(res)
  save_table(res, cache_file)                                         
  res
}

# ---------- PCoA plots --------------------------------------------------------
generate_conqur_plots <- function(tab, meta, batch_col, title, out_pref){
  smp <- intersect(rownames(tab), rownames(meta))
  fac <- meta[smp, batch_col]
  png(paste0(out_pref,"_bray.png"), 8,7,'in',res=300)
  Plot_PCoA(TAX = tab[smp, ], factor = fac,
            main = paste(title,", Bray-Curtis"))
  dev.off()
  png(paste0(out_pref,"_aitch.png"), 8,7,'in',res=300)
  Plot_PCoA(TAX = tab[smp, ], factor = fac, dissimilarity = "Aitch",
            main = paste(title,", Aitchison"))
  dev.off()
}

# ---------- normalise ---------------------------------------------------------
normalize_all <- function(mat, stub){
  dir.create(dirname(stub), showWarnings = FALSE, recursive = TRUE)
  paths <- c(pseudo = paste0(stub,"_clr_pseudo.tsv"),
             czm    = paste0(stub,"_clr_czm.tsv"),
             rclr   = paste0(stub,"_rclr.tsv"))
  out <- list()
  
  ## CLR + pseudo ------------------------------------------------------------
  if (file.exists(paths["pseudo"])) {
       out$CLR_pseudo <- load_table(paths["pseudo"])
  } else {
       pseudo <- min(mat[mat>0]) * 0.1
       out$CLR_pseudo <- decostand(mat + pseudo, "clr")
       save_table(out$CLR_pseudo, paths["pseudo"])
  }
  
  ## CLR + CZM ---------------------------------------------------------------
  if (file.exists(paths["czm"])) {
    out$CLR_CZM <- load_table(paths["czm"])
  } else {
    czm_mat <- cmultRepl(as.matrix(mat), method = "CZM", label = 0,
                         output = "p-counts",
                         z.warning = 1, z.delete  = FALSE,
                         suppress.print = TRUE)
    out$CLR_CZM <- decostand(czm_mat, "clr")   
    save_table(out$CLR_CZM, paths["czm"])        
  }
  
  ## RCLR --------------------------------------------------------------------
  if (file.exists(paths["rclr"])) {
    out$RCLR <- load_table(paths["rclr"])
  } else {
    out$RCLR <- rclr_f(as.matrix(mat)) %>%
                replace(is.na(.), 0)   
    save_table(out$RCLR, paths["rclr"])
  }
  out
}

# ==============================================================================
# SECTION 2  Main Workflow Execution
# ==============================================================================
parser <- ArgumentParser()
parser$add_argument("--otu",        required = TRUE, help = "Path to OTU table")
parser$add_argument("--meta",       required = TRUE, help = "Path to metadata")
parser$add_argument("--prefix",     required = TRUE, help = "Output prefix")
parser$add_argument("--batch_var",  required = TRUE, help = "Batch variable for correction")
parser$add_argument("--covariates", required = TRUE, help = "Comma-separated list of covariates")
parser$add_argument("--tumor_only", action = "store_true", 
                   help = "If set, filters for 'Tumor' samples in the 'Type' column of metadata.")

args <- parser$parse_args()

# Parse covariates
covar_list <- strsplit(args$covariates, ",")[[1]]
covar_list <- trimws(covar_list)  # Remove any whitespace

# Create output directories
dir.create(dirname(args$prefix), showWarnings = FALSE, recursive = TRUE)
corr_dir <- file.path(dirname(args$prefix), "corrected")   ; dir.create(corr_dir, showWarnings=FALSE)
norm_dir <- file.path(dirname(args$prefix), "normalized")  ; dir.create(norm_dir, showWarnings=FALSE)
plot_dir <- file.path(dirname(args$prefix), "pcoa_plots")  ; dir.create(plot_dir, showWarnings=FALSE)

# Load data
# Check file extension and read accordingly
if (grepl("\\.csv$", args$otu)) {
  otu_data <- read_csv(args$otu, show_col_types = FALSE)
} else {
  otu_data <- read_tsv(args$otu, show_col_types = FALSE)
}

# Check for appropriate column name
name_col <- intersect(c("species_name", "name", "species", "X"), colnames(otu_data))[1]
if (is.na(name_col)) {
  stop("ERROR: Could not find appropriate taxon name column in OTU table")
}

otu_raw <- otu_data %>%
           column_to_rownames(name_col) %>% 
           t() %>% 
           as.data.frame()

meta    <- read_tsv(args$meta, show_col_types = FALSE) %>%
           column_to_rownames("sampleid")

# Filter for tumor samples if the flag is provided
if (args$tumor_only) {
  message("INFO: Filtering for 'Tumor' samples only.")
  
  # Identify tumor sample IDs from metadata
  tumor_samples <- rownames(meta[meta$Type == "Tumor", ])
  
  # Filter both metadata and OTU table
  meta    <- meta[tumor_samples, ]
  otu_raw <- otu_raw[rownames(otu_raw) %in% tumor_samples, ]
  
  message("INFO: ", nrow(meta), " tumor samples remaining for analysis.")
}

# Prepare predictors (batch_var + covariates)
predictors  <- c(args$batch_var, covar_list)
rhs_formula <- reformulate(predictors)

## -----------------------  tables (raw + ConQuR) ----------------
tables <- list(raw = otu_raw)
tables$ConQuR_default <- run_conqur(
      otu_raw, meta, args$batch_var, covar_list, "default", corr_dir = corr_dir)

## -----------------------  permanova_summary ---------------------------------
sum_path <- file.path(dirname(args$prefix), "permanova_summary.tsv")
if (file.exists(sum_path)){
  message("SKIP  ▶  permanova_summary.tsv already exists → pipeline exit")
  quit(status = 0)
}

## -----------------------  main loop -----------------------------------------
permsum <- list()
for (tag in names(tables)){
  mat <- tables[[tag]]
  generate_conqur_plots(mat, meta, args$batch_var,
                        tag, file.path(plot_dir, tag))
  permsum[[length(permsum)+1]] <- run_permanova(
        mat, meta, rhs_formula, "bray") %>%
        mutate(Source = tag, Transform = "counts")
  
  norm_res <- normalize_all(mat,
        file.path(norm_dir, paste0(basename(args$prefix),"_",tag)))
  
  for (tr in names(norm_res)){
    permsum[[length(permsum)+1]] <- run_permanova(
          norm_res[[tr]], meta, rhs_formula, "euclidean") %>%
          mutate(Source = tag, Transform = tr)
  }
}

bind_rows(permsum) %>% save_table(sum_path)
message("DONE")
