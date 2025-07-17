#!/usr/bin/env Rscript

# ==============================================================================
# SECTION 1: Load Libraries and Define Helper Functions
# ==============================================================================

# Load necessary libraries
if (!require("pacman")) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(tidyverse, phyloseq, decontam, argparse, tools, ggrepel)

# --- Helper Function 1: Data Loading and Preprocessing ---
load_and_prepare_data <- function(otu_path, meta_path) {
  cat("INFO: Loading and preparing data...\n")
  
  # Read the raw OTU table and metadata from the provided paths
  otu_data_raw <- read.csv(otu_path, sep="\t", check.names = FALSE)
  metadata_raw <- read.csv(meta_path, sep="\t", check.names = FALSE)
  
  # --- OTU table formatting ---
  # Dynamically find the column to be used for row names (taxa names)
  name_col <- intersect(c("name", "species_name", "X", "species"), colnames(otu_data_raw))[1]
  if (is.na(name_col)) {
    stop("ERROR: A suitable column for taxon names was not found in the OTU table.")
  }
  rownames(otu_data_raw) <- otu_data_raw[[name_col]]
  
  # Isolate numeric columns and convert NAs to 0
  numeric_cols <- sapply(otu_data_raw, is.numeric)
  otu_data <- as.data.frame(otu_data_raw[, numeric_cols])
  otu_data[is.na(otu_data)] <- 0
  
  # --- Data Pruning and Quality Control ---
  # Prune empty samples
  original_sample_count <- ncol(otu_data)
  otu_data <- otu_data[, (colSums(otu_data, na.rm=TRUE) != 0)]
  cat("INFO: Pruned", original_sample_count - ncol(otu_data), "empty samples.\n")

  # Prune very rare taxa (singletons)
  original_taxa_count <- nrow(otu_data)
  taxa_sums <- rowSums(otu_data)
  taxa_to_keep <- taxa_sums > 1
  otu_data <- otu_data[taxa_to_keep, ]
  cat("INFO: Pruned", original_taxa_count - nrow(otu_data), "rare taxa (singletons).\n")

  # --- Match OTU table and metadata ---
  # Find and filter for common samples
  common_samples <- intersect(colnames(otu_data), metadata_raw$sampleid)
  if (length(common_samples) == 0) {
    stop("ERROR: No common samples found. Please check sample IDs in both files.")
  }
  
  otu_data <- otu_data[, common_samples]
  metadata <- metadata_raw[metadata_raw$sampleid %in% common_samples, ]
  rownames(metadata) <- metadata$sampleid
  
  cat("INFO:", length(common_samples), "common samples prepared for final analysis.\n")
  
  # Return the processed data as a list
  return(list(otu = otu_data, meta = metadata, name_col = name_col))
}

# --- Helper Function 2: Plotting Function ---
create_and_save_plot <- function(plot_data, blacklist, title, plot_path) {
  # Ensure the contaminant column exists
  if (!"is_contaminant" %in% colnames(plot_data)) {
    cat("WARN: 'is_contaminant' column not found for", title, ". Skipping plot.\n")
    return(NULL)
  }

  # Add a column to label species that are on the blacklist
  blacklist_flexible <- gsub(" ", "[ _]", blacklist)
  blacklist_pattern <- paste(blacklist_flexible, collapse = "|")
  
  # Add a column to label species that match any name in the blacklist.
  plot_data <- plot_data %>%
    mutate(label = ifelse(grepl(blacklist_pattern, name, ignore.case = TRUE), name, ""))

  p <- ggplot(plot_data, aes(x = prevalence_control, y = prevalence_sample)) +
    geom_point(aes(color = is_contaminant), alpha = 0.7, size = 2) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                    box.padding = 0.5, segment.color = 'grey50') +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50"),
                       name = "Contaminant",
                       labels = c("TRUE" = "Yes", "FALSE" = "No")) +
    labs(
      title = title,
      subtitle = paste(sum(plot_data$is_contaminant), "contaminants identified"),
      x = "Prevalence in Control Group (Blood)",
      y = "Prevalence in Sample Group (Tumor)"
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "red"),
      legend.position = "bottom"
    )
  
  # Save the plot
  ggsave(plot_path, plot = p, width = 8, height = 7, dpi = 300)
  cat("INFO: Saved plot to", plot_path, "\n")
}

# --- Helper Function 3: Decontam Package Identification ---
identify_contaminants_decontam <- function(otu_df, meta_df, blacklist_species, plot_dir_path, threshold, batch_var) {
  cat("INFO: Identifying contaminants using decontam package...\n")
  cat("INFO: Using batch variable:", batch_var, "\n")
  
  # Check if batch variable exists in metadata
  if (!batch_var %in% colnames(meta_df)) {
    stop(paste("ERROR: Batch variable", batch_var, "not found in metadata columns."))
  }
  
  master_table <- data.frame(Contaminant = character(0), Batch_Count = character(0))

  for(batch_id in unique(meta_df[[batch_var]])) {
    meta_batch <- meta_df[meta_df[[batch_var]] == batch_id, ]
    
    if(nrow(meta_batch) <= 2 || n_distinct(meta_batch$Type) < 2) {
      cat("WARN: [decontam] Skipping batch", batch_id, "- not enough samples for decontam.\n")
      next
    }
    
    otu_batch <- otu_df[, meta_batch$sampleid]
    meta_batch$is.neg <- meta_batch$Type == "Normal"
    if(sum(meta_batch$is.neg, na.rm=TRUE) == 0){
      cat("WARN: [decontam] Skipping batch", batch_id, "- no negative controls found.\n")
      next
    }
    
    physeq_batch <- phyloseq(otu_table(as.matrix(otu_batch), taxa_are_rows = TRUE), sample_data(meta_batch))
    contam_results <- isContaminant(physeq_batch, method="prevalence", neg="is.neg", threshold=threshold, detailed = TRUE)
    
    # Create data for plotting
    ps_pa <- transform_sample_counts(physeq_batch, function(abund) 1*(abund>0))
    ps_pa_neg <- prune_samples(sample_data(ps_pa)$is.neg, ps_pa)
    ps_pa_pos <- prune_samples(!sample_data(ps_pa)$is.neg, ps_pa)

    plot_data_decontam <- data.frame(
      name = taxa_names(ps_pa),
      prevalence_control = (taxa_sums(ps_pa_neg) / nsamples(ps_pa_neg)),
      prevalence_sample = (taxa_sums(ps_pa_pos) / nsamples(ps_pa_pos)),
      is_contaminant = contam_results$contaminant
    )

    # Generate and save the plot
    plot_title <- paste("Decontam Package - Batch:", batch_id)
    plot_path <- file.path(plot_dir_path, paste0("decontam_pkg_", batch_id, ".png"))
    create_and_save_plot(plot_data_decontam, blacklist_species, plot_title, plot_path)
    
    # Identify and store contaminants
    contaminants_in_batch <- rownames(contam_results)[contam_results$contaminant == TRUE]
    
    if(length(contaminants_in_batch) > 0){
        master_table <- rbind(master_table, data.frame(Contaminant = contaminants_in_batch, Batch_Count = batch_id))
    }
  }
  return(master_table)
}

# --- Helper Function 4: Filter Rare OTUs by Prevalence ---
filter_rare_otu <- function(otu_table, prevalence_threshold = 0.05) {
  cat("\nINFO: Filtering rare OTUs with a prevalence threshold of", prevalence_threshold * 100, "%...\n")

  # Ensure the table is numeric
  otu_table_numeric <- otu_table %>% mutate(across(everything(), as.numeric))

  # Calculate the minimum number of samples an OTU must be present in
  total_samples <- ncol(otu_table_numeric)
  min_sample_count <- total_samples * prevalence_threshold

  # Find which OTUs (rows) meet the prevalence threshold
  presence_counts <- rowSums(otu_table_numeric > 0)
  otus_to_keep <- presence_counts >= min_sample_count

  # Filter the original table
  filtered_otu_table <- otu_table[otus_to_keep, ]

  cat("INFO: Removed", nrow(otu_table) - nrow(filtered_otu_table), "rare OTUs. Kept", nrow(filtered_otu_table), "OTUs.\n")

  return(filtered_otu_table)
}

# --- Helper Function 5: Contaminant Removal ---
remove_contaminants <- function(original_otu, master_contaminant_table, meta_df, batch_var) {
  clean_otu <- original_otu
  if (nrow(master_contaminant_table) > 0) {
    cat("INFO: Removing", nrow(master_contaminant_table), "contaminant entries across all batches.\n")
    for (i in 1:nrow(master_contaminant_table)) {
      contaminant_name <- master_contaminant_table$Contaminant[i]
      batch_val <- master_contaminant_table$Batch_Count[i]
      samples_in_batch <- meta_df$sampleid[meta_df[[batch_var]] == batch_val]
      if (contaminant_name %in% rownames(clean_otu)) {
          clean_otu[contaminant_name, samples_in_batch] <- 0
      }
    }
  } else {
    cat("INFO: No contaminants to remove.\n")
  }
  return(clean_otu)
}

# --- Helper Function 6: Save Decontaminated Table ---
save_decontaminated_table <- function(clean_otu_table, original_name_col, file_prefix, method_suffix) {
  output_file_path <- paste0(file_prefix, ".", method_suffix, "_decontaminated.csv")
  cat("INFO: Writing", method_suffix, "decontaminated OTU table to:", output_file_path, "\n")
  clean_otu_table %>%
    rownames_to_column(var = original_name_col) %>%
    write.csv(., output_file_path, row.names = FALSE)
}

# --- Helper Function 7: Plot Contaminant Proportions ---
plot_contaminant_proportions <- function(original_otu, master_contaminant_table, meta_df, title, plot_path, batch_var) {
  cat("\nINFO: Generating contaminant proportion summary plot for:", title, "...\n")

  # Return early if no contaminants were identified
  if (nrow(master_contaminant_table) == 0) {
    cat("INFO: No contaminants were identified. Skipping proportion plot.\n")
    return(NULL)
  }

  # 1. Calculate total reads per sample
  total_reads <- colSums(original_otu)

  # 2. Calculate contaminant reads per sample, respecting the batch-specific nature
  otu_long_meta <- original_otu %>%
    rownames_to_column(var = "name") %>%
    pivot_longer(cols = -name, names_to = "sampleid", values_to = "reads") %>%
    inner_join(dplyr::select(meta_df, sampleid, all_of(batch_var)), by = "sampleid")

  contaminant_reads_df <- otu_long_meta %>%
    inner_join(master_contaminant_table, by = c("name" = "Contaminant", batch_var = "Batch_Count")) %>%
    group_by(sampleid) %>%
    summarise(contaminant_reads = sum(reads, na.rm = TRUE), .groups = 'drop')

  # 3. Create the final summary dataframe, aggregated by Cohort Type
  summary_df <- data.frame(
    sampleid = names(total_reads),
    total_reads = total_reads
  ) %>%
    left_join(contaminant_reads_df, by = "sampleid") %>%
    mutate(
      contaminant_reads = ifelse(is.na(contaminant_reads), 0, contaminant_reads),
      proportion = ifelse(total_reads > 0, contaminant_reads / total_reads, 0)
    ) %>%
    inner_join(dplyr::select(meta_df, sampleid, Type), by = "sampleid") %>%
    mutate(Type = ifelse(tolower(Type) == "tumor", "Tumor", "Blood"))

  # 4. Create the box plot with jittered points
  p <- ggplot(summary_df, aes(x = Type, y = proportion, fill = Type)) +
    geom_boxplot(alpha = 0.4, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.25, alpha = 0.5, size = 2) +
    stat_summary(fun = mean, geom = "text", aes(label = sprintf("%.3f", after_stat(y))), vjust = -2, size = 4) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = title,
      subtitle = "Proportion of reads from identified contaminants per sample",
      x = "Cohort",
      y = "Contaminant Fraction"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  # 5. Save the plot
  ggsave(plot_path, plot = p, width = 8, height = 7, dpi = 300)
  cat("INFO: Saved contaminant proportion plot to", plot_path, "\n")
}

# ==============================================================================
# SECTION 2: Main Workflow Execution
# ==============================================================================

# --- Define Blacklist ---
blacklist_species <- c("Rhodococcus fascians", "Streptococcus sanguinis", "Kocuria palustris", 
                      "Pseudomonas putida", "Psychrobacter pulmonis", "Geobacillus vulcani", 
                      "Bosea lupini", "Ensifer adhaerens", "Novosphingobium pentaromativorans", 
                      "Methylobacterium aminovorans", "Brevundimonas aurantiaca", 
                      "Caulobacter vibrioides", "Psychrobacter maritimus", "Sphingomonas faeni", 
                      "Sphingomonas koreensis", "Acinetobacter schindleri", "Ensifer meliloti", 
                      "Algoriphagus aquaeductus", "Shinella zoogloeoides", "Sphingopyxis macrogoltabida", 
                      "Pseudomonas fragi", "Ochrobactrum tritici", "Phenylobacterium haematophilum", 
                      "Comamonas aquatica", "Acinetobacter baumannii", "Acinetobacter calcoaceticus", 
                      "Sphingopyxis alaskensis", "Bradyrhizobium liaoningense", 
                      "Flavobacterium lindanitolerans", "Ralstonia insidiosa", "Pseudomonas stutzeri", 
                      "Acinetobacter junii", "Microbacterium oxydans", "Microbacterium chocolatum", 
                      "Brevibacterium epidermidis", "Ochrobactrum intermedium", "Lactococcus lactis", 
                      "Staphylococcus hominis", "Shinella granuli", "Streptococcus salivarius", 
                      "Sphingomonas mucosissima", "Methylobacterium dankookense", 
                      "Haemophilus parainfluenzae", "Sphingobacterium spiritivorum", 
                      "Sphingomonas leidyi", "Caulobacter leidyi", "Acinetobacter parvus", 
                      "Methylobacterium jeotgali", "Methylobacterium adhaesivum", 
                      "Corynebacterium kroppenstedtii", "Pseudomonas brenneri", 
                      "Acinetobacter towneri", "Pelomonas aquatica", "Methylobacterium oryzae", 
                      "Acinetobacter lwoffii", "Methylobacterium oxalidis", "Streptococcus mitis", 
                      "Corynebacterium tuberculostearicum", "Kocuria rhizophila", 
                      "Acidovorax defluvii", "Pseudomonas veronii", "Delftia acidovorans", 
                      "Afipia broomeae", "Pedomicrobium australicum", "Rhizobium radiobacter", 
                      "Bosea vestrisii", "Bradyrhizobium daqingense", "Staphylococcus epidermidis", 
                      "Sphingomonas yabuuchiae", "Delftia tsuruhatensis", "Brevundimonas vesicularis", 
                      "Bradyrhizobium denitrificans", "Escherichia flexneri", "Shigella flexneri", 
                      "Halomonas axialensis", "Methylobacterium radiotolerans", 
                      "Enhydrobacter aerosaccus", "Ralstonia pickettii", "Halomonas aquamarina", 
                      "Micrococcus luteus", "Paracoccus aminovorans", "Halomonas sulfidaeris", 
                      "Bradyrhizobium japonicum", "Pseudomonas fluorescens", "Halomonas meridiana", 
                      "Shewanella algae", "Brevundimonas diminuta", "Acinetobacter johnsonii", 
                      "Bradyrhizobium elkanii", "Sphingomonas echinoides", "Acinetobacter guillouiae", 
                      "Stenotrophomonas maltophilia", "Variovorax paradoxus", 
                      "Propionibacterium acnes", "Pelomonas puraquae")

# Species in this list will be removed from the OTU table before any analysis
hardcoded_removal_list <- c("Homo sapiens", "PhiX")

# --- Argument Parsing ---
parser <- ArgumentParser(description="Perform manual and decontam-based decontamination on an OTU table.")
parser$add_argument("--otu_table", type="character", required=TRUE, help="Path to the input OTU table file.")
parser$add_argument("--metadata", type="character", required=TRUE, help="Path to the metadata file.")
parser$add_argument("--prefix", type="character", required=TRUE, help="Prefix for the output files.")
parser$add_argument("--threshold", type="double", default=0.1, help="Threshold for decontam package [default: 0.1]")
parser$add_argument("--min_prevalence", type="double", default=0.05, help="Minimum prevalence threshold to keep an OTU [default: 0.05]")
parser$add_argument("--var_batch", type="character", required=TRUE, help="Variable for batch correction")
args <- parser$parse_args()

# --- Create Directories ---
output_dir <- dirname(args$prefix)
plots_dir <- file.path(output_dir, "decontamination_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# --- Initial Data Preparation ---
prepared_data <- load_and_prepare_data(args$otu_table, args$metadata)
otu_original <- prepared_data$otu
meta_original <- prepared_data$meta
original_taxon_col_name <- prepared_data$name_col

# --- Pre-filtering Step: Hard-coded Removal ---
cat("\nINFO: Applying hard-coded removal filter...\n")
initial_taxa_count <- nrow(otu_original)
removal_pattern <- paste(hardcoded_removal_list, collapse = "|")
otu_pre_filtered <- otu_original[!grepl(removal_pattern, rownames(otu_original), ignore.case = TRUE), ]
cat("INFO: Removed", initial_taxa_count - nrow(otu_pre_filtered), "species based on the hard-coded list.\n")

# --- Run Decontamination Method : Using 'decontam' Package ---
otu_clean <- filter_rare_otu(otu_pre_filtered, args$min_prevalence)
contaminants_decontam <- identify_contaminants_decontam(otu_clean, meta_original, blacklist_species, plots_dir, threshold = args$threshold, batch_var = args$var_batch)
otu_filtered <- remove_contaminants(otu_clean, contaminants_decontam, meta_original, batch_var = args$var_batch)
save_decontaminated_table(otu_filtered, original_taxon_col_name, args$prefix, "decontam_pkg")

# --- Generate Final Summary Plots ---
plot_contaminant_proportions(
  original_otu = otu_original,
  master_contaminant_table = contaminants_decontam,
  meta_df = meta_original,
  title = "Contaminant Read Proportion (Decontam Pkg Method)",
  plot_path = file.path(plots_dir, "summary_proportions_decontam_pkg.png"),
  batch_var = args$var_batch
)

# --- End of Script ---
cat("INFO: Decontamination script finished successfully.\n")
