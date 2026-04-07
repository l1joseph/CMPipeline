# Pipeline Architecture

CMPipeline is a Nextflow DSL2 metagenomics pipeline for cancer microbiome analysis. It processes unmapped BAM files through host depletion, taxonomic classification, consensus taxa computation, statistical decontamination, and batch correction.

## Workflow Overview

The workflow in `main.nf` is controlled by flags and has 5 sequential steps. Entry points can be configured via `params.start_from` (`beginning`, `decontam`, `batch_correction`).

## Step 1: Host Depletion (`params.skip_host_depletion`)

`extractReads` (samtools) -> `filterReads` (fastp) -> `mapReads` (minimap2)

Reads are extracted from unmapped BAMs, quality-filtered, then sequentially mapped against three human reference databases to remove host reads:
1. **hg38** (GRCh38)
2. **T2T-CHM13v2.0 + PhiX**
3. **Human pangenome** (multiple .mmi files iterated)

Each step keeps only unmapped reads (samtools flag 4). R1 and R2 are processed independently through all three reference genomes. FastQC runs at multiple stages (pre-filter, post-filter, post-mapping per reference).

## Step 2: Taxonomic Classification (`params.skip_classification`)

Pangenome-depleted reads go to two classifiers in parallel:
- **KrakenUniq + Bracken** -> genus & species abundance (MPA format)
- **MetaPhlAn4** -> genus, species, SGB profiles
- **HUMAnN3** (optional, `params.skip_humann3`, default: true) -> functional profiling (gene families, pathway abundance/coverage). Reuses existing MetaPhlAn profiles when available. Per-sample outputs are merged via `merge_humann3` (produces combined tables + normalized CPM/relab). Output: `RESULTS/HUMANN3/` and `RESULTS/HUMANN3/merged/`.

Results are merged per-tool in `preprocess_taxa.nf` via `process_bracken` (combines `.G.mpa.krakenreport.txt` files) and `process_metaphlan` (merges `.profiled_metagenome.txt` files and extracts genus/species/SGB levels).

## Step 3: Consensus Taxa (`params.skip_consensus`)

`compute_consensus_taxa.py` intersects MetaPhlAn and Bracken genus/species lists, filtering to samples with >=100K reads, and retains genera detected by both tools at >=1% prevalence. Outputs a scatter plot (PDF) and filtered genus/species tables. When `skip_consensus=true`, Bracken output feeds directly to downstream steps. Default is `false` (consensus enabled).

## Step 4: Decontamination (`params.run_decontam`)

`251006_decontamination_ver2.R` -- the production decontamination script. Uses the `decontam` R package with `phyloseq` for prevalence-based contaminant identification, run per-batch. Key behavior:
- Auto-detects sample ID, type, and taxon name columns from various naming conventions
- Standardizes sample types to "Tumor" / "Control" using configurable `--tumor_values` / `--control_values` lists
- Expands metadata for samples in the OTU table that lack metadata entries (e.g., TCGA-style IDs with type codes)
- Runs `isContaminant(method="prevalence")` per batch, generates scatter plots per batch, and applies a 2-step rare taxa filter (sample-wise abundance + cohort-wide prevalence)
- Produces decontaminated OTU CSV and a summary text file

## Step 5: Batch Correction (`params.run_batch_correction`)

`2500703_batch_correction_normalization.r` -- auto-phased batch correction using ConQuR (Tune_ConQuR). Key behavior:
- **Phase system**: Phase 1 does a fast grid search (16 combinations, ~10-20 min). If PERMANOVA R2 reduction >= threshold (default 25%), stops. Otherwise runs Phase 3 (72 combinations, ~60-90 min)
- Selects optimal reference batch by Bray-Curtis centroid distance
- Validates covariates at batch level, drops problematic ones automatically
- Outputs: corrected count tables, multiple normalization formats (CLR, rCLR via DEICODE, relative abundance, rarefied), PCoA plots, PERMANOVA summary
- ConQuR is installed from GitHub (`ivartb/ConQuR_par` fork preferred for `batchid` fix in foreach workers)

## Module Pattern

Each Nextflow module follows a consistent pattern:
1. **Skip check**: If output files already exist in `publishDir`, symlink them and exit early
2. **Actual execution**: Run the bioinformatics tool
3. **Conda env**: Each process uses a dedicated conda env from `conda_envs/*.yml`
4. **Scratch**: Most processes use `scratch true` for local temp storage on compute nodes

When modifying modules, preserve the skip-check pattern -- it enables efficient resumption without `-resume`.

## Output Structure

Results are published to `RESULTS/` subdirectories:
- `UNMAPPED_BAM/` -- Extracted and filtered FASTQ files
- `MAPPED_READS/` -- Host-depleted FASTQ files (named `*.hg38.t2t.pangenome.fastq.gz`)
- `FASTQC/`, `MULTIQC/` -- QC reports
- `BRACKEN/` -- KrakenUniq/Bracken classification results and merged MPA reports
- `METAPHLAN4/` -- MetaPhlAn4 profiles and merged abundance tables (genus/species/SGB)
- `HUMANN3/` -- Per-sample functional profiles; `merged/` has combined + normalized tables (CPM, relab)
- `CONSENSUS_TAXA/` -- Intersected taxa tables and scatter plot
- `04_DECONTAMINATION/` -- Decontaminated OTU tables, scatter plots per batch, summary
- `05_BATCH_CORRECTION/` -- Corrected count tables (`corrected/`), normalized tables (`normalized/`), PCoA plots (`pcoa_plots/`), PERMANOVA summary
- `ANTISMASH/` -- AntiSMASH BGC results per sample (standalone workflow)
- `MAG/` -- nf-core/mag assembly + binning output (standalone workflow)
