# Scripts Reference

## Python Scripts (called by Nextflow modules)

| Script | Called By | Purpose |
|--------|-----------|---------|
| `scripts/kreport2mpa.py` | `Bracken` module | Converts KrakenUniq/Bracken kraken-style reports to MPA (MetaPhlAn) format. From KrakenTools. |
| `scripts/combine_mpa.py` | `process_bracken` | Merges multiple MPA-format reports into a single combined table with one column per sample. From KrakenTools. |
| `scripts/merge_metaphlan_tables.py` | `process_metaphlan` | Merges multiple MetaPhlAn profile files into a single abundance table. Validates all inputs are from the same MetaPhlAn version. |
| `scripts/compute_consensus_taxa.py` | `consensus_taxa` | Intersects MetaPhlAn and Bracken genus lists. Filters samples >=100K reads, retains genera at >=1% prevalence in both tools. Outputs filtered genus/species tables and a scatter plot PDF. Dependencies: pandas, matplotlib, seaborn, adjustText. |

## R Scripts (called by Nextflow modules)

| Script | Called By | Purpose |
|--------|-----------|---------|
| `scripts/251006_decontamination_ver2.R` | `Decontamination` module | **Production decontamination script** (~700 lines). Uses `decontam` + `phyloseq`. Auto-detects column naming conventions. Runs per-batch prevalence-based contaminant ID. Generates per-batch scatter plots. Applies 2-step rare taxa filtering. |
| `scripts/2500703_batch_correction_normalization.r` | `BatchCorrection` module | **Production batch correction script** (~800 lines). Auto-phased ConQuR with PERMANOVA stopping criteria. Selects reference batch by centroid distance. Multiple normalizations (CLR, rCLR, relative abundance, rarefied). PCoA visualization. |
| `scripts/decontamination.R` | Not used in pipeline | **Older/simpler decontamination script** (~260 lines). Runs two methods: manual prevalence-based (Fisher's exact test per batch) and `decontam` package. Outputs two decontaminated tables. Has hardcoded interactive-mode paths. |

Both production R scripts auto-detect column names for sample IDs, taxon names, and sample types from multiple naming conventions -- they handle TCGA-style IDs, various metadata column names, and Bracken report suffixes.

## Shell Scripts (setup utilities)

| Script | Purpose |
|--------|---------|
| `scripts/download_references.sh` | Downloads GRCh38, T2T-CHM13v2.0, PhiX, and all HPRC pangenome assemblies to `ref/` |
| `scripts/create_minimap2_indexes.sh` | Builds minimap2 `.mmi` indexes for hg38, T2T+PhiX, and each pangenome assembly |

## SLURM Batch Scripts (`scripts/batch_scripts/`)

| Script | Purpose | Resources |
|--------|---------|-----------|
| `main.sbatch` | Run the full pipeline as a batch job | 32 CPUs, 256GB, 72h |
| `resume_main.sbatch` | Resume a failed pipeline run | 32 CPUs, 256GB, 99h |
| `download_references.sbatch` | Download human reference genomes | 16 CPUs, 64GB, 400min |
| `download_microbial_db.sbatch` | Download KrakenUniq + MetaPhlAn databases | 16 CPUs, 64GB, 800min |
| `unpack_microbial_dbs.sbatch` | Untar KrakenUniq and MetaPhlAn database archives | 16 CPUs, 64GB, 800min |
| `create_minimap2_indexes.sbatch` | Build minimap2 indexes via SLURM | 16 CPUs, 64GB, 400min |
| `build_metaphlan_index.sbatch` | Build MetaPhlAn Bowtie2 index (32 threads) | 32 CPUs, 200GB, 1200min |
| `run_all_birdman.sh` | Run 4 BIRDMAn differential abundance models via QIIME2 | 32 CPUs, 128GB, 48h |

All batch scripts use `platinum` partition, `hcp-ddp302` QOS, `ddp302` account, and the `microbiome` conda environment (except `build_metaphlan_index.sbatch` which uses `metaphlan4_env`). Logs go to `scripts/batch_scripts/log/`.
