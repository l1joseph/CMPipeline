# Configuration Reference

## Key Parameters

All pipeline parameters are defined at the top of `main.nf` (not in `nextflow.config`).

### Workflow Control

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `params.sample` | Path to samples.csv (columns: `patient`, `bam`) | — |
| `params.metadata_file` | Metadata TSV for decontam/batch correction | EAC_GEJ_METADATA.txt |
| `params.start_from` | Entry point: `beginning`, `decontam`, `batch_correction` | `beginning` |
| `params.skip_host_depletion` | Skip host depletion step | `false` |
| `params.skip_classification` | Skip taxonomic classification | `false` |
| `params.skip_consensus` | When true, Bracken output goes directly to decontam | `false` |
| `params.skip_humann3` | Skip HUMAnN3 functional profiling | `true` |
| `params.run_decontam` | Enable decontamination step | `true` |
| `params.run_batch_correction` | Enable batch correction step | `false` |
| `params.consensus_otu_table` | Pre-computed OTU table for `start_from=decontam` | `null` |
| `params.decontam_otu_table` | Pre-computed table for `start_from=batch_correction` | `null` |

### Decontamination Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `decontam_threshold` | 0.1 | Prevalence threshold for contaminant detection |
| `decontam_min_prevalence` | 0.02 | Minimum cohort-wide prevalence for rare taxa filter |
| `decontam_min_abundance` | 5 | Minimum sample-wise abundance for rare taxa filter |
| `decontam_min_batches` | 2 | Minimum batches a taxon must appear in |
| `batch_column` | `shipment_batch` | Metadata column for batch ID |
| `type_column` | `Type` | Metadata column for sample type |
| `tumor_values` | `Tumor` | Values to classify as tumor |
| `control_values` | `Normal` | Values to classify as control |

### Batch Correction Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `batch_corr_covariates` | `age_diag,sex,bmi` | Comma-separated covariate list |
| `tumor_value` | `Tumor` | Tumor label in metadata |
| `tumor_only` | `false` | If true, only correct tumor samples |
| `phase` | `auto` | ConQuR phase: `auto`, `1`, or `3` |
| `r2_threshold` | 0.25 | PERMANOVA R2 reduction threshold for auto-phase stopping |

## Database Paths

| Database | Path |
|----------|------|
| hg38 | `/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi` |
| T2T + PhiX | `/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi` |
| Pangenome | `/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/pangenome_mmi` |
| KrakenUniq | `/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023` |
| MetaPhlAn | `/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/metaphlan` |
| HUMAnN3 nucleotide | `/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/databases/humann3/chocophlan/` |
| HUMAnN3 protein | `/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/databases/humann3/uniref/` |

## Resource Configuration

`conf/base.config` defines SLURM resources per process label. All jobs submit to `platinum` partition with `hcp-ddp302` QOS via `clusterOptions`. Resources scale with `task.attempt` on retry (up to 5 retries). Key labels:

- `process_high`: 16 CPUs, 256GB, 50h (Bracken, MetaPhlAn4, HUMAnN3, BatchCorrection)
- `process_medium`: 4 CPUs, 128GB, 30h (Decontamination, MultiQC)
- `process_high_disk`: 4 CPUs, 256GB, 50h, 100GB disk (Bracken)
- `mapReads`: 4 CPUs, 64GB, 50h (scales with retry)
- `extractReads`/`filterReads`/`fastQC`: 4 CPUs, 64-128GB, 8-10h (fixed, no retry scaling)

## Conda Environments

Each Nextflow process uses a separate conda env defined in `conda_envs/`. Do not consolidate them.

| Environment File | Env Name | Key Packages |
|-----------------|----------|--------------|
| `samtools_env.yml` | -- | samtools (also used for FastQC) |
| `fastp_env.yml` | -- | fastp |
| `minimap2_env.yml` | -- | minimap2, samtools |
| `krakenUniq_bracken_env.yml` | -- | krakenuniq, bracken |
| `metaphlan4_env.yml` | -- | metaphlan (v4), bowtie2 |
| `humann3_env.yml` | -- | humann (v3) |
| `consensus_taxa_env.yml` | -- | pandas, matplotlib, seaborn, adjustText |
| `decontam_env.yml` | `nf-decontam-env` | R 4.3, phyloseq, decontam, vegan, tidyverse, ggrepel |
| `batch_correction_env.yml` | `nf-batch-correction-env` | R 4.3, vegan 2.6.4 (pinned), ConQuR (GitHub install), doParallel, zCompositions, DEICODE |
| `multiqc_env.yml` | -- | multiqc |
| `fastqc_env.yml` | -- | fastqc |
| `antismash_env.yml` | `nf-antismash-env` | antiSMASH (BGC detection) |

**Note on batch correction env**: ConQuR must be installed separately after creating the conda env: `R -e 'devtools::install_github("ivartb/ConQuR_par")'`. The `ivartb/ConQuR_par` fork fixes a `batchid not found` error in foreach workers. Vegan is pinned to 2.6.4 to avoid `adonis` deprecation issues.

## Default Metadata

The default metadata file points to `/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/laura/data/EAC_GEJ_METADATA.txt` -- update `params.metadata_file` in `main.nf` for different datasets.
