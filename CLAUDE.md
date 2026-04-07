# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CMPipeline is a Nextflow DSL2 metagenomics pipeline for cancer microbiome analysis. Processes unmapped BAM files through host depletion, taxonomic classification, consensus taxa, decontamination, and batch correction. Runs on TSCC via SLURM with conda environments.

## Running the Pipeline

```bash
# Interactive node
srun -N 1 -n 1 -c 8 --mem 125G -t 24:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash
conda activate microbiome
nextflow run main.nf

# Batch job
sbatch scripts/batch_scripts/main.sbatch          # Fresh run
sbatch scripts/batch_scripts/resume_main.sbatch    # Resume

# Resume after failure
nextflow run main.nf -resume
```

## Directory Layout

- `main.nf` -- Main workflow and all pipeline parameters
- `run_antismash.nf` -- Standalone antiSMASH workflow
- `nextflow.config` -- Executor config (SLURM), profiles
- `conf/` -- Resource configs (`base.config`, `mag.config`, `antismash.config`)
- `Modules/` -- Nextflow process definitions (one per `.nf` file)
- `scripts/` -- Python/R scripts called by Nextflow processes
- `scripts/batch_scripts/` -- SLURM sbatch scripts (logs in `log/`)
- `conda_envs/` -- Conda environment YAMLs (one per tool)
- `notebooks/` -- BIRDMAn differential abundance analysis (QIIME2)
- `RESULTS/` -- Pipeline output
- `docs/` -- Detailed documentation
- `logs/` -- Session and analysis logs
- `docs/EAC_Microbiome_Presentation.pptx` -- Figures and presentation slides

## Key Documentation

| Document | Contents |
|----------|----------|
| `docs/architecture.md` | 5-step pipeline architecture, module pattern, output structure |
| `docs/configuration.md` | All parameters, database paths, resources, conda envs |
| `docs/scripts_reference.md` | Python, R, shell, and SLURM batch script reference |
| `docs/birdman_analysis.md` | BIRDMAn models, metadata, formulas, SLURM scripts |
| `docs/auxiliary_workflows.md` | AntiSMASH and nf-core/mag standalone workflows |
| `tickets.md` | Open work items and tracking |

## Essential Notes

- Nextflow DSL2 -- all modules use `nextflow.enable.dsl=2`
- Conda environments are per-process; do not consolidate
- All params defined in `main.nf` header, not `nextflow.config`
- Module pattern: skip-check -> execution -> conda env -> scratch. Preserve skip-check when modifying modules.
- `samples.csv` format: `patient,bam` (absolute BAM paths)
- `work/` is Nextflow intermediate storage, can grow very large
- Production decontam script: `scripts/251006_decontamination_ver2.R` (not `decontamination.R`)

## Git Conventions

- No `Co-Authored-By` lines in commits or PRs
- No auto-generated commit trailers
- Keep commit messages concise and descriptive

## Claude Code Preferences

- Keep conversations short and low context
- Reference `docs/` files rather than re-explaining architecture
- Use `EAC_Microbiome_Presentation.pptx` as the central location for figures and slides
- Log notable findings or session outputs to `logs/`
