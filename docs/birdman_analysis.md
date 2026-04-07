# BIRDMAn Differential Abundance Analysis

All BIRDMAn models run via QIIME2 (`conda activate q2-birdman`). Results live under `notebooks/out/`. Each model directory contains `table.qza`, `metadata.tsv`, `results.qza`, `plot.qzv`, and `results_export/metadata.tsv` (the parsed output).

## Metadata

- **Primary metadata**: `RESULTS/04_DECONTAMINATION/bracken_run/EAC_GEJ_METADATA.txt` (75 columns, keyed by `donor_id`, 1040 samples: 520 Tumor + 520 Control)
- **Extended metadata**: `notebooks/metadata/MUT_EAC_GEJ_extended_metadata.csv` (includes mutation data)
- **SBS signatures**: `notebooks/metadata/pruned_attribution_EAC_SBS96_COSMIC_H_test_SBS_abs_mutations.csv` (SBS96 COSMIC signature attributions)

## Model Set 1: Tumor vs Normal (Type effect)

Tests `Type` (Tumor vs Normal) on all 1040 samples. Uses batch-corrected ConQuR counts. The **preferred models** are in `birdman_genus/` and `birdman_species/` (genus/species-level corrected data, unrarefied, with `log_total_reads` offset).

| Directory | Input Data | Models | Status |
|-----------|-----------|--------|--------|
| `notebooks/out/birdman_genus/` | Genus-level ConQuR-corrected (257 taxa, 1040 samples) | type, type_country, type_demographic, type_lifestyle, type_full | **Primary -- use these** |
| `notebooks/out/birdman_species/` | Species-level ConQuR-corrected (351 taxa, 1040 samples) | type, type_country, type_demographic, type_lifestyle, type_full | **Primary -- use these** |
| `notebooks/out/birdman/` | Pre-correction rarefied + unrarefied | type_batch, type_country, type_lifestyle, full + unrarefied variants | Legacy, for comparison only |
| `notebooks/out/birdman_corrected/` | ConQuR-corrected rarefied + unrarefied | type, type_country (rarefied + unrarefied) | Legacy, for comparison only |

**Type model formulas** (from `scripts/batch_scripts/run_birdman_genus.sh` / `run_birdman_species.sh`):
- `type`: `Type`
- `type_country`: `Type + country`
- `type_demographic`: `Type + age_diag + sex`
- `type_lifestyle`: `Type + tobacco + alcohol`
- `type_full`: `Type + country + age_diag + sex`

**Data prep script**: `notebooks/prepare_birdman_filtered.py` -- filters corrected TSV to genus/species, creates BIOM/QZA, sets up model directories with symlinks.

**Visualization**: `notebooks/10_birdman_visualization.ipynb` -- forest plots, Jaccard overlap heatmaps, cross-model LFC concordance, core taxa analysis. Outputs in `notebooks/out/figs_birdman/`.

## Model Set 2: Within-Tumor Clinical Associations (tumor-only)

Tests which microbes associate with clinical features **within tumor samples only** (520 samples), with all other clinical covariates as confounders. Uses unrarefied ConQuR-corrected counts with `log(Total_Reads)` offset (rarefaction not needed -- the offset handles library size).

**Framework**:
1. Prevalence filter taxa to >=10% of tumor samples
2. Handle missing data: categorical -> "Missing" category; `stage` consolidated to 5 levels (0_I, II, III, IV, Missing); `bmi` imputed with median
3. Run VIF check on design matrix before modelling
4. For each feature, confounders = all other base covariates. Three model variants handle batch-country confounding:
   - **batch**: all 520 tumor samples, `shipment_batch` as confounder
   - **country**: all 520 tumor samples, `country` as confounder
   - **multicountry**: 223 samples (exclude UK-only batch), both `country` + `shipment_batch`

**Features tested**: `tumorsite`, `stage`, `GERD_tx`, `cigstatus`, `alcstatus`, `bmi`

**Base formula** (example for tumorsite, batch variant):
```
tumorsite + age_diag + sex + bmi + cigstatus + alcstatus + stage + GERD_tx + shipment_batch + log_total_reads
```

| Directory | Input Data | Models | Status |
|-----------|-----------|--------|--------|
| `notebooks/out/birdman_association_genus/` | Genus ConQuR-corrected, tumor-only (257 taxa / 267 multicountry) | 18 models (6 features x 3 variants) | **Complete** |
| `notebooks/out/birdman_association_species/` | Species ConQuR-corrected, tumor-only (322 taxa / 308 multicountry) | 18 models (6 features x 3 variants) | **Complete** |

Each directory also contains: `table_all.qza` (520 samples), `table_multicountry.qza` (223 samples), `metadata_all.tsv`, `metadata_multicountry.tsv`, `model_manifest.tsv` (all formulas), `vif_report.tsv`, `vif_summary.txt`.

**Scripts**:
- Data prep: `notebooks/prepare_birdman_association.py`
- VIF check: `notebooks/run_vif_check.py`
- SLURM batch: `scripts/batch_scripts/run_birdman_association_genus.sh` / `run_birdman_association_species.sh` (32 CPUs, 128GB, 96h)
- Visualization: `notebooks/15_birdman_association.ipynb` -- forest plots per feature/variant/effect level, cross-variant concordance, cross-feature heatmaps, VIF display, summary CSV. Outputs in `notebooks/out/figs_association/`.

## BIRDMAn SLURM Batch Scripts Summary

| Script | Purpose | Resources |
|--------|---------|-----------|
| `run_birdman_genus.sh` | Type models on genus-level corrected data (5 models) | 32 CPUs, 128GB, 48h |
| `run_birdman_species.sh` | Type models on species-level corrected data (5 models) | 32 CPUs, 128GB, 48h |
| `run_birdman_association_genus.sh` | Clinical association models, genus (18 models) | 32 CPUs, 128GB, 96h |
| `run_birdman_association_species.sh` | Clinical association models, species (18 models) | 32 CPUs, 128GB, 96h |
| `run_all_birdman.sh` | Legacy: Type models on pre-correction data | 32 CPUs, 128GB, 48h |
| `run_birdman_corrected.sh` | Legacy: Type models on rarefied corrected data | 32 CPUs, 128GB, 48h |
| `run_birdman_unrarefied.sh` | Legacy: Type models on unrarefied pre-correction data | 32 CPUs, 128GB, 48h |
| `run_birdman_corrected_unrarefied.sh` | Legacy: Type models on unrarefied corrected data | 32 CPUs, 128GB, 48h |

All BIRDMAn batch scripts use `q2-birdman` conda env, `platinum` partition, `hcp-ddp302` QOS.
