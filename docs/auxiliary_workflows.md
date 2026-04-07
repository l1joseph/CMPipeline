# Auxiliary Workflows

These are standalone workflows that run separately from the main pipeline.

## AntiSMASH -- Biosynthetic Gene Cluster Detection

Standalone workflow (`run_antismash.nf`) for identifying biosynthetic gene clusters (BGCs) in metagenome assemblies or MAGs. Run separately after nf-core/mag assembly.

```bash
# Run on MEGAHIT assembly contigs (recommended for cancer microbiome)
nextflow run run_antismash.nf --mode contigs -profile conda -resume

# Run on binned MAGs from nf-core/mag
nextflow run run_antismash.nf --mode mags -profile conda -resume
```

Two modes: `contigs` (uses prodigal-m for metagenomic gene finding) and `mags` (uses prodigal for single-genome). Module: `Modules/antismash.nf`. Config: `conf/antismash.config`. Output: `RESULTS/ANTISMASH/`. Uses `errorStrategy 'ignore'` for empty/tiny assemblies common in cancer microbiome samples.

## nf-core/mag -- Metagenome Assembly & Binning

Config (`conf/mag.config`) for running nf-core/mag on TSCC with host-depleted reads from CMPipeline's `MAPPED_READS` output. Uses Singularity (cache: `/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/singularity_cache`).

```bash
nextflow run nf-core/mag -profile singularity -c conf/mag.config \
    --input mag_samplesheet.csv --outdir RESULTS/MAG -resume
```

Resource allocations: MEGAHIT (16 CPUs, 128GB, 48h), GTDB-Tk (16 CPUs, 256GB, 48h), binning tools (8 CPUs, 64GB, 24h). Default `errorStrategy 'ignore'` for non-resource failures -- cancer microbiome samples often produce sub-500bp assemblies that downstream tools can't handle.
