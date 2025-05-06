# CMPipeline

## Workflow Introduction
<img src="https://github.com/ammalabbasi/CMpipeline/blob/main/workflow_logo/v0.1.png" width="95%" height="95%">

## How to run CMPipeline
1.  Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) (e.g., using Conda).
2.  Make sure you have the following files and directory structure in your working directory:
    * `main.nf`
    * `nextflow.config` (Note: Based on the provided files, the primary config might be `conf/base.config` rather than `nextflow.config` at the top level)
    * `conf/base.config`
    * `Modules/` (containing all `.nf` script files, such as `extract_reads.nf`, `map_both_reads.nf`, etc.)
    * `sample.csv`

    > You might find template or shared versions of these files in locations like `/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/CMPipeline` (adjust path as needed).
3.  Prepare your `sample.csv` file:
    This file tells the pipeline which samples to process and where their input BAM files are located. It must contain `patient` (sample ID) and `bam` (full path to BAM file) columns. This corresponds to the essential `params.sample` parameter used by the pipeline.
    Example `sample.csv`:
    ```csv
    patient,bam
    PD56137a,/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/test_data/PD56137a.unmapped.viral.bam
    PD56137b,/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/test_data/PD56137b.unmapped.viral.bam
    ```

4.  Request an interactive node (if using a cluster like TSCC) and run Nextflow:

    ```bash
    # Example Node Request (adjust resources and account details as needed)
    srun -N 1 -n 1 -c 8 --mem 125G -t 24:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash

    # Activate your nextflow conda environment (replace env_nf if needed)
    conda activate env_nf

    # Run nextflow from your working directory
    # Ensure main.nf is the correct entry point script name based on your files
    nextflow run main.nf

    # Resume after interruption:
    # If the pipeline terminates unexpectedly, you can resume from the last successful step
    nextflow run main.nf -resume

    # Optional: Receive an email notification on completion (-N flag)
    # nextflow run main.nf -N your_email@example.com
    ```
5.  Results will be stored in the **RESULTS** folder (or specific subdirectories like `UNMAPPED_BAM`, `MAPPED_READS`, `BRACKEN`, etc., as defined by parameters in `main.txt`).

## Pipeline Details

### Overview

This Nextflow pipeline processes sequencing data (starting from BAM files) to identify microbial composition after removing host (human) reads. The main workflow steps are:

1.  **Extract Unmapped Reads:** Extracts paired-end unmapped reads (R1 and R2) from input BAM files.
2.  **Quality Control (QC):** Performs initial QC on raw reads (FastQC).
3.  **Read Filtering:** Filters reads for quality, length, and adapters (fastp).
4.  **Post-filtering QC:** Performs QC on filtered reads (FastQC).
5.  **Host Read Removal:** Sequentially maps filtered reads against human references (hg38, T2T+PhiX) and a pangenome database (minimap2), keeping only unmapped reads. QC (FastQC) is performed after each mapping stage.
6.  **Taxonomic Classification:** Performs classification and abundance estimation on non-host reads (KrakenUniq/Bracken). (Note: MetaPhlAn4 module is included but commented out).
7.  **Reporting (Optional):** Includes modules for result aggregation (MultiQC, consensus taxa merging), potentially requiring activation.

### Paired-End Read Handling

A key feature of this pipeline's workflow (`main.nf`) is how it handles paired-end reads. Instead of managing separate streams for R1 and R2 files, it groups them early on:

1.  **Tuple Output:** The `extractReads` process outputs a channel where each item is a tuple containing the sample ID, the R1 path, and the R2 path: `tuple(sampleID, r1_path, r2_path)`.
2.  **`multiMap` for Flexibility:** This output channel is immediately processed using the `multiMap` operator. This creates different "views" of the channel for downstream processes: one view containing just the paired paths `tuple(r1_path, r2_path)` (e.g., for `FASTQC`), and another containing the ID plus the paths `tuple(sampleID, r1_path, r2_path)` (e.g., for `filterReads`, `mapBothReads`).
3.  **Consistent Pairing:** This approach ensures that R1 and R2 reads for a sample stay linked together throughout the pipeline, simplifying process inputs and reducing the risk of pair mismatches.
4.  **Contrast with Previous Method:** Commented-out code within `main.nf` suggests an older approach where R1 and R2 reads were likely passed as separate inputs to processes like `FASTQC` and required separate mapping steps (`mapReadsR1`, `mapReadsR2`). The current tuple-based method provides a more robust and streamlined way to manage paired-end data in Nextflow.

### Module Descriptions

The pipeline utilizes the following modules located in the `Modules/` directory:

* **`extract_reads.nf` (`extractReads` process):**
    * **Purpose:** Extracts unmapped paired-end reads from a coordinate-sorted BAM file.
    * **Tools:** `samtools view`, `samtools bam2fq`.
    * **Input:** `val(meta)` (containing `meta.patient`, `meta.bam`).
    * **Output:** `tuple val(<sampleID>), path(<sampleID>.R1.UNMAPPED.fastq.gz), path(<sampleID>.R2.UNMAPPED.fastq.gz)`.

* **`fastqc.nf` (`FASTQC` process):**
    * **Purpose:** Assesses FASTQ quality. Used multiple times (aliases: `FASTQC1`, `FASTQC2`, `FASTQCHG38`, etc.).
    * **Tools:** `fastqc`.
    * **Input:** `tuple path(R1.fastq.gz), path(R2.fastq.gz)`.
    * **Output:** `path(*_fastqc.html)`, `path(*_fastqc.zip)`.

* **`filter_reads.nf` (`filterReads` process):**
    * **Purpose:** Filters paired-end reads (quality, length, adapters).
    * **Tools:** `fastp`.
    * **Input:** `tuple val(<sampleID>), path(R1_raw.fastq.gz), path(R2_raw.fastq.gz)`.
    * **Output:** `tuple val(<sampleID>), path(<sampleID>.R1.UNMAPPED.FASTP.FILTERED.fastq.gz), path(<sampleID>.R2.UNMAPPED.FASTP.FILTERED.fastq.gz)`.

* **`map_both_reads.nf` (`mapBothReads` process):**
    * **Purpose:** Sequentially removes host reads via mapping (hg38, T2T, Pangenome).
    * **Tools:** `minimap2`, `samtools fastq`.
    * **Input:** `tuple val(<sampleID>), path(R1_filt.fastq.gz), path(R2_filt.fastq.gz)`, `val(mmi_files)`.
    * **Output:** `tuple val(<sampleID>), path(<R1_hg38_unmapped>), path(<R1_t2t_unmapped>), path(<R1_pan_unmapped>), path(<R2_hg38_unmapped>), path(<R2_t2t_unmapped>), path(<R2_pan_unmapped>)` (representing the 6 FASTQ files).
    * **Note:** Replaces the older `map_reads.nf` by handling R1 and R2 together in one process, simplifying the workflow.

* **`metaphlan4.nf` (`metaphlan4` process):**
    * **Purpose:** Taxonomic profiling using MetaPhlAn 4.
    * **Tools:** `cat`, `gunzip`, `metaphlan`.
    * **Input:** `tuple path(R1_final.fastq.gz), path(R2_final.fastq.gz)`.
    * **Output:** `path(*.metagenome.bowtie2.bz2)`, `path(*.sam.bz2)`, `path(*.profiled_metagenome.txt)`.
    * **Execution Status:** **Currently skipped by default.** The call to this process is commented out in `main.txt`. This is often done for efficiency in low-microbiome environments, as MetaPhlAn typically requires a substantial number of microbial reads (e.g., >100,000) for reliable profiling.
    * **How to Enable:** To run this step, you need to uncomment the relevant line in the `main.txt` file. Find the line containing `// metaphlan4(...)` within the `workflow { ... }` block and remove the leading `//`.

* **`Bracken.nf` (`Bracken` process):**
    * **Purpose:** Taxonomic classification (KrakenUniq) and abundance estimation (Bracken).
    * **Tools:** `cat`, `gunzip`, `krakenuniq`, `bracken`.
    * **Input:** `tuple path(R1_final.fastq.gz), path(R2_final.fastq.gz)`.
    * **Output:** `path(*.krakenuniq.report.txt)`, `path(*.classified.fasta)`, `path(*.unclassified.fasta)`, `path(*.bracken.*.report.txt)`, `path(*.bracken.*.krakenreport.txt)`.

* **`multiqc.nf` (`multiqc` process):**
    * **Purpose:** Aggregates QC results into a single report.
    * **Tools:** `multiqc`.
    * **Input:** `val(input_files)` (channel collecting QC outputs).
    * **Output:** `path(multiqc_report.html)`, `path(multiqc_data.tsv)`.
    * **Execution Status:** **NOTE:** While this module is included with aliases (e.g., `MCR1`) in `main.nf`, there is **no active call** to execute it within the main `workflow` block. To generate a MultiQC report, the relevant process call (e.g., `MCR1(...)` or similar) needs to be added to the `workflow` block in `main.nf`, collecting the desired input channels (e.g., FastQC outputs).

* **`consensus.taxa.nf` (Contains `process_metaphlan4`, `process_bracken` processes):**
    * **Purpose:** Merges taxonomic results across samples.
    * **Tools:** `merge_metaphlan_tables.py`, `kreport2mpa.py`, `combine_mpa.py`, `grep`, `sed`.
    * **Input:** Paths to profile/report files from multiple samples.
    * **Output:** Merged tables (`merged_*.txt`) or combined reports (`crc.bracken.*.report.txt`).
    * **Execution Status:** **NOTE:** The processes defined within this module (`process_metaphlan4`, `process_bracken`) are **not currently called** in the `main.nf` workflow block. To merge results across samples, these process calls need to be added to the `workflow` block in `main.nf`, providing the appropriate collective input channels (e.g., all Bracken reports).

### Skip Logic for Existing Files

Most processes implement logic to skip execution if expected output files already exist in their designated `publishDir`. This check (`if [[ -f "output_file1" && -f "output_file2" ]]`) happens at the beginning of the `script:` block. If files are found, the script prints a skip message and exits successfully (`exit 0`), preventing redundant computation and allowing efficient pipeline resumption. Optionally, existing files may be linked or copied to the current working directory.