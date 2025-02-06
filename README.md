# CMPipeline

## Workflow Introduction
<img src="https://github.com/ammalabbasi/CMpipeline/blob/main/workflow_logo/v0.1.png" width="95%" height="95%">

## How to run EVC pipeline
1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) as a conda environment
2. Make sure you have the following folders and files in your working directory:
  - main.nf
  - nextflow.config
  - conf/base.config
  - sample.csv

  > You can find these shared files and folder in `/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/CMPipeline`

3. Next, download the human reference genomes to be used for filtration. We recommend [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/), [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/), and all currently available pangenomes from the [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org). A download script is provided for convenience.
```bash
bash scripts/download_references.sh
```

4. Next, create Minimap2 indexes for the previously downloaded reference genomes. A script is provided for convenience to build minimap2 indexes.
```bash
bash scripts/create_minimap2_indexes.sh
```

3. Prepare your sample.csv file:
```
patient,bam
PD56137a,/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/test_data/PD56137a.unmapped.viral.bam
PD56137b,/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/test_data/PD56137b.unmapped.viral.bam
```

4. Request an interactive node and run Nextflow in your working directory under an interactive node:

```
# Node requesting
srun -N 1 -n 1 -c 8 --mem 125G -t 24:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash

# Activate your nextflow conda environment
conda activate env_nf

# Run nextflow
nextflow run main.nf

# If your pipeline terminates with external error, or the interactive node is killed, you can resume your task after setting up the previous steps again with the following command:
nextflow run main.nf -resume

# Optionally, you can recieve an notifiction email on completion with -N flag:
nextflow run main.nf -N your_email@gmail.com
```
5. Every process result and report will be stored in the **RESULT** folder
