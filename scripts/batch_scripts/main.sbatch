#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=256G
#SBATCH -t 72:00:00
#SBATCH -p platinum
#SBATCH -q hcp-ddp302
#SBATCH -A ddp302
#SBATCH --mail-user=l1joseph@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=CMP_EAC_GEJ
#SBATCH --output=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy/scripts/batch_scripts/log/CMP_nextflow_%j.out
#SBATCH --error=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy/scripts/batch_scripts/log/CMP_nextflow_%j.err

# Load Conda
source ~/.bashrc

# Activate nextflow environment
conda activate microbiome
echo "Nextflow microbiome environment activated"

# Change to pipeline directory
cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy/
echo "Working directory: $(pwd)"

# Run nextflow - parameters already set in main.nf
echo "Starting Nextflow pipeline..."
nextflow run main.nf \
    -resume \
    -N l1joseph@ucsd.edu

echo "Pipeline completed"