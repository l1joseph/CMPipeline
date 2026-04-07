#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 400
#SBATCH -p platinum
#SBATCH -q hcp-ddp302
#SBATCH -A ddp302
#SBATCH --mail-user=l1joseph@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=CMP_download_genomes
#SBATCH --output=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/scripts/batch_scripts/log/CMP_download_genomes_%j.out
#SBATCH --error=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/scripts/batch_scripts/log/CMP_download_genomes_%j.err


# Load Conda
source ~/.bashrc

conda activate microbiome
echo "conda env activated"

cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/out/
echo "cd ref dir"

bash ../scripts/download_references.sh

echo "DONE"
