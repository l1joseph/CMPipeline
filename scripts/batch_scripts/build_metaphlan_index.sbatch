#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 1200
#SBATCH -p platinum
#SBATCH -q hcp-ddp302
#SBATCH -A ddp302
#SBATCH --mail-user=l1joseph@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=CMP_metaphlan_index
#SBATCH --output=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/scripts/batch_scripts/log/CMP_metaphlan_index_%j.out
#SBATCH --error=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/scripts/batch_scripts/log/CMP_metaphlan_index_%j.err

# Load Conda
source ~/.bashrc

conda activate metaphlan4_env
echo "conda env activated"
echo "Bowtie2 version:"
bowtie2 --version | head -1

cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/out/ref/metaphlan
echo "cd metaphlan dir"

echo "Decompressing .bz2 files..."
# bzip2 -dk mpa_vJun23_CHOCOPhlAnSGB_202307_SGB.fna.bz2
# bzip2 -dk mpa_vJun23_CHOCOPhlAnSGB_202307_VSG.fna.bz2

echo "Building Bowtie2 index with 32 threads..."
bowtie2-build --large-index --threads 32 \
  mpa_vJun23_CHOCOPhlAnSGB_202307_VSG.fna \
  mpa_vJun23_CHOCOPhlAnSGB_202307

echo "DONE - Bowtie2 index built successfully"
echo "Final files:"
ls -lh
