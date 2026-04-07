#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 800
#SBATCH -p platinum
#SBATCH -q hcp-ddp302
#SBATCH -A ddp302
#SBATCH --mail-user=l1joseph@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=CMP_unpack
#SBATCH --output=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/scripts/batch_scripts/log/CMP_unpack_%j.out
#SBATCH --error=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/scripts/batch_scripts/log/CMP_unpack_%j.err


# Load Conda
source ~/.bashrc

conda activate microbiome
echo "conda env activated"

cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/out/ref/
echo "cd out dir"

echo "unpacking kraken"
tar -xzf krakenUniq/kuniq_microbialdb_minus_kdb.20230808.tgz -C krakenUniq/

echo "unpacking metaphlan"
tar -xf metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202307.tar -C metaphlan/
echo "DONE"

