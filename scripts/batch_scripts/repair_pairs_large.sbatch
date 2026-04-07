#!/bin/bash
#SBATCH -p platinum
#SBATCH -q hcp-ddp302
#SBATCH -A ddp302
#SBATCH --job-name=repair_large
#SBATCH --output=scripts/batch_scripts/log/repair_large_%A_%a.out
#SBATCH --error=scripts/batch_scripts/log/repair_large_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --array=1-3
#SBATCH --mail-user=l1joseph@ucsd.edu
#SBATCH --mail-type=END,FAIL

# Re-pair the 3 high-read-count samples that OOM'd with 16GB

set -e
set -o pipefail

PROJ_DIR="/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy"
OUT_DIR="${PROJ_DIR}/RESULTS/MAPPED_READS_PAIRED"

source ~/.bashrc
conda activate microbiome

cd "${PROJ_DIR}"

# The 3 samples that failed
SAMPLES=(PD44801a PD44819a PD60613a)
SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"

R1="RESULTS/MAPPED_READS/${SAMPLE}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz"
R2="RESULTS/MAPPED_READS/${SAMPLE}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz"

echo "=========================================="
echo "Job: ${SLURM_JOB_NAME} (task ${SLURM_ARRAY_TASK_ID})"
echo "Sample: ${SAMPLE} (high read count - 128GB memory)"
echo "Started: $(date)"
echo "=========================================="

python scripts/repair_pairs.py \
    --r1 "${R1}" \
    --r2 "${R2}" \
    --out-dir "${OUT_DIR}"

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
