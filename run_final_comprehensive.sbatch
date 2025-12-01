#!/bin/bash
#SBATCH --job-name=raf_repair
#SBATCH --output=logs/raf_repair_%a.out
#SBATCH --error=logs/raf_repair_%a.err
#SBATCH --array=1-1680
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

module load r/4.4.0

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /home/anadiger_umass_edu/stats535_rafeast

# No need to delete anything. The script checks for files.
Rscript sweep_final_comprehensive.R ${SLURM_ARRAY_TASK_ID}