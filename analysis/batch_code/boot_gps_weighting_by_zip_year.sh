#!/bin/bash
#SBATCH --partition=fasse
#SBATCH --account=dominici_lab
#SBATCH  -c 1
#SBATCH --job-name mqin_gps_weighting
#SBATCH --output boot_gps_weighting_by_zip_year.out
#SBATCH --error boot_gps_weighting_by_zip_year.err
#SBATCH --mem=16GB
#SBATCH --time=00:20:00
#SBATCH --array=1-2
#SBATCH --mail-user=michelleqin@college.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0-fasrc01 R/4.0.5-fasrc02 rstudio/1.1.453-fasrc01
unset R_LIBS_SITE

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_LIBS_USER=$HOME/apps/Causal_ADRD_batch/R_4.0.5:$R_LIBS_USER

# R CMD BATCH --no-save batch_test.R Logs/run${SLURM_ARRAY_TASK_ID}.Rout
# R CMD BATCH --no-save batch_test.R
R CMD BATCH --no-save batch_test.R batch_logs/run0.txt

