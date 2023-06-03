#!/bin/bash
#SBATCH --partition=fasse
#SBATCH --account=dominici_lab
#SBATCH  -c 1
#SBATCH --job-name=mqin_associational
#SBATCH --output=boot_associational.out
#SBATCH --error=boot_associational.err
#SBATCH --mem=8000
#SBATCH --time=00:05:00
#SBATCH --array=1-500
#SBATCH --mail-user=michelleqin@college.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0-fasrc01 R/4.0.5-fasrc02 rstudio/1.1.453-fasrc01
unset R_LIBS_SITE
export R_LIBS_USER=$HOME/apps/Causal_ADRD/R_4.0.5:$R_LIBS_USER

export JOB_ARRAY_ID=$SLURM_ARRAY_TASK_ID
Rscript boot_associational.R $SLURM_ARRAY_TASK_ID
