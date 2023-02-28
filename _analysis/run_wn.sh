#!/bin/bash

#SBATCH --partition=fasse
#SBATCH -c 36 # Number of cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH -J cgps_adrd      # Job name 
#SBATCH --mem 136000       # Memory request (92 Gb)
#SBATCH -t 05-00:00       # Maximum execution time (D-HH:MM)
#SBATCH -o cgps_adrd_%a_%j.out  # Standard output
#SBATCH -e cgps_adrd_%a_%j.err  # Standard error
#SBATCH --array=1-40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nkhoshnevis@g.harvard.edu

module load R/4.2.0-fasrc01
module load GCCcore/7.3.0 LLVM/6.0.0


unset R_LIBS_SITE
export R_LIBS_USER=$HOME/apps/R_420:$R_LIBS_USER

export JOB_ARRAY_ID=`echo $SLURM_ARRAY_TASK_ID`

srun Rscript gps_models_nk_wn.R $SLURM_ARRAY_TASK_ID






