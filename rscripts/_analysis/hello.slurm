#!/bin/bash

#SBATCH -J hello
#SBATCH -p fasse
#SBATCH --array=1-10
#SBATCH -c 1
#SBATCH -t 10:00
#SBATCH -o hello-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nkhoshnevis@g.harvard.edu


module load R/4.2.0-fasrc01

srun Rscript hello.R $SLURM_ARRAY_TASK_ID