#! /bin/sh
#SBATCH -n 1               # (Max) number of tasks per job, for R usually 1
#SBATCH -p airoldi         # Partition to use
#SBATCH --mem=100 # Memory required per CPU, in MegaBytes
#SBATCH -a 1-10           # Array of 500 jobs, with ids 1, 2, ..., 500
#SBATCH -t 20

# And the code to run:
Rscript -e "source('andrewliu_mcmc.R') ; runMCMC(data=Impala) " ${SLURM_ARRAY_TASK_ID}
