#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --cpus-per-task=CPUS
#SBATCH --time=HOURS:30:00
#SBATCH --output=job-%j.log
#SBATCH --partition=regular

echo "SLURM_JOBID: " $SLURM_JOBID

module purge
module load Gurobi/10.0.1-GCCcore-11.3.0

# Call Gurobi:
