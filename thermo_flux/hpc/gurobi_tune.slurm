#!/bin/bash
#SBATCH --job-name=tune
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:05:00
#SBATCH --output=job-%j.log
#SBATCH --partition=regular

echo "SLURM_JOBID: " $SLURM_JOBID

module purge
module load Python
module load Gurobi/10.0.1-GCCcore-11.3.0

gurobi.sh gurobi_tune.py



