#!/bin/bash
#SBATCH --job-name=ye_reg
#SBATCH --cpus-per-task=8
#SBATCH --time=16:30:00
#SBATCH --output=job-%j.log
#SBATCH --partition=regular

echo "SLURM_JOBID: " $SLURM_JOBID

module purge
module load Gurobi/10.0.0-GCCcore-11.3.0

gurobi.sh gurobi_solve_cov.py $seed 



