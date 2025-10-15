#!/bin/bash
#SBATCH --job-name=thermo2
#SBATCH --cpus-per-task=4
#SBATCH --time=3-00:00:30
#SBATCH --output=job-%j.log
#SBATCH --partition=regular

echo "SLURM_JOBID: " $SLURM_JOBID

module purge
module load Gurobi/10.0.0-GCCcore-11.3.0

gurobi.sh gurobi_solve_cov.py $seed $model 



