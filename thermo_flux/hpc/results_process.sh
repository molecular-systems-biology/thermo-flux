#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16000
 
module purge
module load Python/3.9.5-GCCcore-10.3.0
module load Gurobi/9.5.0-GCCcore-10.3.0 
 
source /data/$USER/.envs/thermo_py/bin/activate
 
python3 --version
which python3
 
python3 -u analyse_results.py > test2.out

deactivate