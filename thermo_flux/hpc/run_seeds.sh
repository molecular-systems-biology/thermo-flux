#!/bin/bash
for ((seed=0; seed<=25; seed++))
  do
      sbatch --export=seed=$seed, gurobi_solve_cov.slurm
      echo $model $seed 
  done
  

     
 




