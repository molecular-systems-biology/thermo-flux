#!/bin/bash
#
#Script to setup simulations on the Habrok cluster
#
#José Losa
#2023.05.26


# Store the input arguments in variables with understandable names:
NAME=$1
SEED=$2
CPU=$3
HOURS=$4
MINS=$5
SECS=$6
WRITE_SOLS=$7
MEMORY=$8


# Setup an empty 'todolist' to store all the jobs to be submitted:
rm -f todolist.sh
touch todolist.sh

SLURM_MINS=$((MINS+2)) #add 2 min to job to allow gurobi to get to time limit and save .sol file 

# Setup batch scripts for each optimization:
for file in *.mps; do
  
	# Change the settings for the job submission:
	cat ~/FOLDER_WITH_SCRIPTS/batch_jobSubmission_template.sh | sed "s/JOBNAME/$NAME/g" | sed "s/CPUS/$CPU/g"  | sed "s/HOURS/$HOURS/g" | sed "s/MINS/$SLURM_MINS/g" | sed "s/SECS/$SECS/g" | sed "s/MEMORY/$MEMORY/g" > batch_script_${file%.*}.slurm

	# Append line to call Gurobi:
	echo "gurobi.sh ~/FOLDER_WITH_SCRIPTS/gurobi_optimization.py ${file%.*} $HOURS $MINS $SECS $CPU $SEED $WRITE_SOLS" >> batch_script_${file%.*}.slurm

	# Append job to the 'todolist':
	echo "sbatch  batch_script_${file%.*}.slurm" >> todolist.sh
	
done

# Turn the 'todolist' into an executable:
chmod +x todolist.sh