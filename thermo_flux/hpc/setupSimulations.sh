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
MEMORY=$7

# Setup an empty 'todolist' to store all the jobs to be submitted:
rm -f todolist.sh
touch todolist.sh

SLURM_MINS=$((MINS+1)) #add 1 min to job to allow gurobi to get to time limit and save .sol file 


# Setup batch scripts for each optimization:
for i in $(seq 1 $SEED)
do
	# Change the settings for the job submission:
	cat ~/FOLDER_WITH_SCRIPTS/batch_jobSubmission_template.sh | sed "s/JOBNAME/$NAME/g" | sed "s/CPUS/$CPU/g"  | sed "s/HOURS/$HOURS/g" | sed "s/MINS/$SLURM_MINS/g" | sed "s/SECS/$SECS/g" | sed "s/MEMORY/$MEMORY/g" > batch_script_seed$i.slurm

	# Append line to call Gurobi:
	echo "gurobi.sh ~/FOLDER_WITH_SCRIPTS/gurobi_optimization.py $NAME $HOURS $MINS $SECS $CPU $i 1 $MEMORY" >> batch_script_seed$i.slurm

	# Append job to the 'todolist':
	echo "sbatch  batch_script_seed$i.slurm" >> todolist.sh
done

# Turn the 'todolist' into an executable:
chmod +x todolist.sh