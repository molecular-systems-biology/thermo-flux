### How to set up the jobs for submission to Habrók

1. Create a folder on your main Habrok directory:
	```
	cd ~
	mkdir FOLDER_WITH_SCRIPTS
	```

2. Add the following scripts to the folder created above (just drag and drop from the local computer): `setupSimulations.sh`, `batch_jobSubmission_template.sh`, and `gurobi_optimization.py`

3. Move to the `/scratch` directory and create a folder to store the files of a given optimization
	```
	cd /scratch/$USER
	mkdir FOLDER_WITH_OPTIMIZATION_RESULTS
	```

4. Add the `.mps` file describing the optimization problem (obtained from running the Jupyter Notebook) to the folder above

5. Move into the folder created:
	```
	cd /scratch/FOLDER_WITH_OPTIMIZATION_RESULTS
	```

5. Define the list of jobs to be run by calling the `setupSimulations.sh` script with 4 arguments, as follows:
	```
	bash ~/FOLDER_WITH_SCRIPTS/setupSimulations.sh NAME_OF_MPS_FILE NUMBER_OF_SEEDS NUMBER_OF_CPUS HOURS
	```

6. Submit the list of jobs to be run:
	```
	./todolist
	```

### Note:

If you want to define the starting point for some optimization variables, save an `.mts` file with the same base name as the model (i.e., the `.mps` file). Both files should be uploaded to the cluster.


Upon running `setupSimulations.sh`, the `.mts` file will be automatically recognized and used to set the MIP start.