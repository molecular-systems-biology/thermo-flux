.. _installation:

Installation
============
	
1. To avoid dependency conflicts it is reccomended to use a python environment e.g:

	```conda create -n thermoflux python=3.11```
	 
	```conda activate thermoflux```

2. Now thermo_flux can be safely installed. Clone the thermo_flux repository and navigate to the thermo_flux directory:
  
	```git clone https://github.com/tednsmith/thermo_flux.git```

	```cd thermo_flux```

3. For development use the -e flag (for an editable install), navigate to where you cloned the thermo_flux directory and run:

	```python -m pip install -e .``` 
	
4. For thermodynamic optimisations, ensure a Gurobi license is installed correctly [(free for academics)](https://www.gurobi.com/academia/academic-program-and-licenses/):

	```conda install -c gurobi gurobi```
	
