# thermo_flux
Tools for thermodynamically constrained FBA [[1](#references)]


## Installation

Requirements
- Python >= 3.11
- [Gurobi 11.0](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-)


### Installation on windows 
	
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
	

## Examples 

Example usage notebooks can be found in the examples directory. 


## References 
1 .Niebel B, Leupold S, Heinemann M. An upper limit on Gibbs energy dissipation governs cellular metabolism. Nat Metab. 2019 Jan;1(1):125-132. doi: 10.1038/s42255-018-0006-7. Epub 2019 Jan 7. PMID: 32694810.
