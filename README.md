# Thermo-Flux
Tools for flux balance anlysis with thermodynamic constraints. 

Full documentation on ReadTheDocs (https://thermo-flux.readthedocs.io/en/latest/) and accompanying protocol at [doi:10.1101/2025.11.20.689566](https://doi.org/10.1101/2025.11.20.689566).

## Installation

Requirements
- Python >= 3.11
- [Gurobi 11.0](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-)
	
1. To avoid dependency conflicts it is reccomended to use a python environment e.g:

	```conda create -n thermoflux python=3.11```
	 
	```conda activate thermoflux```

2. Now thermo_flux can be safely installed. Clone the thermo_flux repository and navigate to the thermo_flux directory:
  
	```git clone https://github.com/molecular-systems-biology/thermo-flux```

	```cd thermo_flux```

3. For development use the -e flag (for an editable install), navigate to where you cloned the thermo_flux directory and run:

	```python -m pip install -e .``` 
	
4. For thermodynamic optimisations, ensure a Gurobi license is installed correctly [(free for academics)](https://www.gurobi.com/academia/academic-program-and-licenses/):

	```conda install -c gurobi gurobi```
	
## Examples 

Example usage notebooks can be found in the examples directory. 

## Cite us

If you use Thermo-Flux please cite our protocol paper at: 

Smith EN, Fargier N, Pedro J & Heinemann M (2025) *Thermo-Flux: generation and analysis of comprehensive thermodynamic-stoichiometric metabolic network models.* 2025.11.20.689566 [doi:10.1101/2025.11.20.689566](https://doi.org/10.1101/2025.11.20.689566) [PREPRINT]
  
Thermo-Flux relies on eQuilibrator for underlying thermodynamic calculations so please also cite the eQuilibrator database: 

M. E. Beber, M. G. Gollub, D. Mozaffari, K. M. Shebek, A. I. Flamholz, R. Milo, and E. Noor, *eQuilibrator 3.0: a database solution for thermodynamic constant estimation* Nucleic Acids Research (2021), [DOI:10.1093/nar/gkab1106](http://dx.doi.org/10.1093/nar/gkab1106)

## Acknowledgments
Thanks to Yusuke Himeoka, Moritz Beber, Elad Noor and Mattia Gollub for helpful discussions and advice. 
