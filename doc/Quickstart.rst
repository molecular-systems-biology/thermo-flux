Quick start tutorial
====================
As a first example, we will show how to start from a stoichiometric model and build the thermodynamic model using a toy model (e.coli core metabolism).
For more details, please have a look at the examples/ecoli/textbook/Ecoli_example_textbook.ipynb.

Building the thermo-model object
********************************

Load a Cobra model using the cobra loading function.
``model = cobra.io.load_model('e_coli_core')``

Create the thermodynamic model object from the cobra model 
``tmodel = ThermoModel(model)``

Step 1 : Defining the physical and biochemical parameters
*********************************************************

We use the pint package to define units for parameters using the Q_(value, unit) function
It is important to use the same unit constructor as the equilibrator_api (from equilibrator_api import  Q_)
 
``tmodel.pH = {"c": Q_(7.6), "e": Q_(7)}`` #pH

``tmodel.I = {"c": Q_(0.25,'M'), "e": Q_(0.25,'M')}`` #ionic stength

Steps 2-6 : Calculating Thermodynamic parameters
************************************************

Calculation of all thermodynamic parameters in a model is combined into a single function ``tmodel.update_thermo_info``. 

This function identifies all compounds based on their annotations (Step 2), calculates metabolites formation energies (Step 3), and calculates the Gibbs free energy of reactions (Step 6).

``tmodel.update_thermo_info(fit_unknown_dfG0=True, round_dp=2)``

Step 5. Executing pH-dependent charge and proton balancing 
****************************************************************

Reactions can now be automatically proton and charge balanced

``for rxn in tmodel.reactions: thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False, round_dp=2)``

Step 7. Establishing the thermodynamic-stoichiometric solution space 
***************************************************************************

We can now add the TFBA variables to the model :

``tmodel.m = None #reset the gurobi model object in case you're re-running this cell ``

``tmodel.add_TFBA_variables(gdiss_constraint = True, sigmac_limit = 12, error_type = 'covariance')`` 

This has now created a Gurobi model object accessible via ``tmodel.m`` .The gurobi model object can be interacted with directly via the gurobypy api. Extensive documentation can be found here: https://www.gurobi.com/documentation/10/refman/py_model.html

Set a 15s time limit for the solver

``tmodel.m.Params.TimeLimit = 15``

Start the optimization 

``tmodel.m.optimize()``

``Time limit reached. Best objective 0.545 , best bound 0.521, gap 1.36%``

Now the model is fully parameterized we can run FBA style optimizations with thermodynamic constraints.

