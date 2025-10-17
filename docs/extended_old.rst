step 1 :  Definition of physical and biochemical parameters
***********************************************************
For definition of membrane potentials, ‘ce’ represents the voltage between compartment ‘c’ and compartment 'e’ defined as Phic - Phie.
Table 1 (see old version)


Step 2: Definition of metabolites and chemical species 
******************************************************
‘Thermo-Flux’ automatically recognises the metabolite identifiers and links to the eQuilibrator compound retrieved with metabolite.compound.
 Annotations can be user-defined or updated in the attribute annotation of the metabolite class e.g. 

metabolite.annotation = {'CHEBI’:’11111’, ‘kegg’:’C00000’} 

 In ‘Thermo-Flux’, the average charge, average number of protons, and magnesium ions are returned by the function metabolite.average_charge_protons() which first interrogates the eQuilibrator compound and then uses the physical parameters defined in Step 1 to return the condition specific metabolite information. 
To facilitate the understanding of these latter average calculations, this function also returns the information on each species’ abundance and their charge, number of protons and magnesium ion (Figure 2c).   

Box 1 : additionnal consideration for metabolite definition
-----------------------------------------------------------

Definition of metabolites with non-decomposable or unknown structures
=====================================================================
The formula and charge for these metabolites should be defined using the COBRApy attributes with metabolite.formula and metabolites.charge. 

Local cache to access eQuilibrator compounds 
============================================
When 'Thermo-Flux' queries an eQuilibrator compound for the first time, eQuilibrator will require downloading the latest up-to-date database of eQuilibrator compounds. This local cache is named compound.sqlite and integrates native functions to retrieve compounds or manually add compounds (https://equilibrator.readthedocs.io/en/latest/local_cache.html).  

Step 3. Calculation of Gibbs formation energies  
***********************************************
The function model.update_thermo_info() will automatically calculate the required parameters based on the defined physiochemical conditions (step 1) and the metabolites of the model will now have a defined transformed Gibbs formation energy ([Equation]) and an average charge and number of protons.

box 2 : additional consideration for the calculation of Gibbs formation energies  
-----------------------------------------------------------
uncertainty
(different default uncertainty can be specified with model.rmse_inf = Q_(3000, 'kJ/mol')),
We can also estimate a non-zero Gibbs formation energy for metabolites with non-decomposable or unknown structures (see supplementary section “metabolites with unknown formation energy”). This is implemented by the fit_unknown_dfG0=True argument when estimating Gibbs formation energies.   

redox
In ‘Thermo-Flux’ a formation energy and a standard error can be explicitly defined, and the “redox” attribute set to true to ensure the formation energy is not automatically recalculated e.g. 
cyt_c_red_c.dfG0prime() = Q_(-12.05,’kJ/mol') 
cyt_c_red_c.redox = True 
cyt_c_red_c.dfG_SE = Q_(0,’kJ/mol')  

biomass:
In ‘Thermo-Flux’ the function: thermo_flux.tools.drg_tools.dfGbm() returns the [Equation] given a specified empirical formula of biomass and can be used to explicitly define the biomass formation energy e.g. 
dfGbm = thermo_flux.tools.drg_tools.dfGbm(H=1.613,O=0.557,N=0.158, P=0.012, S=0.003, K = 0.022, Mg = 0.003, Ca= 0.001, units = kJ/g) 
model.metabolites.biomass.dfGprime() = dfGbm  
model.metabolites.biomass.biomass = True 
model.metabolites.biomass.dfG_SE = 0  

Care must be taken when defining the units of the biomass formation energy. To maintain consistency with cellular metabolic reactions, the unit of the formation energy is entered as kJmol-1 like other metabolites, but in reality it is in JgDW-1 . This is because the biomass equation converts mmol of metabolites into gDW of biomass whereas formation energies are defined as kJmol-1 .

 Biomass formation energy is made dependent on the pH of the biomass metabolite’s compartment when transformed based on the number of hydrogen atoms of which it is forme. It is done automatically when building a ‘Thermo-Flux’ model if model.update_biomass_dfG0 is set to True.


step 4.
Defining additional transported protons can be achieved by altering the reaction stoichiometry (Figure 3b). Alternatively additional transported protons can be defined using reaction.transported_h() which represents additional protons transported by a reaction e.g. reaction.transported_h = {'e' : -1, 'c' : 1}, to define an additional proton moving from the extracellular (e) compartment to the cytosol (c).

box3 : 
transporter variants :
Additionally, in case of transport processes, for which at the given pH value no charge neutral transport variant exists, we suggest introducing an additional transport reaction, in which protons balancing the charge are co-translocated together with the respective species, i.e. adding a proton symporter or antiporter. This additional transport variant ensures that for every metabolite, a transport variant exists that does not translocate net charge. Addition of transporter variants can automatically be achieved with the function reaction.add_transporter_variants() which identifies the species transported in the original reaction and adds variants to represent the transport of all alternative species. For example, a model may contain a reaction for phosphate transport, pi_e -> pi_c. At pH 5, this ion exists entirely in the H2PO4- form with a charge of -1 (Figure 3a). Therefore, all the major species of the latter ion are already represented but a charge neutral transporter does not exist. A proton coupled reaction of pi_e  + H_e -> pi_c + H_c is automatically added to the model (Figure 3b). 

Some transport reactions involve chemical transformation of the transported metabolite e.g. phosphotransferase system  (PTS) sugar transporters which phosphorylate sugars during transport (McCoy et al, 2015). In this case it is not possible to automatically determine the specific metabolite that is transported, as it does not appear as both a substrate and product of the reaction. Therefore, it is necessary to manually specify the transported metabolite. using e.g.  
reaction.transported_mets = {Glc_e: -1} to represent extracellular glucose as the metabolite that is transported across the membrane. 

By setting the argument ‘report’ to True, the function model.update_thermo_info() can provide a reporting table as a pandas DataFrame, with information on the stoichiometry, balancing status and transported metabolites/charge/protons of each reaction. In this table, reactions that require inspection by the user will appear in the top rows. 

 In general, this is automatically determined but in some cases is ambiguous. Ambiguous reactions are highlighted to the user for manual curation. Curation consists of specifying manually the number of transported free protons or ions. (e.g. reaction.transported_h = {'e' : -1, 'c' : 1} to represent the transport of one proton from the extracellular to cytosolic compartment). As an example, the reaction of mitochondrial Complex II " Ubiquinone-8_c + succinate_m <=> fumarate_m +   Ubiquinol-8_c” would need the user to specify "tmodel.reactions.ComplexII.transported_h = {'m': -2.0, 'c': 2.0}" as two protons are moved from the mitochondria to the cytosol and are subsequently taken up by the protonation of Ubiquinone-8  into Ubiquinol-8. 


Step 5. pH-dependent charge and proton balancing  
non transport
 The reaction_balance() function can be used to automatically balance the protons in a reaction based on the compartment conditions with the option to also balance magnesium ions if desired.
In the example of ATP hydrolysis, 0.7 protons will be added to have an equal number of protons and charge on both sides of the reaction (protons are positively charged and therefore charge balance is also maintained). 

Transport reactions 
‘Thermo-Flux’ first identifies the most abundant species (using metabolite.major_microspecies automatically)

Mg ions : 
 Analogously to protons, Mg2+ ions can also be balanced, and this option is available to the user  by setting balance_mg = True.

Step 6. Calculation of Gibbs energy of reactions  
To calculate the standard reaction energy of all reactions in the model, the function model.update_thermo_info() can be used. Once it has been run, the standard reaction energy and the standard transformed reaction energy (calculated using standard transformed formation energies) can be retrieved for each reaction with reaction.drG0 and reaction.drG0prime, respectively. 

Step 7. Establishment of the thermodynamic-stoichiometric solution space  
metabolite concentration bounds : 
 In practice metabolite concentration bounds are defined by setting the lower_upper and upper_bound attributes and a user defined unit e.g.   metabolite.lower_bound = Q_(10, ‘µM’). The concentration values will then be automatically converted to mol/L before applying thermodynamic constraints. 


The function model.add_TFBA_variables() sets up a thermodynamic FBA optimisation problem using the Gurobi optimiser that can be optimised using model.m.optimize(). Implementation of the constraints in the linear program is detailed in the methods see: implementing conditional constraints in a linear program. 

Box 4 : 
Compartmented metabolite concentrations  and whole cell concentrations
The function model.total_cell_conc() will add whole cell metabolite concentration constraints on the compartmented metabolic concentrations, based on whole cell metabolite data and the relative compartment volumes which must be provided as an input to the function, respectively as a pandas DataFrame and a python dictionary. 

relax 2nd law with The user can relax the second law constraint for any specific reaction by setting reaction.ignore_snd = True. 
ignoring metabolite concentrations : 
The concentration of pseudo metabolites that are often added to stoichiometric models as a convenient way to add constraints should also be ignored. using the function metabolite.ignore_conc = True.

Variability analysis 
In ‘Thermo-flux’ variability analysis is implemented with the function solver.gurobi.variability_analysis(), which sets the optimization problem for any variables provided as an argument to the function. Specifically, the function uses the Gurobi multi-scenario optimization feature, with two scenarios for each variable (one minimizes the variable and the other maximizes it). The results are retrieved with solver.gurobi.variability_results() and both functions can still be used if the optimization is solved using a high-performance computing (HPC) cluster. 

Step 8. Regression: fitting models to experimental data 
the function model.regression() can be used to add regression constraints and objectives to the previously constructed thermodynamic FBA problem. Data can be provided for any flux or metabolite concentration, in the pandas DataFrame format. , 

Box 5 : 
Model starting points
The function thermo_flux.solver.gurobi.model_start has been built to allow MIP start from only non-computed values and reduce the probability of multiplying numerical issues between them. This function can even enable the start from a set of specific variables which are known to not cause numerical issues (for example, starting from only metabolite concentrations). The user can provide starting points in either .sol or .mst format : thermo_flux.solver.gurobi.model_start(tmodel,'filename.sol’,ignore_vars=['all'],fix_vars=['qm','ln_conc'],fix='start').




