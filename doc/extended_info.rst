Extended Info
=============

Introduction
************
This page provides detailed technical guidance on the steps required to build a complete and functional model with Thermo-Flux. Going beyond the quickstart and the two tutorials (*e.coli* core model and iMM904 model), here we present examples of code implementation that are relevant to augment your own stoichiometric model with thermodynamic constraints. For the background information and the modelling assumptions, we encourage you to refer to the protocol paper.

Step 1: Definition of physical and biochemical parameters
*********************************************************
Thermo-Flux models rely on a well-defined thermodynamic environment. This includes pH, ionic strength, temperature, and compartment-specific membrane potentials.

Membrane Potentials
-------------------
Membrane potential difference between compartment c and e is defined as
:math: `\Phi_{ce} = \Phi_c - \Phi_e`
with ``Phi_c`` representing the potential in compartment ``c``.

Table 1 : Physical and biochemical parameters, their units and examples of how to define them in ‘Thermo-Flux’. 

+-------------------------------+-----------------------+-----------------------------------------+
| Parameter                     | Unit                  | Code example                            |
+===============================+=======================+=========================================+
| Potential hydrogen (pH)       | Dimensionless         | .. code-block:: python                  |
|                               |                       |                                         |
|                               |                       |    model.pH = {"c": Q_(7.6), "e": Q_(7)}|
+-------------------------------+-----------------------+-----------------------------------------+
| Ionic strength (I)            | mol L-1 (M)           | .. code-block:: python                  |
|                               |                       |                                         |
|                               |                       |    model.I = {"c": Q_(0.25, 'M'),       |
|                               |                       |               "e": Q_(0.25, 'M')}       |
+-------------------------------+-----------------------+-----------------------------------------+
| Membrane potential differences| Volt (V)              | .. code-block:: python                  |
| (phi)                         |                       |                                         |
|                               |                       |    model.phi = {'ce': Q_(0.15, 'V')}    |
+-------------------------------+-----------------------+-----------------------------------------+
| Potential magnesium (pMg)     | Dimensionless         | .. code-block:: python                  |
|                               |                       |                                         |
|                               |                       |    model.pMg = {'c': Q_(3), 'e': Q_(3)} |
+-------------------------------+-----------------------+-----------------------------------------+
| Temperature (T)               | Kelvin (K)            | .. code-block:: python                  |
|                               |                       |                                         |
|                               |                       |    model.T = Q_(298, 'Kelvin')          |
+-------------------------------+-----------------------+-----------------------------------------+



Step 2: Definition of metabolites and chemical species
******************************************************
Connection with eQuilibrator
----------------------------
``Thermo-Flux`` automatically recognises the metabolite identifiers and links to the eQuilibrator compound retrieved with ``metabolite.compound``.

Annotations can be user-defined or updated in the attribute ``annotation`` of the metabolite class e.g.:

::

    metabolite.annotation = {'CHEBI':'11111', 'kegg':'C00000'}

Metabolite species
------------------

In ``Thermo-Flux``, the average charge, average number of protons, and magnesium ions are returned by the function ``metabolite.average_charge_protons()``, which first interrogates the eQuilibrator compound and then uses the physical parameters defined in Step 1 to return the condition-specific metabolite information.

To facilitate the understanding of these average calculations, this function also returns the information on each species’ abundance and their charge, number of protons and magnesium ion (Figure 2c).


Definition of metabolites with non-decomposable or unknown structures
---------------------------------------------------------------------
The formula and charge for these metabolites should be defined using the COBRApy attributes with ``metabolite.formula`` and ``metabolite.charge``.
::
     model.metabolite.formula = 'C6H12O6'
     model.metabolite.charge = -1

Local cache to access eQuilibrator compounds
--------------------------------------------
When ``Thermo-Flux`` queries an eQuilibrator compound for the first time, eQuilibrator will require downloading the latest up-to-date database of eQuilibrator compounds. This local cache is named ``compound.sqlite`` and integrates native functions to retrieve compounds or manually add compounds (see `eQuilibrator local cache <https://equilibrator.readthedocs.io/en/latest/local_cache.html>`_).

Step 3: Calculation of Gibbs formation energies
***********************************************

The function ``model.update_thermo_info()`` will automatically calculate the required parameters based on the defined physiochemical conditions (Step 1) and the metabolites of the model will now have a defined transformed Gibbs formation energy (``\Delta_f G^\prime``) and an average charge and number of protons.

.. rubric:: Box 1: Additional considerations for Gibbs energy of formation calculation
--------------------------------------------------------------------------------------

**Uncertainty**

Different default uncertainty can be specified with ``model.rmse_inf = Q_(3000, 'kJ/mol')``.

We can also estimate a non-zero Gibbs formation energy for metabolites with non-decomposable or unknown structures (see supplementary section “metabolites with unknown formation energy”). This is implemented by the ``fit_unknown_dfG0=True`` argument when estimating Gibbs formation energies.

**Redox**

In ``Thermo-Flux`` a formation energy and a standard error can be explicitly defined, and the ``redox`` attribute set to true to ensure the formation energy is not automatically recalculated.
For example, the midpoint potential of cytochrome C is 250 mV (Lennarz & Lane, 2013). Applying the Nernst equation at equilibrium yields ∆𝑟𝐺° ′=− 1 × 96.5 × 0.250 = −24.125 𝑘𝐽𝑚𝑜𝑙−1. If the oxidised and reduced cytochrome C always appear as a pair in reactions, the formation energy of the reduced form can be defined as −24.125/2 = −12.05 𝑘𝐽𝑚𝑜𝑙−1 and the oxidised form as 12.05 𝑘𝐽.𝑚𝑜𝑙−1. e.g.:

::

    cyt_c_red_c.dfG0prime() = Q_(-12.05, 'kJ/mol')
    cyt_c_red_c.redox = True
    cyt_c_red_c.dfG_SE = Q_(0, 'kJ/mol')


**Biomass**

In ``Thermo-Flux``, the function ``thermo_flux.tools.drg_tools.dfGbm()`` returns the biomass formation energy given a specified empirical formula of biomass and can be used to explicitly define the biomass formation energy, e.g.:

::

    dfGbm = thermo_flux.tools.drg_tools.dfGbm(H=1.613, O=0.557, N=0.158, P=0.012,
                                              S=0.003, K=0.022, Mg=0.003, Ca=0.001,
                                              units='kJ/g')
    model.metabolites.biomass.dfGprime() = dfGbm
    model.metabolites.biomass.biomass = True
    model.metabolites.biomass.dfG_SE = 0

Care must be taken when defining the units of the biomass formation energy. To maintain consistency with cellular metabolic reactions, the unit of the formation energy is entered as ``kJ mol^{-1}`` like other metabolites, but in reality it is in ``J gDW^{-1}``. This is because the biomass equation converts mmol of metabolites into gDW of biomass whereas formation energies are defined as ``kJ mol^{-1}``.

Biomass formation energy is made dependent on the pH of the biomass metabolite’s compartment when transformed based on the number of hydrogen atoms of which it is formed. It is done automatically when building a ``Thermo-Flux`` model if ``model.update_biomass_dfG0`` is set to True.

Step 4: Delineation of transporter characteristics
**************************************************

For each transport reaction, ‘Thermo-Flux’ will automatically determine the transported metabolite, the transported charge, and the transported protons, depending on the defined physiological parameters of the compartments and the reaction stoichiometry. Additional transported protons can be achieved by altering the reaction stoichiometry (Figure 3b). Alternatively, additional transported protons can be defined using ``reaction.transported_h()``, which represents additional protons transported by a reaction, e.g.:

::

    reaction.transported_h = {'e': -1, 'c': 1}

to define an additional proton moving from the extracellular (``e``) compartment to the cytosol (``c``).

.. rubric:: Box 2: additional considerations for transport reactions
--------------------------------------------------------------------

**Adding transporter variants**
Additionally, in case of transport processes, for which at the given pH value no charge-neutral transport variant exists, we suggest introducing an additional transport reaction, in which protons balancing the charge are co-translocated together with the respective species, i.e., adding a proton symporter or antiporter. This additional transport variant ensures that for every metabolite, a transport variant exists that does not translocate net charge.

Addition of transporter variants can automatically be achieved with the function ``reaction.add_transporter_variants()``, which identifies the species transported in the original reaction and adds variants to represent the transport of all alternative species.

For example, a model may contain a reaction for phosphate transport, ``pi_e -> pi_c``. At pH 5, this ion exists entirely in the ``H_2PO_4^-`` form with a charge of -1 (Figure 3a). Therefore, all the major species of the latter ion are already represented but a charge-neutral transporter does not exist. A proton coupled reaction of ``pi_e + H_e -> pi_c + H_c`` is automatically added to the model (Figure 3b).

**Transporters with simultaneous chemical transformation of the transported metabolite**
Some transport reactions involve chemical transformation of the transported metabolite, e.g., phosphotransferase system (PTS) sugar transporters which phosphorylate sugars during transport (McCoy et al., 2015). In this case it is not possible to automatically determine the specific metabolite that is transported, as it does not appear as both a substrate and product of the reaction. Therefore, it is necessary to manually specify the transported metabolite using e.g.:

::

    reaction.transported_mets = {Glc_e: -1}

to represent extracellular glucose as the metabolite that is transported across the membrane.

**Reporting**

By setting the argument ``report`` to True, the function ``model.update_thermo_info()`` can provide a reporting table as a pandas DataFrame, with information on the stoichiometry, balancing status, and transported metabolites/charge/protons of each reaction. In this table, reactions that require inspection by the user will appear in the top rows.

**Ambiguous proton or ion transporters**

It is important to distinguish between free protons that are transported as part of the transport mechanism (e.g. in proton symporters) and protons which are bound/released from metabolites as part of a chemical reaction.
In general, this is automatically determined but in some cases is ambiguous. Ambiguous reactions are highlighted to the user for manual curation. Curation consists of specifying manually the number of transported free protons or ions, e.g., ``reaction.transported_h = {'e': -1, 'c': 1}`` to represent the transport of one proton from the extracellular to cytosolic compartment.

As an example, the reaction of mitochondrial Complex II

::

    Ubiquinone-8_c + succinate_m <=> fumarate_m + Ubiquinol-8_c

would need the user to specify:

::

    tmodel.reactions.ComplexII.transported_h = {'m': -2.0, 'c': 2.0}

as two protons are moved from the mitochondria to the cytosol and are subsequently taken up by the protonation of Ubiquinone-8 into Ubiquinol-8.

Step 5: pH-dependent charge and proton balancing
************************************************

Non-transport reactions
-----------------------

The function ``reaction_balance()`` can be used to automatically balance the protons in a reaction based on the compartment conditions with the option to also balance magnesium ions if desired.

In the example of ATP hydrolysis, 0.7 protons will be added to have an equal number of protons and charge on both sides of the reaction (protons are positively charged and therefore charge balance is also maintained).

Transport reactions
-------------------

To balance transport reactions, ``Thermo-Flux`` first identifies the most abundant species (using ``metabolite.major_microspecies`` automatically), then considers it as being transported. The balancing then occurs by comparing what is in the inner compartment, what is being transported and what will be in the outer compartment.

Magnesium ions
--------------
Analogously to protons, Mg2+ ions can also be balanced, and this option is available to the user by setting ``balance_mg=True``.

Step 6: Calculation of Gibbs energy of reactions
************************************************

To calculate the standard reaction energy of all reactions in the model, the function ``model.update_thermo_info()`` can be used. Once it has been run, the standard reaction energy and the standard transformed reaction energy (calculated using standard transformed formation energies) can be retrieved for each reaction with ``reaction.drG0`` and ``reaction.drG0prime``, respectively.

Step 7: Establishment of the thermodynamic-stoichiometric solution space
************************************************************************

Metabolite concentration bounds
-------------------------------

In practice metabolite concentration bounds are defined by setting the ``lower_bound`` and ``upper_bound`` attributes and a user defined unit e.g.:

::

    metabolite.lower_bound = Q_(10, 'µM')

The concentration values will then be automatically converted to mol/L before applying thermodynamic constraints.

The function ``model.add_TFBA_variables()`` sets up a thermodynamic FBA optimisation problem using the Gurobi optimiser that can be optimised using ``model.m.optimize()``. Implementation of the constraints in the linear program is detailed in the methods (see: implementing conditional constraints in a linear program).

.. rubric:: Box 3: additional considerations for the formulation of the thermodynamic/stoichiometric solution space
-------------------------------------------------------------------------------------------------------------------

**Compartmented metabolite concentrations and whole cell concentrations**

The function ``model.total_cell_conc()`` will add whole cell metabolite concentration constraints on the compartmented metabolic concentrations, based on whole cell metabolite data and the relative compartment volumes which must be provided as an input to the function, respectively as a pandas DataFrame and a Python dictionary.

**Relaxing the second law constraint**

The user can relax the second law constraint for any specific reaction by setting ``reaction.ignore_snd = True``.

**Ignoring metabolite concentrations**

The concentration of pseudo metabolites that are often added to stoichiometric models as a convenient way to add constraints should also be ignored by setting

::

    metabolite.ignore_conc = True

**Variability analysis**

In ``Thermo-Flux`` variability analysis is implemented with the function ``solver.gurobi.variability_analysis()``, which sets the optimization problem for any variables provided as an argument to the function. Specifically, the function uses the Gurobi multi-scenario optimization feature, with two scenarios for each variable (one minimizes the variable and the other maximizes it). The results are retrieved with ``solver.gurobi.variability_results()`` and both functions can still be used if the optimization is solved using a high-performance computing (HPC) cluster.

Step 8: Regression: fitting models to experimental data
*******************************************************

The function ``model.regression()`` can be used to add regression constraints and objectives to the previously constructed thermodynamic FBA problem. Data can be provided for any flux or metabolite concentration, in the pandas DataFrame format.

The Dataframe for the fluxes and the metabolite data needs to be in the following format :  

+--------------+-------------+-------+-----+
| condition    |   rxn/met   | mean  | sd  |
+==============+=============+=======+=====+
| condition 1  | rxn/met  A  |  XXX  | YYY |
+              +-------------+-------+-----+
|              | rxn/met  B  |  ZZZ  | WWW |
+--------------+-------------+-------+-----+

Note the pandas.MultiIndex (condition,rxn/met).

.. rubric:: Box 4: additional considerations for regressions
------------------------------------------------------------

**Model starting points**

The function ``thermo_flux.solver.gurobi.model_start`` has been built to allow MIP start from only non-computed values and reduce the probability of multiplying numerical issues between them. This function can even enable the start from a set of specific variables which are known to not cause numerical issues (for example, starting from only metabolite concentrations). The user can provide starting points in either ``.sol`` or ``.mst`` format:

::

    thermo_flux.solver.gurobi.model_start(tmodel, 'filename.sol',
                                         ignore_vars=['all'],
                                         fix_vars=['qm','ln_conc'],
                                         fix='start')

**Multiple starts with different random seeds**

As Gurobi is using a branch-and-cut approach to solve the MILP problem, it can sometimes face performance variability issues. An effective way of tackling this problem is to run several optimizations with different values of the seed parameter ``GRBmodel.params.Seed``.
