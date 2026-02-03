Infeasible gurobi model : some tips
===================================

Simulations using a thermodynamic–stoichiometric model rely on the **Gurobi** solver. Once the stoichiometric model has been converted into a thermodynamic model, curated, and fully balanced, the Gurobi model object can be initialized (see Step 7 of the workflow).

At this stage, the thermodynamic simulations may still be infeasible or return a zero flux for the objective. The following guidelines and tips can help you diagnose and resolve these issues.


1. Computing an Irreducible Inconsistent Subsystem (IIS).
**********************************************************
An IIS is a subset of the constraints and variable bounds with the following properties:
*It is still infeasible, and
*If a single constraint or bound is removed, the subsystem becomes feasible.
Most of the debugging can be done using the built-in function that identifies constraints responsible for infeasibility:

```python
thermo_flux.solver.gurobi.compute_IIS(tmodel)
```

For more details, see the official Gurobi documentation:
[https://support.gurobi.com/hc/en-us/articles/15656630439441-How-do-I-use-compute-IIS-to-find-a-subset-of-constraints-that-are-causing-model-infeasibility](https://support.gurobi.com/hc/en-us/articles/15656630439441-How-do-I-use-compute-IIS-to-find-a-subset-of-constraints-that-are-causing-model-infeasibility)

Note : an infeasible model is different than zero growth / zero objective flux.To better diagnose the infeasibility issue, you may temporarily **force biomass production** by setting a lower bound on the biomass reaction. This can intentionally render the model infeasible, making it easier for the IIS to identify problematic constraints.


2. Check reaction directions and flux bounds
*********************************************
Reaction directions that are inconsistent with thermodynamic laws (and metabolome data if applied) are a common cause of infeasibility.

* Inspect any **predefined reaction directions**. If a reaction is forced in a direction that is thermodynamically infeasible, the optimization will fail.
* As a first test, try allowing all reactions to be **reversible**.
* Verify the **model objective and bounds**:

  * Is the biomass lower bound achievable?
  * Are flux bounds (e.g., from physiological data) too restrictive?


3. Identifying the problematic reactions by not applying the second law constraint
***********************************************************************************
You can identify problematic reactions by temporarily ignoring the second-law constraints:

1. Set `ignore_snd = True` for all reactions and test whether the model becomes feasible.
2. Gradually set `ignore_snd = False` for individual reactions or groups of reactions to isolate the source of infeasibility. To narrow down the reactions which are causing the issue, start by ignoring the second law for the following categories:

* **proton or charge transport reactions** between compartments
* **all metabolite transport reactions**
* **complex reactions** that generated warnings during the balancing step

Once the model becomes feasible, set `ignore_snd=True` back gradually until the problematic reactions are identified.


