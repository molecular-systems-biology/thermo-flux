Infeasible gurobi model : some tips
===================================

Simulations using a thermodynamic-stoichiometric model rely on the gurobi solver. Once the stoichiometric model has been converted to a thermodynamic model, curated and fully balanced, the gurobi model object can be initialized as seen in step 7.
At this step, it could be that the thermodynamic simulations are not feasible yet or return a zero flux for the objective. 

To debug a bit more this model, here's a few tips to try :

-  Most of the debugging can be done using the built in function to find the constraints making the model infeasible "thermo_flux.solver.gurobi.compute_IIS(tmodel)". See https://support.gurobi.com/hc/en-us/articles/15656630439441-How-do-I-use-compute-IIS-to-find-a-subset-of-constraints-that-are-causing-model-infeasibility for more info.  
 Identify whether the model is infeasible or cannot produce biomass : a growth rate of 0 is different from an infeasible model.  It is possible to force the model to produce biomass by setting a lower bound on the biomass reaction. The model could be made infeasible like this, and the IIS could be able to run more easily. 

-  inspect any predefined reaction directions. If a reaction is forced in a direction that is now found to be thermodynamically infeasible the optimisation may be infeasible. Try allowing all reactions to be reversible in the first instance. In the same idea :  check model objective and bounds e.g. is there a biomass lower bound that cannot be achieved? 
-  test if a specific reaction is thermodynamically infeasible. For every reaction set ignore_snd = True and see if the model is now feasible. Then set ignore_snd = False for each reaction or for blocks of reactions to try and identify the problem.
To further narrow down and identify which reactions are causing trouble, ignore the second law for different kinds of transport reactions or complex reactions :
       -- protons/charge transport reactions between each compartments
      --  all metabolites transport
       -- complex reactions that give a warning in the balancing step



# Debugging Infeasible Gurobi Models

Simulations using a thermodynamic–stoichiometric model rely on the **Gurobi** solver. Once the stoichiometric model has been converted into a thermodynamic model, curated, and fully balanced, the Gurobi model object can be initialized (see Step 7 of the workflow).

At this stage, the thermodynamic simulations may still be infeasible or return a zero flux for the objective. The following guidelines and tips can help you diagnose and resolve these issues.

---

## 1. Identify the Source of Infeasibility

Most of the debugging can be performed using the built-in function that identifies constraints responsible for infeasibility:

```python
thermo_flux.solver.gurobi.compute_IIS(tmodel)
```

For more details, see the official Gurobi documentation:
[https://support.gurobi.com/hc/en-us/articles/15656630439441-How-do-I-use-compute-IIS-to-find-a-subset-of-constraints-that-are-causing-model-infeasibility](https://support.gurobi.com/hc/en-us/articles/15656630439441-How-do-I-use-compute-IIS-to-find-a-subset-of-constraints-that-are-causing-model-infeasibility)

### Important distinction

* **Infeasible model:** No solution satisfies all constraints.
* **Zero growth / zero objective flux:** The model is feasible but cannot produce biomass (or the selected objective).

To better diagnose the issue, you may temporarily **force biomass production** by setting a lower bound on the biomass reaction. This can intentionally render the model infeasible, making it easier for the IIS procedure to identify problematic constraints.

---

## 2. Check Reaction Directionality and Bounds

Incorrect reaction constraints are a common cause of infeasibility.

* Inspect any **predefined reaction directions**. If a reaction is forced in a direction that is thermodynamically infeasible, the optimization will fail.
* As a first test, try allowing all reactions to be **reversible**.
* Verify the **model objective and bounds**:

  * Is the biomass lower bound achievable?
  * Are flux bounds (from physiological data for instance) too restrictive?

---

## 3. Test for Thermodynamically Infeasible Reactions

You can identify problematic reactions by temporarily ignoring the second-law constraints:

1. Set `ignore_snd = True` for all reactions and test whether the model becomes feasible.
2. Gradually set `ignore_snd = False` for individual reactions or groups of reactions to isolate the source of infeasibility.

### Recommended order of investigation

To narrow down the issue more efficiently, start by ignoring the second law for the following categories:

* **Proton or charge transport reactions** between compartments
* **Metabolite transport reactions**
* **Complex reactions** that generated warnings during the balancing step

Once the model becomes feasible, reintroduce constraints incrementally until the offending reactions are identified.


