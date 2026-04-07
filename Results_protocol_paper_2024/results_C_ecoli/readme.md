Folder organization :
- Figures used in the manuscript are in Figures/
- plotting_datafolder/ contains processed datasets that can be used to recreate the figures.
- Data used to define and constrain the thermodynamic model (physiological data, metabolome data, files to annotate the model) are in datafiles/.
- ecoli_{condition}_metabolome{bool}.csv are the flux bounds resulting from thermodynamically constrained flux variability analysis

Notebooks : 
- 1a_tfva_ecoli_Glc_AC : thermodynamically constrained flux variability analysis for Glucose and Acetate condition
- 1b_tfba_ecoli : thermodynamically constrained FBA in low-glucose conditions (Vemuri 2007)
- 2_ecoli_analysis_directions : Comparing  thermodynamically constrained FVA reaction directions with predefined directions from the BiGG reconstruction iJR904
- 3_ecoli_tfba_gurscan : compare predictions from 1b (thermodynamics+ metabolome data) vs standard predefined direction fba and vs loopless FBA with fully reversible model

