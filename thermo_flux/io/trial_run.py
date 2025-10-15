from excel import *



if __name__ == "__main__":
    model = create_cobra_model("yeast", model_excel="..\\..\\examples\\yeast\\yeast_v3.SV_ENS_03.xlsx", keggids_csv="..\\..\\examples\\yeast\\yeast_kegg_id.csv")
    print("-------------")
    print("Some checks:")
    print(model)
    print(model._parameters)
    print(model.metabolites[4].annotation)
    print([m.charge for m in model.metabolites if m.name.split("_")[0] in "ficytc"])
    print([m.formula for m in model.metabolites if m.name.split("_")[0] in "h"])
    print([m.formula for m in model.metabolites if m.name.split("_")[0] in "electron"])