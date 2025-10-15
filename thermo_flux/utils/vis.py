"""visualisation utilities"""

import sys
from cobra.util.array import create_stoichiometric_matrix
import pandas as pd

def progressbar(it, prefix="", size=60, out=sys.stdout, item_label_attribute = None): 
    """Displays a progress bar when running loops
    from https://stackoverflow.com/questions/3160699/python-progress-bar""" 

    count = len(it)

    def show(j, item_label):
        x = int(size*j/count)
        print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count} {item_label}",
              end='                 \r', file=out, flush=True)
        
    show(0, '')
    for i, item in enumerate(it):
        if item_label_attribute is not None:
            item_label = item.__getattribute__(item_label_attribute)
        else:
            item_label = None
        show(i+1, item_label)
        yield item

    print("\n", flush=True, file=out)


def compare_met(model_1, model_2, met, model_names):
    #create S 
    S1 = create_stoichiometric_matrix(model_1) 
    S2 = create_stoichiometric_matrix(model_2) 
    
    df1 = pd.DataFrame(S1, columns=[rxn.id for rxn in model_1.reactions], index = [met.id for met in model_1.metabolites])
    df2 = pd.DataFrame(S2, columns=[rxn.id for rxn in model_2.reactions], index = [met.id for met in model_2.metabolites])

    #make df 
    df = pd.concat([df1.loc[met],df2.loc[met]], axis=1)
    df.columns = model_names

    df['diff'] = df.iloc[:,0] - df.iloc[:,1]
    

    return(df)
    