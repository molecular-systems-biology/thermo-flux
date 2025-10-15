print('start')

import numpy as np
import pandas as pd
import cobra 
import os
from os import listdir
import gurobipy

print('modules imported')

arr = os.listdir(os.curdir) #directory you have downloaded .sol files


sol_files = [file for file in arr if file.endswith('.sol')]

time_files = [file for file in arr if file.startswith('res_plant')]

model = cobra.io.read_sbml_model("leaf_model.xml")
print(model)

rxn_ids = [rxn.id for rxn in model.reactions]
met_ids = [met.id for met in model.metabolites]

#times
def get_sol_progress(result_in):
    results = []

    with open(result_in, 'r') as f:
        for line in f:
            line = line.strip()
            data = line.split('_')
            if len(data) > 1:
                seed = data[1]
                time = data[2]
                objval = data[3]

                results.append([seed, time, objval])


    df = pd.DataFrame(results)
    df.columns = ['seed','time','objval']
    df['time'] = df['time'].astype(float)/60/60
    df['objval'] = df['objval'].astype(float) 

    return df

time_dfs = {}
for time_file in time_files:
    time_df = get_sol_progress(time_file)
    time_df = time_df.sort_values(['seed','time'])

    model_n = time_file.split('_')[2]

    time_dfs[model_n] = time_df


def get_results_sol(sol_in,model):

    with open(sol_in, 'r') as f: 
        v  =[]
        drG = []
        drG0 = []
        drG_error = []
        drG_conc = []
        lnc = []
        dfG_biomass = []
        res_dfG_biomass = []
        qm = []
        for line in f:
            line.strip()
            if len(line.split()) ==2:
                (names, xn) = line.split()
                xn = float(xn)
                if "flux" in names:
                    v.append(xn)

                if "drG[" in names:
                    drG.append(xn)

               # if "drG0" in names:
                #    drG0.append(xn)

                if "drG_error[" in names:
                    drG_error.append(xn)

                if "drG_conc[" in names:
                    drG_conc.append(xn)

                if "metabolite_log_concentration" in names:
                    lnc.append(xn)
                    
                if 'Gdiss_ex_conc' in names:
                    Gdiss_ex_conc = xn
                    
                if 'g_1' in names:
                    g_1 = xn

                if 'covariance_degrees_of_freedom' in names:
                    qm.append(xn)
                    
               
        drG_df = pd.DataFrame(drG, index=rxn_ids, columns=['drG'])

        Gdiss = np.array(v)*np.array(drG)

        drG_df['flux'] = v
        drG_df['Gdiss'] = Gdiss

        #drG_df['drG0'] = drG0
        drG_df['drG_error'] = drG_error
        drG_df['drG_conc'] = drG_conc
        drG_df['reaction'] = [rxn.reaction for rxn in model.reactions]

        drG_df.index.name = 'Reaction'

        lnc_df = pd.DataFrame(lnc, index=met_ids, columns=['lnc'])

        qm_df = pd.DataFrame(qm,columns=['qm'])

        
     #   g2 = g_1+Gdiss_ex_conc

        return drG_df, lnc_df, qm_df


#identify final solutions to reduce data size
max_sol = {}
for sol_file in sol_files:
    print(sol_file.split('_'))
    photons = sol_file.split('_')[1]
    seed = sol_file.split('_')[2]
    sol = sol_file.split('_')[3][:-4]
    print(photons, seed, sol)
    
    if sol != '':
        if (photons,seed) not in max_sol:
            max_sol[(photons,seed)] = int(sol)
        else:
            if int(sol) > max_sol[(photons,seed)]:
                max_sol[(photons,seed)] = int(sol)


drG_dfs = {}
lnc_dfs = {}
qm_dfs = {}
i = 0
for sol_file in sol_files:
    print(sol_file)
    i +=1
    print(sol_file.split('_'))
    photons = sol_file.split('_')[1]
    #drGlim = sol_file.split('_')[2]
    #Qlim = sol_file.split('_')[3]
    seed = sol_file.split('_')[2]
    sol = sol_file.split('_')[3][:-4]
    print(photons, seed, sol)

    
    if sol != '':
        sol = int(sol)
        time_df = time_dfs[photons]
        times = time_df.loc[time_df['seed']==seed]['time'].values
        sol_time = times[sol]
        rel_time = times[sol]/max(times)

        if sol == max_sol[(photons, seed)]:

            #print(model_number, sample, seed, sol)
            if not os.stat(sol_file).st_size == 0:
                drG_df, lnc, qm = get_results_sol(sol_file, model)
                print(len(drG_df))
                drG_dfs[photons, seed, sol, sol_time, rel_time] = drG_df
                lnc_dfs[photons, seed, sol, sol_time, rel_time] = lnc
                qm_dfs[photons, seed, sol, sol_time, rel_time] = qm

     

    print('sol files processed '+str(i)+' of '+str(len(sol_files)))


total_frame = pd.concat(list(drG_dfs.values()), keys=list(drG_dfs.keys()))
total_frame.index = total_frame.index.set_names(['photons', 'seed', 'sol', 'sol_time', 'rel_time','rxn'])
total_frame.reset_index(inplace=True)
total_frame.to_pickle('total_frame_16.pkl')

lnc_frame = pd.concat(list(lnc_dfs.values()), keys=list(lnc_dfs.keys()))
lnc_frame.reset_index(inplace=True)
lnc_frame.columns = ['photons', 'seed', 'sol', 'sol_time', 'rel_time','met','lnc']
lnc_frame.to_pickle('lnc_frame_16.pkl')

qm_frame = pd.concat(list(qm_dfs.values()), keys=list(qm_dfs.keys()))
qm_frame.reset_index(inplace=True)
qm_frame.columns = ['photons', 'seed', 'sol', 'sol_time', 'rel_time','dof','qm']
qm_frame.to_pickle('qm_frame_16.pkl')
