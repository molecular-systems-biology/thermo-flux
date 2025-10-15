
from cobra.util.array import create_stoichiometric_matrix
import numpy as np
import pandas as pd
import scipy.stats
import gurobipy as gp
from gurobipy import GRB
from equilibrator_api import R, Q_


def add_TFBA_variables(tmodel, m, conds=[''], error_type='linear',
                       qnorm=1, alpha=0.95, epsilon=0.5, nullspace=None,
                       gdiss_constraint=False, sigmac_limit=12.3, split_v=False, big_M = False):
    
    """
    Generates gurobi model from thermo model. Contains constraints for FBA optimization
    """
    # RT term 
    RT = (R*tmodel.T.m_as('K')).m
    
    # max flux
    max_flux = abs(np.array([(rxn.upper_bound, rxn.lower_bound) for rxn in tmodel.reactions])).max()

    # max_drG
    max_drG = tmodel.max_drG.m_as('kJ/mol')
    
    # generate stoichiometric matrix
    S = create_stoichiometric_matrix(tmodel) 
    #remove very small stoichiometries to help numerics
    S[np.abs(S) < 1e-9] = 0
    Nc, Nr = S.shape  # number of compouds, number of reactions
    nconds = len(conds)  # number of conditions (for regression to multiple conditions)
    
    #array of drG0_prime and standard error
    #drGt0tr = np.array([round((rxn.drG0prime.m + rxn.drGtransport.m),5) for rxn in tmodel.reactions])
    drGt0tr = np.array([(rxn.drG0prime.m + rxn.drGtransport.m) for rxn in tmodel.reactions])
    drG_SE = np.array([rxn.drG_SE.m for rxn in tmodel.reactions])


    # Create variables
    b = m.addMVar(shape = (nconds,Nr), vtype=GRB.BINARY, name="reaction_directions")
    v = m.addMVar(lb=-500,ub=500,shape = (nconds,Nr), name="fluxes")
    drGp = m.addMVar(lb=0,ub=max_drG,shape = (nconds,Nr), name="drG_positive")
    drGn = m.addMVar(lb=0,ub=max_drG,shape = (nconds,Nr), name="drG_negative")
    ln_conc = m.addMVar(lb=np.log(1e-8),ub=np.log(1),shape = (nconds,Nc),name="ln_conc")
    drG = m.addMVar(lb=-max_drG, ub = max_drG,shape = (nconds,Nr),name="drG") #ToDo! update to default max drG value 
    drG_error =  m.addMVar(lb=-1e6, ub = 1e6,shape=Nr,name="drG_error")
    drG_conc =  m.addMVar(lb=-GRB.INFINITY, ub = GRB.INFINITY,shape = (nconds,Nr),name="drG_conc")

    if split_v:
        vp = m.addMVar(lb=0,ub=max_flux,shape = (nconds,Nr), name="fluxes_positive")
        vn = m.addMVar(lb=0,ub=max_flux,shape = (nconds,Nr), name="fluxes_negative")
      
    #dictionary of matrix variables 
    mvars = {'b':b, 'v':v, 'drGp':drGp, 'drGn':drGn, 'ln_conc':ln_conc, 'drG':drG, 'drG_error':drG_error, 'drG_conc':drG_conc}

    if split_v:
        mvars['vp'] = vp
        mvars['vn'] = vn
    
    #import bounds from tmodel 
    for cond_index, cond in enumerate(conds):
        ln_conc[cond_index,:].lb = [np.log(met.lower_bound.m_as('M')) for met in tmodel.metabolites]
        ln_conc[cond_index,:].ub = [np.log(met.upper_bound.m_as('M')) for met in tmodel.metabolites]
        
        v[cond_index,:].lb = [rxn.lower_bound for rxn in tmodel.reactions]
        v[cond_index,:].ub = [rxn.upper_bound for rxn in tmodel.reactions]

        if split_v:
            vp[cond_index,:].ub = [abs(rxn.upper_bound) for rxn in tmodel.reactions]
            vn[cond_index,:].ub = [abs(rxn.lower_bound) for rxn in tmodel.reactions]
        
        for indx, met in enumerate(tmodel.metabolites):
            if met.ignore_conc:
                ln_conc[cond_index,indx].lb = 0
                ln_conc[cond_index,indx].ub = 0
                
        m.update()
    
    
        #add constriants 
        #drG concnetration term
        m.addConstr(drG_conc[cond_index] == (S.T@ln_conc[cond_index])*RT, name = 'drG_conc_constraint_'+str(cond))
        
        #drG term 
        m.addConstr(drG[cond_index] == drGt0tr + drG_error + drG_conc[cond_index], name = 'drG_constraint_'+str(cond)) #+ epsilon 

        #mass balance ## add a feastol variable to allow the MIP start to be feasible when loading fluxes ##doesn't not have a significant impact on the solution
        feastol=m.params.FeasibilityTol
        #it's a variable so it can be minimized in the objective function ..
        mass_balance = 'else'
        feasvar = m.addVar(lb=0, ub=0, name = 'feasvar')

        if mass_balance=='strict':
            m.addConstr(S@v[cond_index] == 0, name = 'mass_balance_'+str(cond))
        else:
            m.addConstr(S@v[cond_index] <= 0+feasvar, name = 'mass_balance_neg'+str(cond))
            m.addConstr(S@v[cond_index] >= 0-feasvar, name = 'mass_balance_pos'+str(cond))


        #for charge and proton transporters the model is infeasible if the membrane potential is 0 and the pH is the same in the two compartments
        #in this case we will not add the drG constraints for these reactions allowing them to be freely reversible (unless otherwise constrained)
        #this issue arrises becuse pH and membrane potentials are fixed in the model and not allowed to vary unlike other metabolite concentrations
        charge_transport_rxns = tmodel.get_charge_transporters()
        proton_transport_rxns = tmodel.get_proton_transporters()
        charge_proton_rxns = charge_transport_rxns + proton_transport_rxns
        
        for i, rxn in enumerate(tmodel.reactions):
            if not rxn.ignore_snd:
                if not all([rxn.drG0prime == 0, rxn.drGtransport==0, rxn in charge_proton_rxns]): #catch infeasible proton and charge transport 
                        #mixed integer constraints for flux direction 1000 is max flux 
                        if split_v is False:
                            m.addConstr(v[cond_index,i] <= max_flux*b[cond_index,i], name='integer_rev_'+rxn.id+'_'+str(cond))
                            m.addConstr(v[cond_index,i] >= -max_flux*(1-b[cond_index,i]), name = 'integer_fwd_'+rxn.id+'_'+str(cond))
                        if split_v is True:
                            m.addConstr(v[cond_index,i] == vp[cond_index,i]-vn[cond_index,i],name='v_net'+rxn.id+'_'+str(cond))
                            m.addConstr(vp[cond_index,i] <= max_flux*b[cond_index,i], name='integer_rev_'+rxn.id+'_'+str(cond))
                            m.addConstr(vn[cond_index,i] <= max_flux*(1-b[cond_index,i]), name = 'integer_fwd_'+rxn.id+'_'+str(cond))

                        #mixed integer constraint for drG to set direction based on + or - 
                        if big_M is False:
                            max_drG = min(max_drG, max(abs(drG[cond_index,i].lb),abs(drG[cond_index,i].ub)))
                        else:
                            max_drG = big_M

                        m.addConstr(drG[cond_index,i] == drGp[cond_index,i]-drGn[cond_index,i],name='drG_net'+rxn.id+'_'+str(cond))
                        m.addConstr(drGn[cond_index,i] <= max_drG*b[cond_index,i],name = 'drG_neg_'+rxn.id+'_'+str(cond))
                        m.addConstr(drGp[cond_index,i] <= max_drG*(1-b[cond_index,i]), name = 'drG_pos_'+rxn.id+'_'+str(cond))
                        m.addConstr((drGn[cond_index,i]+drGp[cond_index,i]) >= epsilon, name = 'drG_epsilon_'+rxn.id+'_'+str(cond))
                         
        m.update()
                
    #condition independent variables and constraints i.e. drG errors
    if error_type == 'covariance':
        #calculate reaction covariance matrix 
        standard_dgr_Q = tmodel._drG0_cov_sqrt

        

        Nq = standard_dgr_Q.shape[1]
        m_bound = scipy.stats.chi2.ppf(alpha, Nq) #chi2 bound for individual degrees of freedom

        qm = m.addMVar(lb=-m_bound,ub=m_bound,shape=Nq,name="qm")

        mnorm = m.addVar(lb=0, ub = m_bound**0.5,name='mnorm')
        
        m.addConstr(drG_error == standard_dgr_Q @ qm, name = 'drG_error_constraint')

        if type(qnorm) is int: #qorm constraint can be ignored - essentially comletely free (but still correlated) errors - can be useful to make regression easier 
            m.addConstr(mnorm == gp.norm(qm,qnorm) , name = 'quad_m') 
            #an alternative formulation of the constraint - not sure which is faster?
            #m.addConstr((qm[known_dof]@qm[known_dof]) <= scipy.stats.chi2.ppf(alpha, Nq), name = 'quad_m')
        
        elif qnorm == 'sep_norm':## as described in equilibrator's doc, we can replace l2 norm on qm with upper and lower bounds for each qm iseparately, corresponding to a confidence interval on each individual degree of freedom in the uncertainty:
            mnorm = m.addMVar(lb=0, ub = m_bound,shape=Nq,name='mnorm')
            m.addConstrs((mnorm[i] == gp.abs_(qm[i]) for i in range(Nq)))
            deg1_bound = scipy.stats.chi2.ppf(alpha, 1) ## chi2 ppf     for 1 degree of freedom
            m.addConstrs((mnorm[i] <= deg1_bound**0.5 for i in range(Nq)))
            
        mvars['qm'] = qm
        mvars['mnorm'] = mnorm
    else:
        #linear independent error
        drG_GAMS_error = m.addMVar(lb=-1.95,ub=1.95,shape=Nr,name="drG_GAMS_error") #should this be 1.96 or 1.98?
        mvars['drG_GAMS_error'] = drG_GAMS_error
        m.addConstrs((drG_error[i] == drG_SE[i]*drG_GAMS_error[i] for i in range(Nr)), name = 'box_error_constraint')
        
        #nullspace constraint must also be added when using linear errors
        if nullspace is not None:   
            drG_null = m.addMVar(lb=-100000,ub=100000,shape=Nr,name="drG_null")
            mvars['drG_null'] = drG_null
            dfG_null = np.array([met.dfG0.m for met in tmodel.metabolites])
            #remove protons and charge for nullspace calculation
            S_null = S.copy()
            for i, met in enumerate(tmodel.metabolites):
                if met.id in ['h_b','h_c','h_m','h_e','charge_c','charge_m','charge_e']:
                    S_null[i] = 0  
            drG0_null = dfG_null@S_null
        
            m.addConstr(drG_null == drG0_null+drG_error, name ="drG_null_cons" ) #use untransformed drG values from multi_dG
            m.addConstrs((0 == drG_null @ nullspace.T[i] for i in range(nullspace.shape[1])),name = 'Nullspace_constraint')


        
    #gibbs enery dissipation constraints
    if gdiss_constraint:

        #gibbs energy dissipation limit
        g_diss_lim = tmodel.T.m_as('K')*sigmac_limit #convert units to kJ/mol
        
        int_index = [idx for idx, rxn in enumerate(tmodel.reactions) if rxn not in tmodel.boundary]
        ex_index = [idx for idx, rxn in enumerate(tmodel.reactions) if rxn in tmodel.boundary]

        if error_type == 'covariance':
            Gdiss = m.addMVar(shape=(nconds,len(ex_index)),lb=-100000, ub = 1000000, name = "Gdiss")
            g_2 = m.addMVar(lb=0, ub = g_diss_lim, shape = nconds, name="g_2")
            mvars['Gdiss'] = Gdiss
            mvars['g_2'] = g_2
            for cond_index, cond in enumerate(conds):
                m.addConstrs((Gdiss[cond_index,ii] == drG[cond_index,i] * v[cond_index,i] for ii, i in enumerate(ex_index)), name = 'Gdiss_constraint_'+str(cond))
                m.addConstr(g_2[cond_index] == gp.quicksum(Gdiss[cond_index]), name = 'Gdiss_sum_constraint_'+str(cond))



        elif error_type == 'linear':
            

            #Gibbs energy balance 
            #drG term omitting concentration - this is condition independant
            drG2 = m.addMVar(lb=-10000, ub = 10000,shape=Nr,name="drG2")
            m.addConstr(drG2 == drGt0tr + drG_error, name = 'drG2_constraint')
            mvars['drG2'] = drG2

            #Gdiss term is condition dependent
            Gdiss_ex = m.addMVar(shape=(nconds,len(ex_index)),lb=-100000, ub = 1000000, name = "Gdiss_ex")
            g_2 = m.addMVar(lb=0, ub = 2*g_diss_lim, shape = (nconds), name="g_2")
            mvars['Gdiss_ex'] = Gdiss_ex
            mvars['g_2'] = g_2

            Gdiss_int = m.addMVar(shape=(nconds,len(int_index)),lb=-100000, ub = 1000000, name = "Gdiss_int")
            g_1 = m.addMVar(lb=-2*g_diss_lim, ub = 0,shape=(nconds), name="g_1")
            mvars['Gdiss_int'] = Gdiss_int
            mvars['g_1'] = g_1

            #total Gdiss needs to include concentration term 
            #can just fix limit on exchange reactions because the internal reactions are linked via the stoichiometry 
            Gdiss_ex_conc = m.addMVar(lb=-g_diss_lim, ub = g_diss_lim, shape=nconds, name="Gdiss_ex_conc")
            mvars['Gdiss_ex_conc'] = Gdiss_ex_conc

            for cond_index, cond in enumerate(conds):
                m.addConstrs((Gdiss_ex[cond_index,ii] == drG2[i] * v[cond_index,i] for ii, i in enumerate(ex_index)), name= 'Gdiss_ex_constraint_'+str(cond))   
                m.addConstr(g_2[cond_index] == gp.quicksum(Gdiss_ex[cond_index]))

                m.addConstrs((Gdiss_int[cond_index,ii] == drG2[i] * v[cond_index,i] for ii, i in enumerate(int_index)), name = 'Gdiss_int_constraint_'+str(cond))
                m.addConstr(g_1[cond_index] == gp.quicksum(Gdiss_int[cond_index]))

                m.addConstr((g_1[cond_index]+g_2[cond_index]) == 0, name = 'GEB_p_'+str(cond)) #ensure gibbs energy balance for internal and exhcange reactions 
                #m.addConstr((g_1[cond_index]+g_2[cond_index]) >= -1e-3, name = 'GEB_n_'+str(cond)) #ensure gibbs energy balance for internal and exhcange reactions 

                #concnetration term can be calulcated for just exchange metabolites
                m.addConstr(((ln_conc[cond_index]@S[:,ex_index])@v[cond_index,ex_index])*RT == Gdiss_ex_conc[cond_index], name='Gdiss_lnc_'+str(cond))
                m.addConstr((g_2[cond_index]+Gdiss_ex_conc[cond_index]) <= g_diss_lim, name= 'Gdiss_lim_ex_'+str(cond))

                #Gdiss limit not needed for internal reactions? should be constrained already by external... 
                #m.addConstr((g_1[cond_index]+Gdiss_ex_conc[cond_index]) >= g_diss_lim, name= 'Gdiss_lim_int_'+str(cond))

    # combined objective to maximise over all conditions if more than one condition is defined
    if tmodel.objective.direction == 'min':
        direction = GRB.MINIMIZE
        feastol_sign = 1
    else:
        direction = GRB.MAXIMIZE 
        feastol_sign = -1
    m.setObjective(gp.quicksum(v[:, [i for i, rxn in enumerate(tmodel.reactions) if rxn.objective_coefficient == 1][0]])+feastol_sign*feasvar, direction)


    m.params.NonConvex = 2
    m.params.TimeLimit =10

    m.update()
    
    return m, mvars



def total_cell_conc(tmodel,conds = [''], metabolites = [], volume_data = None, extracellular = None, GAMS_style = False):
    '''add total cell concentrations constraints'''

    mvars = tmodel.mvars
    m=tmodel.m
    
    #make vector of cell compartment volumes (extracellular compartment may be excluded)

    comps = [i for i in list(tmodel.compartments.keys()) if i in volume_data.columns] ## only consider compartments for which we have data on

    if extracellular is not None:
        comps.remove(extracellular)

    ncomps = len(comps)
    
    #compartment volume array based on conditions (order must match comps list)
    comp_vol_array = volume_data[comps].values

    
    nconds = len(conds)
    #number of meausred metabolites or metabolites where we want to add total cell concentration constraints 
    Nmetmeas = len(metabolites) 
    
    #define total cell concentration variable - these will be regressed to data
    cell_conc = m.addMVar(lb=1e-9,ub=0.8,shape = (nconds,Nmetmeas),name="cell_conc") #in M 
    mvars['cell_conc'] = cell_conc
    
    #vector of metabolite concentrations across compartments
    # make this as narrow as possible to help with exp and log linear aproximation?
    met_conc_dist = m.addMVar(lb=0,ub=1,shape = (nconds,Nmetmeas,ncomps),name="met_conc_dist")

    mvars['met_conc_dist'] = met_conc_dist
    m.update()
    
    #load variables
    ln_conc = mvars['ln_conc']
      
    for cond_index, cond in enumerate(conds):
        #loop through each measurement and add required constraints 
        for met_indx, met_id in enumerate(metabolites):
            comp_indxs = [] #index of all the comparmtents the metabolite is found in 
            for model_indx, met in enumerate(tmodel.metabolites):
                if met.id[:-2] == met_id:
                    if met.compartment in comps:
                        comp_indx = comps.index(met.compartment)
                        comp_indxs.append(comp_indx)

                        #link constrain the concentration of each metaboite to the met_conc_dist vector
                        m.addGenConstrExp(ln_conc[cond_index,model_indx], met_conc_dist[cond_index,met_indx,comp_indx],name="exp_ln_conc_["+met.id+"]"+'_'+str(cond))                     
            #for each measurement add a variable for cellular concentraiton 
            if GAMS_style:
                if comp_indxs == [comps.index('c')]:
                    #print('just_c',met_id, met_conc_dist[cond_index,met_indx,comp_indxs],comp_vol_array[cond_index][comp_indxs])
                    m.addConstr(met_conc_dist[cond_index,met_indx,comp_indxs]@np.array([1]) == cell_conc[cond_index,met_indx], name = 'cell_conc_constraint_['+met_id+']'+'_'+str(cond))
                else:
                    if met_id != 'nad':
                        m.addConstr(met_conc_dist[cond_index,met_indx,comp_indxs]@comp_vol_array[cond_index][comp_indxs] == cell_conc[cond_index,met_indx], name = 'cell_conc_constraint_['+met_id+']'+'_'+str(cond))
                if met_id == 'nad': #nad measuremnt is just cytosolic 
                    comp_indxs = [comps.index('c')]
                    m.addConstr(met_conc_dist[cond_index,met_indx,comp_indxs]@np.array([1]) == cell_conc[cond_index,met_indx], name = 'cell_conc_constraint_['+met_id+']'+'_'+str(cond))


            else:
                m.addConstr(met_conc_dist[cond_index,met_indx,comp_indxs]@comp_vol_array[cond_index][comp_indxs] == cell_conc[cond_index,met_indx], name = 'cell_conc_constraint_['+met_id+']'+'_'+str(cond))

    m.update()

    return m, mvars

def regression(tmodel,m, mvars,conds, flux_data, metabolite_data, volume_data,
               conc_fit=True, flux_fit=True, drG_fit=True, resnorm=1, qm_resnorm = 2,
               error_type = 'covariance', conc_units = None,extracellular=None):
    '''Add variables and constraints for regression to data.
    
    Parameters
    ----------
    tmodel: ThermoModel
        ThermoModel object to be used for regression
    m: gurobipy.Model
        gurobipy model object to add variables and constraints to
    mvars: dict
        dictionary of variables to be updated with new variables
    conds: list
        list of conditions to be used for regression
    flux_data: pandas.DataFrame
        dataframe of flux data to be regressed to
    metabolite_data: pandas.DataFrame
        dataframe of metabolite data to be regressed to
    volume_data: pandas.DataFrame
        dataframe of volume data to be used for total cell concentration constraints
    conc_fit: bool
        if True add constraints for metabolite concentration regression
    flux_fit: bool
        if True add constraints for flux regression
    drG_fit: bool
        if True add constraints for drG regression
    resnorm: int
        1 for linear residual formulation (sum of absolute differences), 2 for quadratic residual formulation (sum of squared differences) 
    error_type: str
        'covariance' for covariance based error, 'linear' for linear error
    conc_units: str
        units of concentration data to be regressed to

    Returns
    -------
    m: gurobipy.Model
        gurobipy model object with added variables and constraints
    mvars: dict
        dictionary of variables now updated with new variables
        '''

    
    nconds = len(conds) #number of conditions (for regression to multiple conditions)

    if flux_fit:
        measured_fluxes = list(flux_data.unstack()['mean'].T.index) #list of measured fluxes for indexing
        Nmeas = len(measured_fluxes) #number of fluxes measured
        resflx = m.addMVar(lb=0, ub = GRB.INFINITY,shape=(nconds,Nmeas),name="resflx") #flux residual
        mvars['resflx'] = resflx

    if conc_fit:
        measured_mets = list(metabolite_data.unstack()['mean'].T.index) #list of measured mets for indexing
        Nmetmeas = len(measured_mets) #number of metabolites measured
        resconc = m.addMVar(lb=0, ub = GRB.INFINITY, shape=(nconds,Nmetmeas),name="resconc") #concentration residual
        mvars['resconc'] = resconc

    
    #load variables
    ln_conc = mvars['ln_conc']
    b = mvars['b']
    v = mvars['v']
    
    #set empty onjective
    obj_conc = 0
    obj_flux = 0
    obj_drG = 0
    
    reaction_id = [rxn.id for rxn in tmodel.reactions]

    if conc_fit:
        total_cell_conc(tmodel,conds=conds, metabolites=measured_mets,
                        volume_data=volume_data, extracellular=extracellular)

        cell_conc = mvars['cell_conc']
    
    for cond_index, cond in enumerate(conds):
       
        if conc_fit:

            if conc_units is None:
                raise Exception('No concentration units defined, please specifiy units of metabolite concentration data to be regressed to')

            # define metabolite data to fit
            dmetsmeas = metabolite_data.loc[cond].dropna()

            for i, row in enumerate(dmetsmeas.iterrows()):
                met_id = row[1].name
                met_idx = measured_mets.index(met_id)
                meas_conc = Q_(row[1]['mean'], conc_units).to('M').m  #convert data from user defined units to M for fitting 
                sd = Q_(row[1]['sd'], conc_units).to('M').m
         
                if resnorm == 2:
                    m.addConstr(resconc[cond_index,met_idx] >= ((cell_conc[cond_index,met_idx]-meas_conc)/sd)* ((cell_conc[cond_index,met_idx]-meas_conc)/sd), name = 'resconc_cond_'+met_id+'_'+str(cond))
                else:
                    #linear fomulation of resconc
                    m.addConstr(resconc[cond_index,met_idx] >= ((cell_conc[cond_index,met_idx]-meas_conc)/sd), name = 'resconc_p_'+met_id+'_'+str(cond))
                    m.addConstr(resconc[cond_index,met_idx] >= (-(cell_conc[cond_index,met_idx]-meas_conc)/sd), name = 'resconc_n_'+met_id+'_'+str(cond))
           
            #concnetration fit objective term 
           # obj_conc = (resconc.sum()/(Nmetmeas*nconds))
            obj_conc = (resconc.sum())

        if flux_fit:
            dvmeas = flux_data.loc[cond] 
            sd_array=[]
            for i, row in enumerate(dvmeas.iterrows()):
                rxn_id = row[1].name
                flux = row[1]['mean']
                sd = row[1]['sd']
                model_idx = reaction_id.index(rxn_id)
                flx_idx = measured_fluxes.index(rxn_id)
                
                #set hint that solution is near a good fit
                v[cond_index, model_idx].VarHintVal = flux
                
                if abs(flux) >= 0: #ignore unmeasured fluxes dont assume they are 0
                  
                    #fix at least the direction of measured fluxes to match those measured 
                    if flux > 0:
                        b[cond_index,model_idx].lb = 1
                        b[cond_index,model_idx].ub = 1
                        v[cond_index,model_idx].lb = 1e-5 #set non-zero flux bound to force flux 
                    if flux < 0:
                        b[cond_index,model_idx].ub = 0
                        b[cond_index,model_idx].lb = 0
                        v[cond_index,model_idx].ub = -1e-5 #set non-zero flux bound to force flux 


                    if resnorm == 2:
                        m.addConstr(resflx[cond_index,flx_idx] >= ((v[cond_index,model_idx]-flux)/sd)*((v[cond_index,model_idx]-flux)/sd), name = 'resflx_cond_'+rxn_id+'_'+str(cond_index))
                    else:
                        #linear formulation of resflx 
                        m.addConstr(resflx[cond_index,flx_idx] >= ((v[cond_index,model_idx]-flux)/sd), name = 'resflx_p_'+rxn_id+'_'+str(cond_index))
                        m.addConstr(resflx[cond_index,flx_idx] >= (-(v[cond_index,model_idx]-flux)/sd), name = 'resflx_n_'+rxn_id+'_'+str(cond_index))
            
            #flux fit objective term 
            #obj_flux = (resflx.sum()/(Nmeas*nconds)) 
            obj_flux = (resflx.sum())
                        
    if drG_fit:
        drG_SE = [rxn.drG_SE.m for rxn in tmodel.reactions]
        drGt0tr = np.array([(rxn.drG0prime.m + rxn.drGtransport.m) for rxn in tmodel.reactions])
        unknown_drG = [i for i, SE in enumerate(drG_SE) if SE >= 100]
        known_drG = [i for i, SE in enumerate(drG_SE) if SE < 100]   
        
        if error_type != 'covariance':
            drG_GAMS_error = mvars['drG_GAMS_error']

            norm_drG_GAMS_error = m.addMVar(lb=0, ub = GRB.INFINITY, shape=len(tmodel.reactions),name="drG_error_GAMS_norm") 
            m.addConstrs((norm_drG_GAMS_error[i] >= (drGt0tr[i]+(drG_GAMS_error[i]*drG_SE[i])) for i in unknown_drG), name = 'drG_error_GAMS_norm_P')
           # m.addConstrs((norm_drG_GAMS_error[i] >= (-1*(drGt0tr[i]+(drG_GAMS_error[i]*drG_SE[i]))) for i in unknown_drG), name = 'drG_error_GAMS_norm_N')


            adrGerr = m.addMVar(lb=0, ub = GRB.INFINITY,shape=(len(tmodel.reactions)),name="adrGerr") 
            m.addConstrs((adrGerr[i] == drG_GAMS_error[i]**2 for i in known_drG), name = 'adrGerr')
            #m.addConstrs((adrGerr[i] == (drGt0tr[i]+(drG_GAMS_error[i]*drG_SE[i])) for i in unknown_drG), name = 'adrGerr')

            known_res = m.addVar(name='known_res')
            unknown_res = m.addVar(name='unknown_res')
            m.addConstr(known_res == gp.quicksum(adrGerr[known_drG]))
            #m.addConstr(unknown_res == gp.quicksum(adrGerr[unknown_drG]))
            m.addConstr(unknown_res == gp.quicksum(norm_drG_GAMS_error[unknown_drG]))

            alpha = 0.05
           # obj_drG = ((known_res/len(known_drG))+(alpha*(unknown_res)/len(unknown_drG))) #0.05 is alpha value from Niebel et.al 2019
            obj_drG = (known_res + alpha*unknown_res)
            print(alpha)
            mvars['adrGerr'] = adrGerr
            mvars['known_res'] = known_res
            mvars['unknown_res'] = unknown_res
            mvars['norm_drG_GAMS_error'] = norm_drG_GAMS_error
            
        else:
            if qm_resnorm == 2:
                obj_drG = mvars['qm']@mvars['qm']
            elif qm_resnorm == 1:
                obj_drG = mvars['mnorm']
           

    m.setObjective(obj_conc + obj_flux + obj_drG, GRB.MINIMIZE)

    return m, mvars

def regression_legacy(tmodel,m, mvars,conds, flux_data, metabolite_data, volume_data,
               conc_fit=True, flux_fit=True, drG_fit=True, resnorm=1, 
               error_type = 'linear'):

    '''legacy function to set up regession optimization to replicated GAMS results
    included incorect assigning of cytosolic volume to entire cell'''
    
    nconds = len(conds) #number of conditions (for regression to multiple conditions)

    measured_fluxes = list(flux_data.unstack()['mean'].T.index) #list of measured fluxes for indexing
    Nmeas = len(measured_fluxes) #number of fluxes measured
    
    measured_mets = list(metabolite_data.unstack()['mean'].T.index) #list of measured mets for indexing
    Nmetmeas = len(measured_mets) #number of metabolites measured

    #variable to define the regression residual 
    resflx = m.addMVar(lb=0, ub = GRB.INFINITY,shape=(nconds,Nmeas),name="resflx") #flux residual
    resconc = m.addMVar(lb=0, ub = GRB.INFINITY, shape=(nconds,Nmetmeas),name="resconc") #concentration residual
    mvars['resflx'] = resflx
    mvars['resconc'] = resconc
    
    #define total cell concentration variable - these will be regressed to data
    cell_conc = m.addMVar(lb=1e-6,ub=600,shape = (nconds,Nmetmeas),name="total_cell_conc") #in mM 
    mvars['cell_conc'] = cell_conc

    #these bounds can affect the accuracy of the general constraint for log to exp - set as narrow as possible 
    max_conc = 1 #M 
    min_conc = 1e-8 #M 
    conc = m.addMVar(lb=min_conc,ub=max_conc,shape = (nconds,len(tmodel.metabolites)),name="met_conc") #in M
    mvars['conc'] = conc
    m.update()
    
    #load variables
    ln_conc = mvars['ln_conc']
    b = mvars['b']
    v = mvars['v']
    
    #set empty onjective
    obj_conc = 0
    obj_flux = 0
    obj_drG = 0
    
    reaction_id = [rxn.id for rxn in tmodel.reactions]
    
    for cond_index, cond in enumerate(conds):
        #define metabolite data to fit
        dmetsmeas = metabolite_data.loc[cond].dropna()


        if conc_fit:

            #loop through each measurement and add required constraints 
            for i, row in enumerate(dmetsmeas.iterrows()):
                met_id = row[1].name
                met_idx = measured_mets.index(met_id)
                meas_conc = row[1]['mean']
                sd = row[1]['sd']


                #add conversion from log to normal space for measured metabolite
                c_conc = 0
                m_conc = 0
                for ii, met in enumerate(tmodel.metabolites):
                    if met.id[:-2] == met_id:
                        if met.compartment != 'e':
                            m.addGenConstrExp(ln_conc[cond_index,ii],conc[cond_index,ii],name="exp_ln_conc_["+met.id+"]"+'_'+str(cond))
                        if met.compartment == 'c':
                            c_conc = conc[cond_index,ii]
                            c_lnconc = ln_conc[cond_index,ii]

                        if met.compartment == 'm':
                            m_conc = conc[cond_index,ii]
                            m_lnconc = ln_conc[cond_index,ii]

                #calculate total cell concentration
                #to match GAMS (which I think is wrong) if metabolite is cytosolic assume measurment is 100% cytosolic
                if type(m_conc) == gp.MVar:
                    c_vol = volume_data.loc[cond,'c']
                    m_vol = volume_data.loc[cond,'m']

                    m.addConstr(cell_conc[cond_index,met_idx] == (1000*c_conc*c_vol + 1000*m_conc*m_vol), name = 'cell_conc_'+met_id+'_'+str(cond)) #cell conc is in mM to match data
                else:
                    m.addConstr(cell_conc[cond_index,met_idx] == 1000*c_conc, name = 'cell_conc_'+met_id+'_'+str(cond))

                if resnorm == 2:
                    m.addConstr(resconc[cond_index,met_idx] >= ((cell_conc[cond_index,met_idx]-meas_conc)/sd)* ((cell_conc[cond_index,met_idx]-meas_conc)/sd), name = 'resconc_cond_'+met_id+'_'+str(cond))
                else:
                    #linear fomulation of resconc
                    m.addConstr(resconc[cond_index,met_idx] >= ((cell_conc[cond_index,met_idx]-meas_conc)/sd), name = 'resconc_p_'+met_id+'_'+str(cond))
                    m.addConstr(resconc[cond_index,met_idx] >= (-(cell_conc[cond_index,met_idx]-meas_conc)/sd), name = 'resconc_n_'+met_id+'_'+str(cond))


                  
           
            #concnetration fit objective term 
            obj_conc = (resconc.sum()/(Nmetmeas*nconds))

        if flux_fit:
            dvmeas = flux_data.loc[cond] 
            sd_array=[]
            for i, row in enumerate(dvmeas.iterrows()):
                rxn_id = row[1].name
                flux = row[1]['mean']
                sd = row[1]['sd']
                model_idx = reaction_id.index(rxn_id)
                flx_idx = measured_fluxes.index(rxn_id)
                
                #set hint that solution is near a good fit
                v[cond_index, model_idx].VarHintVal = flux
                
                if abs(flux) >= 0: #ignore unmeasured fluxes dont assume they are 0
                  
                    #fix at least the direction of measured fluxes to match those measured 
                    if flux > 0:
                        b[cond_index,model_idx].lb = 1
                        b[cond_index,model_idx].ub = 1
                        v[cond_index,model_idx].lb = 1e-4 #set non-zero flux bound to force flux 
                    if flux < 0:
                        b[cond_index,model_idx].ub = 0
                        b[cond_index,model_idx].lb = 0
                        v[cond_index,model_idx].ub = -1e-4 #set non-zero flux bound to force flux 


                    if resnorm == 2:
                        m.addConstr(resflx[cond_index,flx_idx] >= ((v[cond_index,model_idx]-flux)/sd)*((v[cond_index,model_idx]-flux)/sd), name = 'resflx_cond_'+rxn_id+'_'+str(cond_index))
                    else:
                        #linear formulation of resflx 
                        m.addConstr(resflx[cond_index,flx_idx] >= ((v[cond_index,model_idx]-flux)/sd), name = 'resflx_p_'+rxn_id+'_'+str(cond_index))
                        m.addConstr(resflx[cond_index,flx_idx] >= (-(v[cond_index,model_idx]-flux)/sd), name = 'resflx_n_'+rxn_id+'_'+str(cond_index))
            
            #flux fit objective term 
            obj_flux = (resflx.sum()/(Nmeas*nconds)) 
                        
    if drG_fit:
        drG_SE = [rxn.drG_SE.m for rxn in tmodel.reactions]
        unknown_drG = [i for i, SE in enumerate(drG_SE) if SE >= 1000]
        known_drG = [i for i, SE in enumerate(drG_SE) if SE < 1000]   
        

        if error_type != 'covariance':
            drG_GAMS_error = mvars['drG_GAMS_error']

            adrGerr = m.addMVar(lb=0, ub = GRB.INFINITY,shape=(len(tmodel.reactions)),name="adrGerr") 
            m.addConstrs((adrGerr[i] == drG_GAMS_error[i]**2 for i in known_drG), name = 'adrGerr')
            m.addConstrs((adrGerr[i] == (drG_GAMS_error[i]**2)*drG_SE[i] for i in unknown_drG), name = 'adrGerr')

            known_res = m.addVar(name='known_res')
            unknown_res = m.addVar(name='unknown_res')
            m.addConstr(known_res == gp.quicksum(adrGerr[known_drG]))
            m.addConstr(unknown_res == gp.quicksum(adrGerr[unknown_drG]))
        
            obj_drG = ((known_res/len(known_drG))+(0.05*(unknown_res)/len(unknown_drG))) #0.05 is alpha value from Niebel et.al 2019
            
        else:
            obj_drG = mvars['qm']@mvars['qm']
    

    m.setObjective(obj_conc + obj_flux + obj_drG, GRB.MINIMIZE)

    return m, mvars


def variability_analysis(tmodel, vars = []):
    '''Set up a multiscenario optimisation problem to perform variability analysis on a given variable

    Parameters
    ----------
    tmodel: ThermoModel
        ThermoModel object to be used for variability analysis
    vars: list
        list of variables to perform variability analysis on

    Returns
    -------
    m: gurobipy.Model
        gurobipy model object with added variables and constraints
    '''

    m = tmodel.m
    m.NumScenarios = 0 #reset any existing scenarios
    #reset any current objective 
    for variable in m.getVars():
        variable.Obj = 0 
    

    #in case var is a gurobi mvar then loop through and convert to a list of individual variables
    vars = [var for var in vars]

    m.NumScenarios = (len(vars))*2
    print(m.NumScenarios)
    m.update()
    
    #clear any previous scenario objectives 
    for i in range(m.NumScenarios):
        m.params.ScenarioNumber=i
        for variable in m.getVars():
            variable.ScenNObj = 0
        
    min_i = 0 
    max_i = 1
        
    for i, var in enumerate(vars):
        
        m.params.ScenarioNumber=min_i
        var.ScenNObj = -1
        min_i += 2
        m.params.ScenarioNumber=max_i
        var.ScenNObj = 1
        max_i += 2
    
    m.update()

    return m


def variability_results(tmodel):

    m = tmodel.m

    no_scenarios = m.NumScenarios
    if no_scenarios > 0:
        obj_val = {}
        obj_bound = {}
        optimal_bounds = {}
        MIPGaps = {}
        for i in range(0, no_scenarios, 2):
            var_idx = int(i / 2)

            # Minimization:
            m.params.ScenarioNumber = i
            m.update()
            var_name = [v.varName for v in m.getVars() if v.ScenNObj != 0][0]
            ObjBound = m.ScenNObjBound
            ObjVal = m.ScenNObjVal
            if ObjVal != 0:
                MIPGap = abs((ObjBound-ObjVal)/ObjVal)
            else:
                MIPGap = 0

            obj_val[var_name] = [(-1) * ObjVal]
            obj_bound[var_name] = [(-1) * ObjBound]
            MIPGaps[var_name] = [MIPGap]

            if MIPGap <= m.params.MIPGap:
                optimal_bounds[var_name] = [(-1) * ObjBound]
            else:
                optimal_bounds[var_name] = [float('nan')]

            # Maximization:
            m.params.ScenarioNumber = i + 1
            m.update()
            var_name = [v.varName for v in m.getVars() if v.ScenNObj != 0][0]

            ObjBound = m.ScenNObjBound
            ObjVal = m.ScenNObjVal
            if ObjVal != 0:
                MIPGap = abs((ObjBound-ObjVal)/ObjVal)
            else:
                MIPGap = 0
            obj_val[var_name].append((+1) * m.ScenNObjVal)
            obj_bound[var_name].append((1) * m.ScenNObjBound)
            MIPGaps[var_name].append(MIPGap)
            if MIPGap <= 0.0001:
                optimal_bounds[var_name].append((1) * ObjBound)
            else:
                optimal_bounds[var_name].append(float('nan'))

        return optimal_bounds, obj_val, obj_bound, MIPGaps

def read_bounds(bounds_file):
    '''Read in bounds text file from hpc optimisation and returns a dataframe of the bounds.  
    Parameters
    ----------
        bounds_file: str 
            Path to bounds file with format variable_name: [lb, ub]
        
    Returns
    -------
        drG_bounds: pandas.DataFrame 
            DataFrame of bounds for each drG
    '''
    drG_bounds = {}
    with open(bounds_file, 'r') as f:  
        for line in f:
            line.strip()
            var_name = line.split(':')[0]
            #index = int((line.split()[0][0:-1]).split(',')[1][:-1])
            lb = float(line.split()[1][1:-1])
            
            ub = float(line.split()[2][0:-1])
            drG_bounds[var_name] = {'lb': lb, 'ub': ub}

    return pd.DataFrame(drG_bounds).T



    return pd.DataFrame(drG_bounds).T
  
def compute_IIS(tmodel):

    m=tmodel.m
    m.computeIIS() #compute irreducible inconsistent subset 

    IIS_constr_list = []
    IIS_bound_list = []


    #loop through model constraints
    for c in m.getConstrs():
        if c.IISConstr ==1: # if the constraint has IIS = 1 flag then print it
            print(c, c.RHS) 
            IIS_constr_list.append(c.ConstrName)  

    for c in m.getQConstrs():
        if c.IISQConstr ==1:
            print(c, c.QCRHS) 
            IIS_constr_list.append(c)

    for c in m.getGenConstrs():
        if c.IISGenConstr ==1:
            print(c) 
            IIS_constr_list.append(c)

    for bnd in m.getVars(): #loop through the variables 
        if bnd.IISLB ==1:
            print(bnd, bnd.LB,'LB') #if the lower bound is in the IIS print it
            IIS_bound_list.append(bnd)
        if bnd.IISUB ==1:
            print(bnd, bnd.UB, 'UB') #same for upper bound 
            IIS_bound_list.append(bnd)


    infeas_drg_rxns = []
    infeas_drg_mets = []
    for varr in IIS_bound_list:
        var = varr.varName
        print(var)
        if any(var_id in var for var_id in ['drG_']):
            #use re to get number form []
            if ',' in var:
                rxn = int(var.split(',')[1].split(']')[0])
            else:
                rxn = int(var.split('[')[1].split(']')[0])
                
            infeas_drg_rxns.append(tmodel.reactions[rxn])

        if any(var_id in var for var_id in ['ln_conc']):
            if ',' in var:
                met = int(var.split(',')[1].split(']')[0])
            else:
                met = int(var.split('[')[1].split(']')[0])
            infeas_drg_mets.append(tmodel.metabolites[met])
            
    return IIS_constr_list, IIS_bound_list, list(set(infeas_drg_rxns)), list(set(infeas_drg_mets))


def model_start(tmodel, sol_file, ignore_vars = [],fix_vars=[], fix='start' ):
    """Set the model start point from a saved gurobi solution
    fix_vars: the variables that you want to import the starting points from
    fix: ['start', hint, bounds], 
    Note that for some solutions numerical issue make an imported start point ifeasible. Therefore it 
    is beneficial to only fix some start points and allow the remining unknwon variables to be calcualted"""
    
    #never try and fix the solution of calcualted variables as this can cause numerical issues
    #these values are calcualted form other fixed paremters so they do not need to be fixed
    ignore_vars.extend(['met_conc','resflx','cell_conc','resconc','drG_conc','met_conc_dist'])
    
    m = tmodel.m
    
    with open(sol_file, 'r') as f:  
        for line in f:
            line.strip()
            if len(line.split()) ==2:
                (names, xn) = line.split()
                
                if 'all' not in ignore_vars:
                
                    if not any(a in names for a in ignore_vars):
                        xn = float(xn)
                        var = m.getVarByName(names)
                        if var is not None:
                            if 'start' in fix:
                                var.Start = xn
                            if 'hint' in fix:
                                var.VarHintVal = xn
                            if 'bound' in fix:
                                var.LB = xn
                                var.UB = xn
                            
                if  any(a in names for a in fix_vars):
                    xn = float(xn)
                    var = m.getVarByName(names)
                    if var is not None:
                        if 'start' in fix:
                            var.Start = xn
                        if 'hint' in fix:
                            var.VarHintVal = xn
                        if 'bound' in fix:
                            var.LB = xn
                            var.UB = xn


def gdiss_var(tmodel, var, verbose=False):
    '''Return the gibbs energy dissiaption of a specific variable'''

    ex_index = [idx for idx, rxn in enumerate(tmodel.reactions) if rxn in tmodel.boundary]
    int_index = [idx for idx, rxn in enumerate(tmodel.reactions) if rxn not in tmodel.boundary]

    # print(var)
    ex_error = (tmodel.mvars['v'].x * var)[0,ex_index].sum()
    int_error = (tmodel.mvars['v'].x * var)[0,int_index].sum()

    if verbose:
        print('Ex:',ex_error)
        print('Int:',int_error)
        print('Diff:',int_error+ex_error)

    return ex_error, int_error
        
                 

def gdiss_model(tmodel):
    '''Print the individual gibbs energy balance for different components of drG. This can be used to identify 
    why the Gibbs enery balance might not be met'''
    

    print("drG error term")
    var = (tmodel.mvars['drG_error'].x)
    gdiss_var(tmodel,var, verbose=True )


    print("\n")
    print("RTlnC concentration term")
    var = (tmodel.mvars['drG_conc'].x)
    gdiss_var(tmodel,var, verbose=True )

    print("\n")
    print("drG0' term")
    var = [rxn.drG0prime for rxn in tmodel.reactions]
    var = np.array([val.m if val is not None else 0 for val in var])
    gdiss_var(tmodel,var, verbose=True )

    print("\n")
    print("drG0 term non-transformed")
    var = ([rxn.drG0 for rxn in tmodel.reactions])
    var = np.array([val.m if val is not None else 0 for val in var])
    gdiss_var(tmodel,var, verbose=True)

    print("\n")
    print("drG0' transform term")
    var = ([rxn.drGtransform for rxn in tmodel.reactions])
    var = np.array([val.m if val is not None else 0 for val in var])
    gdiss_var(tmodel,var, verbose=True)

    print("\n")
    print("drG0' total transport term")
    var = ([rxn.drGtransport for rxn in tmodel.reactions])
    var = np.array([val.m if val is not None else 0 for val in var])
    gdiss_var(tmodel,var, verbose=True)

    print("\n")
    print("drG0' charge transport term")
    var = ([rxn.drG_c_transport for rxn in tmodel.reactions])
    var = np.array([val.m if val is not None else 0 for val in var])
    gdiss_var(tmodel,var, verbose=True)

    print("\n")
    print("drG0'proton transport term")
    var = ([rxn.drG_h_transport for rxn in tmodel.reactions])
    var = np.array([val.m if val is not None else 0 for val in var])
    gdiss_var(tmodel,var, verbose=True)


def get_solution(tmodel):

    # create a pandas dataframe to store the solution
    #define number of conditions and loop through them
    sols = []
    for cond_index in range(tmodel.mvars['v'].shape[0]):
        
        sol = pd.DataFrame(columns = ['reaction', 'v', 'drG','Gdiss','drG0_prime', 'drG_error', 'drG_conc', 'b'])
        sol['reaction'] = [rxn.id for rxn in tmodel.reactions]
        sol['condition'] = cond_index
        sol['v'] = tmodel.mvars['v'][cond_index].x
        sol['drG'] = tmodel.mvars['drG'][cond_index].x
        sol['Gdiss'] = sol['v']*sol['drG']
        sol['drG0_prime'] = [rxn.drG0prime.m+rxn.drGtransport.m for rxn in tmodel.reactions]
        sol['drG_error'] = tmodel.mvars['drG_error'].x
        sol['drG_conc'] = tmodel.mvars['drG_conc'][cond_index].x
        sol['b'] = tmodel.mvars['b'][cond_index].x

        sol.set_index('reaction', inplace=True)
        sols.append(sol)

    return pd.concat(sols)


def variable_scan(tmodel, scan_range, var):
    '''Scan a variable over a range of values and return the solution for each value'''

    if tmodel.m is None:
        print('No tmodel has been created, cannot perform variable scan')

        return None

    else:       
        tmodel.m.NumScenarios = 0 #reset any existing scenarios
        tmodel.m.NumScenarios = len(scan_range)
        
        
        for i, val in enumerate(scan_range):
            tmodel.m.params.ScenarioNumber=i
            tmodel.m.update()

            var.ScenNUB = val
            var.ScenNLB = val
            tmodel.m.update()
                    
        tmodel.m.update()


def multi_scenario_sol(tmodel, var):
    '''Return the solution for a variable from a multiple scenario optimization'''

    cond_vars = []

    m = tmodel.m
    mvars = tmodel.mvars
    
    #number of conditions
    n_conds = mvars[var].shape[0]

    # loop through the conditions
    for cond_index in range(n_conds):
        var_vals = []
        # loop through the scenarios
        for scen in range(m.NumScenarios):
            m.params.ScenarioNumber = scen
            m.update()
            try:
                var_vals.append(mvars[var][cond_index].ScenNX)
            except:
                var_vals.append(np.zeros((mvars[var][cond_index]).shape))
                
        
        flux_sol =np.array(var_vals)
        cond_vars.append(flux_sol)
        
    return np.array(cond_vars)
        

def drG_bounds(tmodel, concentration = True,drG_error = True, alpha=0.95, condition_index = 0):
    '''calculate the drg bounds for a given model
    
    modified from equilibrator-api examples https://equilibrator.readthedocs.io/en/latest/equilibrator_examples.html '''

    if tmodel.m is None:
        print('No tmodel has been created, cannot calculate drG bounds')
        
    import cvxpy
    import scipy
    from ..utils.vis import progressbar

    S = create_stoichiometric_matrix(tmodel) 

    standard_dgr_Q = tmodel._drG0_cov_sqrt


    Nc, Nr = S.shape
    Nq = standard_dgr_Q.shape[1]
    lb = tmodel.mvars['ln_conc'][condition_index].lb
    ub = tmodel.mvars['ln_conc'][condition_index].ub

    ln_conc = cvxpy.Variable(shape=Nc, name="metabolite log concentration")
    m = cvxpy.Variable(shape=Nq, name="covariance degrees of freedom")

    constraints = [
        lb <= ln_conc,  # lower bound on concentrations
        ln_conc <= ub,  # upper bound on concentrations
        cvxpy.norm2(m) <= scipy.stats.chi2.ppf(alpha, Nq) ** (0.5)  # quadratic bound on m based on confidence interval
    ]

    drG_mean = np.array([rxn.drG0prime.m+rxn.drGtransport.m for rxn in tmodel.reactions])

    if (concentration) & (drG_error):
        dg_prime = (drG_mean + 
            (R*tmodel.T).m_as("kJ/mol") * S.T @ ln_conc +
            standard_dgr_Q @ m
        )
    elif drG_error:
        dg_prime = (drG_mean + 
            standard_dgr_Q @ m
        )
    elif concentration:
        dg_prime = ((R*tmodel.T).m_as("kJ/mol") * S.T @ ln_conc)

    drG_bounds = {}

    for i, rxn in enumerate(progressbar(tmodel.reactions,'', 40, item_label_attribute = 'id')):
        if (rxn.drG_SE.m < tmodel._rmse_inf.m): #ignore unbounded reaction drG
            try:
                prob_max = cvxpy.Problem(cvxpy.Maximize(dg_prime[i]), constraints)
                prob_max.solve(verbose=False)
                max_df = prob_max.value

                prob_min = cvxpy.Problem(cvxpy.Minimize(dg_prime[i]), constraints)
                prob_min.solve(verbose=False)
                min_df = prob_min.value

                drG_bounds[rxn.id] = [min_df, max_df]

                
            except:
                print(i, rxn.id, 'cvxpy error')

                try:
                    prob_max = cvxpy.Problem(cvxpy.Maximize(dg_prime[i]), constraints)
                    prob_max.solve(verbose=True)
                    max_df = prob_max.value

                    prob_min = cvxpy.Problem(cvxpy.Minimize(dg_prime[i]), constraints)
                    prob_min.solve(verbose=True)
                    min_df = prob_min.value

                    drG_bounds[rxn.id] = [min_df, max_df]

                    print(i, rxn.id, f"{min_df:.1f}, {max_df:.1f} kJ/mol")
                except:
                    print(i, rxn.id, 'could not be solved')
                    pass

    bounds_df = pd.DataFrame.from_dict(drG_bounds, orient='index', columns=['lb','ub'])

    return bounds_df

def calc_conc_bounds(tmodel, conds,  metabolite_data = None, extracellular_data = None, volume_data= None, extracellular = 'e', CI = 2.6, conc_units = None):
    '''Calcualte concentration bounds for individual metabolites based on data from whole cell metabolomics and extracellular measurements.

    Parameters
    ----------
    tmodel: ThermoModel
        Thermodynamic model
    conds: list
        List of conditions to calculate bounds for (must match data)
    metabolite_data: pd.DataFrame
        Dataframe of metabolite concentrations (mean and sd) for each condition
    extracellular_data: pd.DataFrame
        Dataframe of extracellular concentrations (lo and up) for each condition (to match legacy GAMS data)
    volume_data: pd.DataFrame
        Dataframe of volume fractions for each compartment for each condition
    extracellular: str
        Name of extracellular compartment (default = 'e')
    CI: float
        Zvalue for Confidence interval for data (default = 2.6 to match GAMS legacy ~99.9% CI)
    conc_units: str
        Units of concentration data (default = mM)

    Returns
    -------
    data_bounds_conc_df: pd.DataFrame
        Dataframe of bounds for each metabolite for each condition
    '''
    if conc_units is None:
        raise Exception('No concentration units defined, please specifiy units of metabolite concentration data')

    data_bounds_conc = {}

    for cond in conds:
        data_bounds_conc[cond] = {met.id:{'lb':0, 'ub':0} for met in tmodel.metabolites}  

    comps = list(tmodel.compartments.keys())

    for cond in conds:

        metabolites = list(metabolite_data.loc[cond].index)
        for met_id in metabolites:
            comps = [] #list of all the comparmtents the metabolite is found in 
            compartmented_measurement = False #flag to indicate if the measurement was for a specific compartment
            for met in tmodel.metabolites:
                 #if metabolite data has _compartment then assume its a compartment specific measurement
                if met.id == met_id:
                    if met.compartment != extracellular:
                        comps.append(met.compartment)
                        compartmented_measurement = True                      

                elif "_".join((met.id).split('_')[:-1]) == met_id:  #string split to avoid ambiguity when using met.id[:-2] e.g. nadph[:-2] == nad
                    if met.compartment != extracellular:
                        comps.append(met.compartment)

            meas_conc = Q_(metabolite_data.loc[cond,met_id]['mean'], conc_units).to('M').m  #convert data from user defined units to M for fitting 
            sd = Q_(metabolite_data.loc[cond,met_id]['sd'], conc_units).to('M').m

            #if a metabolite is only found in one compartment then the bounds are set to the data +/- CI*sd
            if len(comps) == 1:
                #if the measurement was for a specific compartment then the bounds are set to the data +/- CI*sd
                if compartmented_measurement:
                    data_bounds_conc[cond][met_id]['lb'] = meas_conc - CI*sd
                    data_bounds_conc[cond][met_id]['ub'] = meas_conc + CI*sd
                else: #whole cell measuremnt but metabolite is only found in one compartment in the model

                    data_bounds_conc[cond][met_id+'_'+comps[0]]['lb'] = meas_conc - CI*sd
                    data_bounds_conc[cond][met_id+'_'+comps[0]]['ub'] = meas_conc + CI*sd

            #if a metabolite is found in multiple compartments then the bounds are set to the data +/- CI*sd and then converted to the compartment of interest using the volume fraction data
            #only ub can be set assuming all of a measured metabolite is in a single compartment, lb is 0 or default lb
            else:
                for comp in comps:
                    if comp != extracellular:
                        data_bounds_conc[cond][met_id+'_'+comp]['ub'] = (1/volume_data.loc[cond,comp]*(meas_conc + CI*sd))       

        if extracellular_data is not None:
            ex_mets = list(set(extracellular_data.loc[cond]['met'].values))
            for met_id in ex_mets:
                lb = Q_(extracellular_data.loc[(extracellular_data['met']==met_id) & (extracellular_data['bound']=='lo') ].loc[cond]['Value'], conc_units).to('M').m
                ub = Q_(extracellular_data.loc[(extracellular_data['met']==met_id) & (extracellular_data['bound']=='up') ].loc[cond]['Value'], conc_units).to('M').m
                data_bounds_conc[cond][met_id+'_'+extracellular]['lb'] = lb
                data_bounds_conc[cond][met_id+'_'+extracellular]['ub'] = ub

    #melt data bounds conc into dataframe with condition, metabolite, lo and up columns
    data_bounds_conc_df = pd.DataFrame.from_dict({(i,j): data_bounds_conc[i][j] 
                            for i in data_bounds_conc.keys() 
                            for j in data_bounds_conc[i].keys()},
                        orient='index')

    #drop rows of data_bounds_conc_df that have 0 for both lo and up
    data_bounds_conc_df = data_bounds_conc_df.loc[(data_bounds_conc_df!=0).any(axis=1)]

    return data_bounds_conc_df

