"""
Evaluate reaction directions based on flux variability on carbon minimal medium
"""
import numpy as np
import pandas as pd

import cobra

from growth_coupling_suite.reaction_directionality \
    import config_reaction_directionality_default as config_default

from growth_coupling_suite.util.util import check_config

from os.path import dirname, abspath
import os

from concurrent import futures
import multiprocessing as mp


"""
- integrate heterologous reactions in the host model
- conduct FVA on every possible carbon source (minimal medium) at near optimal growth
    - 20 mmol/gDW/h substate uptake rate (Feist et al. 2007)
    - 20 mmol/gDW/h O2 uptake rate
- evaluate flux directions
    - if a reaction can carry flux, is there a fix direction for the tested media?
    - is flux restricted to positive or negative rates under any tested condition?

"""


# %% get hme directory of this module
module_dir = dirname(abspath(config_default.__file__))


# %%


def get_reaction_direction(model_in, rxn_list=[], config=config_default):
    """
    evaluate flux variability of reactions and suggest reaction direction
    only consider minimal media with a single carbon source

    Parameters
    ----------
    model_in : cobra.core.model
        metabolic model in COBRA format.
    rxn_list : list of strings or cobra.core.reaction, optional
        list of reactions (string or cobra.core.reaction) considered for 
        flux variability analysis. The default is [].
    config : module, optional
        module containing optional parameters and settings. The default is config_default.

    Returns
    -------
    reaction_direction_conditions : pandas DataFrame
        reaction directions suggested by flux variability analysis on different substrates
        (1): irreversible in forward direction
        (-1): irreversible in backward direction
        (0): reversible
        (-2): blocked reaction
        (-3): evaluation failed
    flux_variability_conditions : list
        list of results

    """

    model = model_in.copy()
    # check configuration file
    check_config(config, config_default)
    
    # check and transform reaction list
    if not(rxn_list):
        rxn_list_cobra = list(model.reactions)
    else:
        rxn_list_cobra = []
        for rxn in rxn_list:
            try:
                # add reaction in cobra format
                rxn_list_cobra.extend(model.reactions.get_by_any(rxn))
            except:
                pass
        
     
    # load conditions to be tested
    if not(config.filename_carbon_exchange_reactions):
        filename_ex_rxns = "{0}/carbon_exchange_reaction_list.xlsx".format(module_dir)
    else:
        filename_ex_rxns = config.filename_carbon_exchange_reactions
     # check file
    if not(os.path.isfile(filename_ex_rxns)):
        raise TypeError("No exchange reactions provided for carbon substrates!")
    
    carbon_uptake_reactions = load_carbon_upake_reactions(filename_ex_rxns)
    # check existence of uptake reactions
    for ex_rxns in carbon_uptake_reactions:
        try:
            model.reactions.get_by_id(ex_rxns)
        except:
            print(ex_rxns + " not in model")
            carbon_uptake_reactions.remove(ex_rxns)
    
    if len(carbon_uptake_reactions)<=0:
        raise TypeError("No exchange reactions provided for carbon substrates!")
    else:
        print("Test flux variabilities for " + str(len(carbon_uptake_reactions))
              + " substrates ...")
    
    # deactivate uptake of any specified carbon substrate 
    for ex_rxns in carbon_uptake_reactions:
        model.reactions.get_by_id(ex_rxns).lower_bound = 0
    
    # conduct flux variability analysis for specified conditions
    
    # prepare model
    model.reactions.get_by_id("EX_o2_e").lower_bound = -20
           

    # set up and execute FVA in multiple processes
    args = ((model, rxn_list_cobra, config.fraction_of_optimum, b) \
        for b in carbon_uptake_reactions)
     
    if config.processes >= mp.cpu_count():
        processes = mp.cpu_count()-1
    else:
        processes = config.processes
    if processes > len(carbon_uptake_reactions):
        processes = len(carbon_uptake_reactions)
        
    with futures.ProcessPoolExecutor(max_workers=processes) as e:
        res_fva_all = e.map(flux_variability_analysis_condition, args)
    
    # evaluate FVA results
    reaction_direction_conditions = []
    flux_variability_conditions = []
    for res_fva in res_fva_all:
        if len(res_fva)<3:
            print("Flux variability analysis failed for " + res_fva[1])
            # no valid solution for carbon uptake reaction
            continue
        
        flux_variability_conditions.append(res_fva)
        reaction_direction_conditions.append(evaluate_flux_variability_result(res_fva[0],
                                                                              res_fva[1],
                                                                              res_fva[2]))
    reaction_direction_conditions = pd.concat(reaction_direction_conditions,
                                              axis=1, sort=False)  
    
    
    return reaction_direction_conditions, flux_variability_conditions
        
def evaluate_flux_variability_result(res_fva, res_fba, carbon_uptake_reaction=[]):
    # evaluate flux variability solution for expected flux direction (DataFrame object)
    
    flux_threshold = 1e-4
    
    if not(carbon_uptake_reaction):
        col_name = "direction"
    else:
        col_name = carbon_uptake_reaction
    
    reaction_direction = pd.DataFrame({col_name: [n for n in range(len(res_fva))]},
                                      index= [r.name for r in res_fva.iloc])
    for rxn in res_fva.iloc:
        # compare minimum to maximum value
        if (rxn.minimum < -flux_threshold) and (rxn.maximum > flux_threshold):
            # flux can be positive and negative -> reversible reaction
            reaction_direction.loc[rxn.name, col_name] = 0
        elif rxn.minimum >= -flux_threshold and rxn.maximum > flux_threshold:
            # flux can only be positive -> irreversible in forward direction
            reaction_direction.loc[rxn.name, col_name] = 1
        elif rxn.maximum <= flux_threshold and rxn.minimum < -flux_threshold:
            # flux can only be negative -> irreversible in backward direction
            reaction_direction.loc[rxn.name, col_name] = -1
        elif rxn.minimum >= -flux_threshold and  rxn.maximum <= flux_threshold:
            # blocked reaction
            reaction_direction.loc[rxn.name, col_name] = -2
        else:
            reaction_direction.loc[rxn.name, col_name] = -3
            print("WARNING: Unknown flux direction! Min: " + str(rxn.minimum)
                  + " Max: " + str(rxn.maximum))
            
    return reaction_direction
        
       
        
        

def flux_variability_analysis_condition(args):
    # conduct flux variability analysis for a single condition
    
    
    
    # get parameter
    model = args[0]
    rxn_list_cobra = args[1]
    fraction_of_optimum = args[2]
    carbon_uptake_reaction = args[3]
    
    
    # set uptake rates
    model.reactions.get_by_id(carbon_uptake_reaction).lower_bound = -20
    
    # test model
    sol = model.optimize()
    if sol.status != "optimal":
        print("Model is infeasible on " +  carbon_uptake_reaction + ". Try to rescue ...")
        for rxn in rxn_list_cobra:
            with model:
                rxn_model = model.reactions.get_by_id(rxn.id)
                if not(rxn_model.lower_bound <= 0.1 and rxn_model.upper_bound >= 0.1):
                    # release reaction
                    rxn_model.lower_bound = -1000
                    rxn_model.upper_bound = 1000
                    # test model
                    sol_release = model.slim_optimize()
                else:
                    continue
                    
            if sol_release>1e-2:
                rxn_model = model.reactions.get_by_id(rxn.id)
                rxn_model.lower_bound = -1000
                rxn_model.upper_bound = 1000
                # calculate complete solution
                sol = model.optimize()
                break
            
    if sol.status != "optimal":
        print("Model is still infeasible!")
        return [-1, carbon_uptake_reaction]
        
    print("Conduct FVA for substrate: " + carbon_uptake_reaction + 
          " (Optimal objective value " +  str(round(sol.objective_value, 2)) + ")")
    
    sol_fba = pd.DataFrame()
    if True:
        # manual FVA
        # constrain objective function
        for rxn in model.reactions:
            if rxn.objective_coefficient==1:
                objective_rxn = rxn
                break
        objective_rxn.lower_bound = sol.objective_value*fraction_of_optimum
        objective_rxn.upper_bound = sol.objective_value
        
        sol_fva = pd.DataFrame()
        for rxn in rxn_list_cobra:
            # save physiological flux rate
            sol_fba.loc[rxn.id, carbon_uptake_reaction] = sol.get_primal_by_id(rxn.id)
            
            with model:
                # release bounds on reaction
                model.reactions.get_by_id(rxn.id).lower_bound = -1000
                model.reactions.get_by_id(rxn.id).upper_bound = 1000
                # FVA
                sol_fva_rxn = cobra.flux_analysis.variability.flux_variability_analysis(model, 
                                        reaction_list=[rxn],
                                        fraction_of_optimum=fraction_of_optimum,
                                        processes=1)
            
                sol_fva = pd.concat([sol_fva, sol_fva_rxn])
        
    else:
        # FVA
        sol_fva = cobra.flux_analysis.variability.flux_variability_analysis(model, 
                                        reaction_list=rxn_list_cobra,
                                        fraction_of_optimum=fraction_of_optimum,
                                        processes=1)

    return [sol_fva, sol_fba, carbon_uptake_reaction]


    
def load_carbon_upake_reactions(filename):
    # get list of exchange reactions for carbon source uptake
    if ".xlsx" in filename:
        # load Excel file as data frame
        df = pd.read_excel(filename, sheet_name="exchange_reactions")
    elif ".csv" in filename:
        # load .csv file as data frame
        df = pd.read_csv(filename)
    else:
        raise TypeError("Unknown file format. Use .xlsx or .csv format")

    # save reaction IDs of substrate exchange reactions
    carbon_uptake_reactions = list(df.loc[:, "reaction_id"])
    
    return carbon_uptake_reactions

