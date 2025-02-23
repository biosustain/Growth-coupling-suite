"""Provide methods for the deterination of significant design variables of strain designs


"""

from growth_coupling_suite.util import util

from warnings import warn
from typing import TYPE_CHECKING, Optional, Dict, List
import pandas as pd
import numpy as np
from itertools import combinations

if TYPE_CHECKING:
    from cobra import Model, Reaction
    

def design_significance(
        model: 'Model',
        design_dict: Dict[str, dict],
        reverse_design_dict: Dict[str, dict],
        target_reaction: str,
        objective: Optional[str] = None,
        objective_direction: Optional[str] = None,
        additional_constraints: Dict[str, str] = {},
        eval_gpr: bool = False,
        minimum_relative_significance: float = 0.1
        ) -> (pd.DataFrame, dict):
    
    
    """test which design objects (deletions, addins etc) are relevant for enforcing design objective
      
    Parameters
    ----------
    model : cobra.Model
        Metabolic model in COBRA format
    design_dict : dict
        All design variables (knockouts, carbon source etc) with additional information
    reverse_design_dict : dict
        Specifies bounds of design variables for the wildtype   
    target_reaction : str, optional
        target reaction which flux is used to test significance. The default is None.
    objective : str, optional
        objective function reaction. The default is None.
    objective_direction : str, optional
        direction for objective function optimization. The default is None.
    additional_constraints : dict, optional
        apply additional constraints to the model. The default is {}.
    minimum_relative_significance : float, optional
        minimum allowable fraction of objective function value for reduced design

    Returns
    -------
    combination_results : DataFrame
        Significance of all possible combinations of intervention objects in the design.
    design_small : dict
        Smallest strain design that has the same effect on the target reaction flux as the full design.

    """
    
    
    # parameter
    # minimum_relative_significance = 0.1
    minimum_absolute_significance = 1e-4
    
    # save each design in a dict and append to a list
    combination_results_dicts = [] 
    combination_results_insig_dicts = []
    
    with model:
        # set objective
        if objective: 
            model.objective = objective
        if objective_direction:
            model.objective_direction = objective_direction
            
        # apply additional constraints
        # dictionary with tuples of lower and upper bound
        for r_id, cnstr in additional_constraints.items():
            if r_id in model.reactions:
                r = model.reactions.get_by_id(r_id)
                # set bounds
                r.lower_bound = cnstr[0]
                r.upper_bound = cnstr[1]
            else:
                warn("Reaction " + r_id + " not found in model", UserWarning)
         
        # simulate objective with full design
        err_val = "infeasible"
        if objective == target_reaction:
            target_flux_full = model.slim_optimize(error_value=err_val)
        else:
            sol_full = model.optimize()
            if sol_full.status == "infeasible":
                target_flux_full = err_val
            else:
                target_flux_full = sol_full.fluxes[target_reaction]
                
            
        # check objective function value of full design
        if (target_flux_full == err_val) or (target_flux_full == 0):
            warn("Full strain design is infeasible! Try to rescue design", UserWarning)
            # print("Full strain design is infeasible! Try to rescue design")
            full_design_feasible = False
            target_flux_full = 0
            
        else:
            full_design_feasible = True
                   
        # prepare testing of design object combinations
        num_i = len(design_dict)
        obj_keys = list(design_dict.keys())
        idx_i = [i for i in range(num_i)]
        
        # specify keys for source/cofeed variables and all other targets
        obj_keys_medium = [
            o_key
            for o_key, o in design_dict.items()
            if o['type']=='source' or o['type']=='cofeed'
            ]    
        num_i_medium = len(obj_keys_medium)
        idx_i_medium = [obj_keys.index(o_key) for o_key in obj_keys_medium]
        
        obj_keys_genes = [
            o_key
            for o_key, o in design_dict.items()
            if o['type']!='source' and o['type']!='cofeed'
            ]    
        num_i_genes = len(obj_keys_genes)
        idx_i_genes = [obj_keys.index(o_key) for o_key in obj_keys_genes]
        
        
        # create combinations of medium compositions (if there are more than 2)
        obj_comb_medium = [()]
        for num in range(num_i_medium-1):
            obj_comb_medium.extend(combinations(idx_i_medium, num+1))  
    
        # evaluate GPR and exclude infeasible combinations   
        # no reaction targets are added here due to a GPR. Only couples of reactions from the design are identified
        if eval_gpr:
            ID2Idx = {design_dict[obj_keys[i]]["ID"]: i for i in idx_i} # link reaction IDs to their indices
      
            gpr_couples = [] # obligatory reaction deletion couples
            for i in idx_i:
                
                # skip source or cofeed variables -> not gene related
                if design_dict[obj_keys[i]]["type"] in ['source', 'cofeed']:
                    # gpr_couples.append([i])

                    continue
      
                # determine genes of reaction target
                # genes = model.reactions.get_by_id(design_dict[obj_keys[i]]["ID"]).genes
                # determine reactions affected by target reaction deletion
                rxns_linked = util.GPR_linked_reaction_knockouts(
                                model,
                                [design_dict[obj_keys[i]]["ID"]],
                                eval_gpr=True
                                )
                rxns_linked.append(design_dict[obj_keys[i]]["ID"]) # add base reaction to couple
                
                
          
                # create gpr couple
                gpr_couple = [ID2Idx[r] for r in rxns_linked if r in ID2Idx]
                # check if couple is unique
                in_list = False
                for couple in gpr_couples.copy():
            
                    if set(couple) == set(gpr_couple):
                        in_list = True
                        break
                    
                    elif (len(set(couple) - set(gpr_couple)) == 0):
                        gpr_couples.append(gpr_couple)
                        del gpr_couples[gpr_couples.index(couple)]
                        in_list = True
                        break
                         
                    elif (len(set(gpr_couple) - set(couple)) == 0):
                        in_list = True
                        break
                    
                  
                
                if not(in_list):
                    gpr_couples.append(gpr_couple)
                    
               
    
                # gpr_couples.append([ID2Idx[r] for r in rxns_linked if r in ID2Idx])
  
            # create combinations of gpr couples
            couple_obj_comb = []
            for num in range(len(gpr_couples)-1):
                couple_obj_comb.extend(combinations([i for i in range(len(gpr_couples))], num+1))
                
            # translate gpr couples to reaction indices
            obj_comb_genes = [()]
            for couple_obj in couple_obj_comb:
                obj = []
                for obj_idx in couple_obj:
                    obj.extend(gpr_couples[obj_idx])
                obj_comb_genes.append(tuple(np.unique(obj)))
            
                
    
                        
        else:

            obj_comb_genes = [()]
            for num in range(num_i_genes-1):
                obj_comb_genes.extend(combinations(idx_i_genes, num+1))  

        
        # check all combinations of design objects
        # loop through combinations of medium compositions
        # - if a cofeed is the only carbon source, the strain design must include one or more insertions (growth enabler)
        for comb_medium in obj_comb_medium:
            
            # obj_keys_medium idx_i_medium
            
            
            # is single active cofeed the only carbon source?
            obj_id_medium_inactive = [design_dict[obj_keys[i]]['ID'] for i in comb_medium]
            medium = model.medium
            carbon_sources_active = [rxn_id for rxn_id in medium.keys()
                              for ex_met in model.reactions.get_by_id(rxn_id).metabolites.keys()
                              if 'C' in ex_met.elements and rxn_id not in obj_id_medium_inactive]

            cofeed_is_carbon_source = False
            if len(carbon_sources_active) == 1 and carbon_sources_active[0].startswith('cofeed_'):
                # only cofeed is the carbon source
                cofeed_is_carbon_source = True
                
 
            obj_comb_insign = [] 
            for comb_genes in obj_comb_genes:
                
                # construct combination of variables to be excluded from design
                comb = comb_genes + comb_medium
                
                # if cofeed is the only carbon source an insertion needs to be in the design
                obj_type_active = [design_dict[obj_keys[i]]['type'] for i in idx_i_genes if i not in comb_genes]

                if cofeed_is_carbon_source and not(np.any(np.array(obj_type_active)=='addin')):
                    # dismiss design solution
                    combination_results_insig_dicts.append({
                        "interventions": [design_dict[obj_keys[j]]["ID"] for j in [i for i in idx_i if i not in comb]], # design objects still active
                        "excluded": [design_dict[obj_keys[j]]["ID"] for j in comb],
                        "excluded_index": [j for j in comb],
                        "objective_value": -3,
                        "difference_objective_value": -3,
                        "significance": -3
                        })
                    
                    continue
     
            
                # if combination of excluded design variables contains an 
                # insignificant subset jump to next combination
                is_insignificant = False
                for o in obj_comb_insign:
                    if not(set(o) - set(comb)):
                        is_insignificant = True
                        break
                    
                if is_insignificant:
                    # dismiss design solution
                    combination_results_insig_dicts.append({
                        "interventions": [design_dict[obj_keys[j]]["ID"] for j in [i for i in idx_i if i not in comb]], # design objects still active
                        "excluded": [design_dict[obj_keys[j]]["ID"] for j in comb],
                        "excluded_index": [j for j in comb],
                        "objective_value": -2,
                        "difference_objective_value": -2,
                        "significance": -2
                        })

                    continue 

                # simulate model with reduced design combinations
                obj_keys_restore = [obj_keys[idx] for idx in comb]
                target_flux, target_flux_diff, significance = simulate_reduced_design(
                    model, reverse_design_dict, obj_keys_restore, design_dict,
                    target_reaction, target_flux_full
                    )
                                
                combination_results_dicts.append({
                    "interventions": [design_dict[obj_keys[j]]["ID"] for j in [i for i in idx_i if i not in comb]], # design objects still active
                    "excluded": [design_dict[obj_keys[j]]["ID"] for j in comb],
                    "excluded_index": [j for j in comb],
                    "objective_value": target_flux,
                    "difference_objective_value": target_flux_diff,
                    "significance": significance
                    })

                if significance < minimum_relative_significance and significance >= 0:
                    # design is feasible but insignificant
                    # infeasible subsets of designs may still lead to a smaller feasible subset when knockouts are removed
                    obj_comb_insign.append(comb)
       
                
    # significant designs
    combination_results = pd.DataFrame(combination_results_dicts)
                            
    # add insignificant designs     
    combination_results_insign = pd.DataFrame(combination_results_insig_dicts)

    # merge designs
    combination_results = pd.concat(
        [combination_results, combination_results_insign],
        ignore_index=True
        )
             
                
    # sort results from low to high difference to target flux
    combination_results = combination_results.sort_values(
                                by=["significance"],
                                axis=0,
                                ascending=False)
    
    # get smallest significant design
    if full_design_feasible:
        comb_sign = combination_results.loc[
            # (combination_results.loc[:, "significance"]>0.99), :] #& (combination_results.loc[:, "significance"]<1.01), :]
            (combination_results.loc[:, "significance"]>minimum_relative_significance), :] #& (combination_results.loc[:, "significance"]<1.01), :]

    else:
        comb_sign = combination_results.loc[
            (combination_results.loc[:, "significance"]>minimum_absolute_significance), :]
        
    if len(comb_sign) == 0:
        # no significant subset of the design
        return combination_results, None
    
    design_len = []
    design_i = []
    for i, d in comb_sign.iterrows():
        design_len.append(len(d.loc["interventions"]))
        design_i.append(i)
    comb_small = comb_sign.loc[design_i[design_len.index(min(design_len))], "interventions"]
    
    design_small = design_dict.copy()
    for key, obj in design_dict.items():
        if obj["ID"] not in comb_small:
            # delete design object
            del(design_small[key])
        
                
    return combination_results, design_small




def simulate_reduced_design(
        model: 'Model',
        reverse_design_dict: Dict[str, dict],
        # idx_obj_restore: List[int],
        obj_keys_restore: List[str],
        design_dict: Dict[str, dict],
        target_reaction: str,
        target_flux_full: pd.Series,
        # idx_obj_all: List[int],
        # combination_results: pd.DataFrame
        ) -> (float, float, int):
    
    """ 
    Restore design objects and compute the objective function to evaluate their
    significance to the design objective
    
    :param cobra.core.model model: M-model with applied full design
    :param list idx_obj_restore: Design objects to be resoterd
    :param dict design_dict: design solutions
    :param str target_reaction: evaluation reaction, does not need to be
        included in the objective function
    :param pd.Series target_flux_full: Simulated flux distribution with 
        applied full strain design
    :param pd.DataFrame combination_results: Table with all simulated combinations
        of design objects
    :param list idx_obj_all: indices of all objects from the full design
    
    :return pd.DataFrame combination_results: Table with all simulated combinations
        of design objects
    """
    
    with model as model_comb:
        # revert design objects
        for obj in obj_keys_restore:
            
            obj_i = design_dict[obj] # get design object
            r = model_comb.reactions.get_by_id(obj_i["ID"])
            
            if obj_i["type"] in ["deletion", "mediareduction"]:
                # set original bounds
                r.lower_bound = reverse_design_dict[obj]["lower_bound"]
                r.upper_bound = reverse_design_dict[obj]["upper_bound"]
            elif obj_i["type"] in ["source", "cofeed", "addin"]:
                # deactivate target object
                r.lower_bound = 0
                r.upper_bound = 0
                
        # calculate objective 
        if model.objective == target_reaction:
            target_flux = model_comb.slim_optimize(error_value="infeasible")
            feasible = target_flux != "infeasible"
        else:                            
            sol_comb = model_comb.optimize()
            target_flux = sol_comb.fluxes[target_reaction]
            feasible = sol_comb.status == "optimal"
        
        
        if feasible:
            # calculate significance
            target_flux_diff = target_flux_full - target_flux
            if target_flux_full == 0:
                warn ("Objective of full strain design is zero! Not accounted for in determination of significant designs", UserWarning)
                significance = target_flux
            else:  
                
                significance = 1-(target_flux_diff/target_flux_full)
        else:
            target_flux = None
            target_flux_diff = None
            significance = -1

    return (target_flux, target_flux_diff, significance)   
    # return combination_results










