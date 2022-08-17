"""Methods for reducing model size and target reaction space
"""

import pandas as pd


def blocked_reactions(model, reaction_list=[], remove_blocked_reactions=False):
    """
    identify and remove blocked reactions via flux variability anaylsis, i.e. 
    reactions which cannot carry any flux.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.
    reaction_list : list, optional
        IDs of reactions to be assessed. The default is [].
    remove_blocked_reactions : bool, optional
        If true, remove blocked reactions from the model. The default is False.

    Returns
    -------
    model : cobra.Model
        Metbaolic model without blocked reactions.
    blocked_reactions : list
        IDs of blocked reactions.

    """
   

   # populate list of reactions to be checked
    if len(reaction_list) > 0:
        rxn_list = reaction_list
    else:
        rxn_list = model.reactions
        
    # manual fva (stumbled over problems with copied model)
    fva_res = []
    rxn_id = []
    for rxn in rxn_list:
        with model:
            rxn_id.append(rxn.id)
            model.objective = rxn.id
            # maximize
            model.objective_direction = "max"
            sol_max = model.slim_optimize()
            # minimize
            model.objective_direction = "min"
            sol_min = model.slim_optimize()
            # concatenate
            fva_res.append([sol_min, sol_max])
            
    # create frame
    sol_fva = pd.DataFrame(fva_res, index=rxn_id, columns=["minimum", "maximum"])
        
    # analyze results
    blocked_reactions = []
    for rxn in sol_fva.index:
        if sol_fva.loc[rxn, "minimum"]==0 and sol_fva.loc[rxn, "maximum"]==0:  
            blocked_reactions.append(rxn)
            
    # remove blocked reactions
    if remove_blocked_reactions:
        model.remove_reactions(blocked_reactions, remove_orphans=False)
    
    
    
    return model, blocked_reactions


def essential_reactions(model, reaction_list=[], objective_cutoff=1e-03):
    """
    Identify essential reactions, i.e. reactions which flux cannot be zero

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.
    reaction_list : list, optional
        IDs of reactions to be assessed. The default is [].
    objective_cutoff : float, optional
        Maximum value of an invalid objective. The default is 1e-03.

    Returns
    -------
    essential_reactions : list
        Reaction IDs of all essential reactions.

    """
    
    
    # populate reaction list to be analyzed
    if not reaction_list:
        reaction_list = [rxn for rxn in model.reactions]
    else:
        reaction_list_valid = []
        for rxn in reaction_list:
            try:
                reaction_list_valid.append(model.reactions.get_by_any(rxn)[0])
            except:
                continue
        reaction_list = reaction_list_valid
    
    # identify essential reactions
    essential_reactions = []
    for rxn in reaction_list:
        with model:
            rxn = model.reactions.get_by_id(rxn.id)
            # knockout reaction
            rxn.bounds = (0, 0)
            # print(model.reactions.get_by_id(rxn_id).upper_bound)
            # optimize model
            if model.slim_optimize()  < objective_cutoff:
                # essential reaction
                essential_reactions.append(rxn.id)
                
    return essential_reactions
        
    