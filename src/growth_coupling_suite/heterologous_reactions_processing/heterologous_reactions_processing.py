"""
Methods for deriving heterologous reactions for a model from a database of BIGG models

Population of extracted heterologous reactions for a model are stored in
"heterologous_reactions_database"
"""


from  growth_coupling_suite.heterologous_reactions_processing \
    import model_database as md
    
from  growth_coupling_suite.heterologous_reactions_processing \
    import heterologous_reactions_database as hrd  
    
from growth_coupling_suite.heterologous_reactions_processing \
    import config_heterologous_reactions_default as config_default
      
# from growth_coupling_suite.reaction_directionality \
#     import reaction_fluxvariability_evaluation as rfva
    
    
from growth_coupling_suite.util.util import check_config

from growth_coupling_suite.external_APIs import biggAPI
    
import cobra
import pandas as pd
import numpy as np
    
from os.path import abspath, dirname
import ntpath
from pathlib import Path
import pickle

# %%

hrd_dir = dirname(abspath(hrd.__file__))   
md_dir =  Path(dirname(abspath(md.__file__)))

# load BIGG database API
bigg = biggAPI.BiggAPI()


# %%


def get_heterologous_reactions(model, config=config_default, reprocess=False,
                               save_models=True, model_path=hrd_dir,
                               heterologous_models_id=[]):
    """
    Extract heterologous reactions from various COBRA models and assess their directions

    Parameters
    ----------
    model : cobra.core.model.Model
        original model.
    config : module, optional
        contains optional parameters. The default is config_default.
    reprocess : boolean, optional
        if True create a new database model and overwriting old ones. The default is False.
    save_models : boolean, optional
        Pickle and save heterologous reaction database models 
    model_path : string
        Directory for saving the heterologous reaction database models
    heterologous_models_id : list
        IDs of heterologous models from which reactions are drawn. The default is []
        
    Returns
    -------
    hrd_model_origin : cobra.core.model.Model
        COBRA model containing all identified heterologous reactions. Reaction directions were not assessed
    hrd_model : COBRA model
        COBRA model containing  all identified heterologous reactions. Reaction directions assessed
        

    """
    
    # check config file
    check_config(config, config_default)
    
    # load heterologous model ids
    h_model_ids = []
    for h_model in md_dir.glob("*.json"): # glob(md_dir + "/*.json"):
        # h_model_id = str(h_model).split("\\")[-1].replace(".json", "")
        h_model_id = ntpath.split(h_model)[1].replace(".json", "")
        if h_model_id == model.id:
            continue
        elif len(heterologous_models_id)>0:
            if h_model_id in heterologous_models_id:
                h_model_ids.append(h_model_id)
        else:
            h_model_ids.append(h_model_id)
            
    
    
    hrd_model_name = model.id + "_hr_database"
   
    
    if reprocess:
        hrd_model_origin = []
    else:
        hrd_model_origin = load_heterologous_database_model(
            model_name=hrd_model_name,
            model_path=model_path,
            reaction_direction_assessed=False
            )
    
    
    
    if reprocess or not(hrd_model_origin):
        print("Create a new heterologous reaction database")
        # set up new database model in a cobra format
        hrd_model_origin = cobra.Model(hrd_model_name)
        #setattr(hrd_model_origin, "heterologous_model_sources", [])
     
    # load heterologous models   
    h_models = []    
    for h_model_id in h_model_ids:
        print("\t" + h_model_id)
        h_models.append(cobra.io.load_json_model(str(md_dir.joinpath(h_model_id + ".json"))))

            
    # load heterologous models if new database is constructed
    print("Update heterologous reaction database...")
    # h_models = []
    # for h_model_id in h_model_ids:
    #     if not(h_model_id in hrd_model_origin.heterologous_model_sources):
    #         print("\t" + h_model_id)
    #         h_models.append(cobra.io.load_json_model(str(md_dir.joinpath(h_model_id + ".json"))))
    #     else:
    #         print("\t" + "Reactions of " + h_model_id + " already in heterologous reaction database")

    # if not(hrd_model_origin):
    #     print(h_model_ids)
    #     for h_model_id in h_model_ids:
    #         print("\t" + h_model_id)
    #         h_models.append(cobra.io.load_json_model(str(md_dir.joinpath(h_model_id + ".json"))))
   
    # get list of reactions in database
    rxns_in_database = []
    if len(hrd_model_origin.reactions) > 0:
        for rxn in hrd_model_origin.reactions:
            rxns_in_database.append(rxn.id)
     
    # extract heterologous reactions        
    rxns_to_add = {}
    rxns_to_add_id = []    
    for h_model in h_models:
        print("Extract heterologous reactions from: " + h_model.id)
        rxns_to_add_from_model \
            = extract_heterologous_reactions_from_model(model, h_model, 
                                                        currency_mets=config.currency_mets,
                                                        num_carbon_threshold=config.num_carbons_threshold,
                                                        exclude_list=[])
        
        print("\t" + str(len(rxns_to_add_from_model)) + " heterologous reactions found")
        
        # get IDs of heterolgous reactions
        # rxns_to_add_id.extend(list(rxns_to_add_from_model.keys()))
        # add ID of model origin
        for rxn_id in rxns_to_add_from_model:
            if rxn_id in rxns_in_database:
                # already extracted, add origin
                #hrd_model_origin.reactions.get_by_id(rxn_id).origin.append(h_model.id)
                pass
            elif rxn_id in rxns_to_add_id:
                # already extracted, add origin
                # print("\tDuplicate heterologous reaction: "+ rxn_id)
                #rxns_to_add[rxn_id].origin.append(h_model.id)
                pass
            else:
                # new heterologous reaction
                rxn = h_model.reactions.get_by_id(rxn_id)
                # print("\tNew heterologous reaction: " + rxn_id)
                setattr(rxn, "origin", [h_model.id])
                # save reaction
                rxns_to_add[rxn_id] = rxn
                rxns_to_add_id.append(rxn_id)
    
    # exclude reactions with irreasonable effect on growth
    sol_wt = model.slim_optimize()
    rxns_to_add_iter = rxns_to_add.copy()
    for rxn_id, rxn in rxns_to_add_iter.items():
        model.add_reactions([rxn])
        sol_h = model.slim_optimize()
        if sol_h > (sol_wt*2):
            del(rxns_to_add[rxn_id])
            del(rxns_to_add_id[rxns_to_add_id.index(rxn_id)])
            print("Disregard " + rxn_id + " (Objective value: " + str(round(sol_h, 3)) + ")")
            
        # remove reaction again and continue
        model.remove_reactions([model.reactions.get_by_id(rxn_id)])
    
    
    # save heterologous reaction database, unassessed
    extend_database_model(hrd_model_origin, rxns_to_add)
    
    # save unassessed heterologous database model
    if save_models:
        save_heterologous_database_model(
            hrd_model_origin,
            hrd_model_name,
            model_path,
            reaction_direction_assessed=False
            )
 
        
    if config.directionality_assessment:
        # continue with reaction direction assessment         
        hrd_model = model_directionality_assessement(
            hrd_model_origin,
            config=config,
            model_name=hrd_model_name
            )
        
        # save heterologous database model
        if save_models:
            save_heterologous_database_model(
                hrd_model,
                hrd_model_name,
                model_path,
                reaction_direction_assessed=True
                )
    
    else:
        hrd_model = None
        
        


    return hrd_model, hrd_model_origin
    
   
    
   
    # # set up heterologous reaction database model with assessed directions
        
    # # check if heterologous reactions database exist for model, directions assessed
    # # hrd_model_feature = "_dir_assessed"
    
    # if reprocess:
    #     hrd_model = []
    # else:
    #     hrd_model = load_heterologous_database_model(
    #         model_name=hrd_model_name,
    #         model_path=model_path,
    #         reaction_direction_assessed=True
    #         )
          
    # if reprocess or not(hrd_model):
    #     # set up new database model in a cobra format
    #     hrd_model = cobra.Model(hrd_model_name)
    #     #setattr(hrd_model, "heterologous_model_sources", [])
        
    # # get original reactions from database
    # rxns_to_add = {}
    # rxns_to_add_id = []
    # for rxn in hrd_model_origin.reactions:
    #     rxn_id = rxn.id
    #     rxns_to_add[rxn_id] = rxn
    #     rxns_to_add_id.append(rxn_id)
    
    
    # # check which reactions are not assessed yet
    # for rxn in hrd_model.reactions:
    #     if rxn.id in rxns_to_add_id:
    #         # reaction already assessed for direction, delete
    #         del(rxns_to_add[rxn.id])
    #         del(rxns_to_add_id[rxns_to_add_id.index(rxn.id)])
    
    # # evaluate reactions in terms of reaction directionalities
    # if len(rxns_to_add) > 0:
    #     rxns_to_add, directions = evaluate_reaction_directionalities(model, rxns_to_add,
    #                                                       config)
      
    #     # extend database models with heterologous reactions
    #     extend_database_model(hrd_model, rxns_to_add)
         
    #     # save reaction direction (annotated and from thermodynamic analysis)
    #     if hasattr(hrd_model, "reaction_directions"):
    #         hrd_model["reaction_directions"] = \
    #             pd.concat([hrd_model["reaction_directions"], directions])
    #     else:
    #         setattr(hrd_model, "reaction_directions", directions)
    
    # # save database
    # if save_models:
    #     save_heterologous_database_model(
    #         hrd_model,
    #         hrd_model_name,
    #         model_path,
    #         reaction_direction_assessed=True
    #         )
    #     # model_dir = Path(model_path).joinpath(
    #     #     hrd_model_format.format(hrd_model_name, hrd_model_feature))
    #     # with open(str(model_dir), "wb") as f:
    #     #     pickle.dump(hrd_model, f)

                            
    # print("\tNumber of heterologous reactions: " + str(len(hrd_model.reactions)))
        
    # return hrd_model, hrd_model_origin


def model_directionality_assessement(model, config=config_default, model_name="hr_database"):
    # write the docstrings for the reaction
    """
    assess reaction directions of heterologous reactions in a model

    Parameters
    ----------
    model : cobra.Model
        COBRA model.
    config : module, optional
        contains optional parameters. The default is config_default.
    model_name : str, optional
        name of the model. The default is "hr_database".
    
    Returns
    -------
    hrd_model : cobra.Model
        COBRA model containing all identified heterologous reactions. Reaction directions assessed

    """
 
    
    hrd_model = cobra.Model(model_name)
    
    # get unassessed reactions from database
    rxns_to_add = {}
    rxns_to_add_id = []
    for rxn in model.reactions:
        rxn_id = rxn.id
        rxns_to_add[rxn_id] = rxn
        rxns_to_add_id.append(rxn_id)
           
    if len(rxns_to_add) > 0:
        rxns_to_add, directions = evaluate_reaction_directionalities(model, rxns_to_add,
                                                          config)
      
        # extend database models with heterologous reactions
        extend_database_model(hrd_model, rxns_to_add)
        
    return hrd_model


def extend_database_model(hrd_model, rxns_to_add):
    """
    exend heterologous reaction database model with new reactions

    Parameters
    ----------
    hrd_model : cobra.core.model.Model
        heterologous reaction database model.
    rxns_to_add : list
        contains reactions in COBRA format to be added to hrd_model.

    Returns
    -------
    None.

    """

    # load reaction IDs in the heterologous reaction database model
    hrd_model_rxns_list = []
    hrd_model_rxns = {}
    for rxn in hrd_model.reactions:
        hrd_model_rxns_list.append(rxn.id) 
        hrd_model_rxns[rxn.id] = rxn
    
       
    for rxn_id, rxn in rxns_to_add.items():
        if rxn_id in hrd_model_rxns_list:
            # add additional model origin
            hrd_model_rxns[rxn_id].origin.append(rxn.origin)
        else:
            # add reaction to reaction database
            hrd_model.add_reactions([rxn])
            # save reaction
            hrd_model_rxns[rxn_id] = rxn
         
        # extend database with origin models 
        #new_origin = list(set(rxn.origin).difference(set(hrd_model.heterologous_model_sources)))
        #hrd_model.heterologous_model_sources.extend(new_origin)
        # for origin in rxn.origin:
        #     if not(origin in hrd_model.heterologous_model_sources):
        #         hrd_model.heterologous_model_sources.append(origin)

    


def extract_heterologous_reactions_from_model(model, h_model, currency_mets=config_default.currency_mets,
                                              num_carbon_threshold=config_default.num_carbons_threshold,
                                              exclude_list=[],
                                              only_host_compartments=config_default.only_host_compartments):
    """Add cytosolic reactions model if meets all of the following criteria

    1) Not a boundary reaction
    2) Not the ATP synthase, biomass equation
    3) Has a gene reaction rule, and is not spontaneous
    4) Takes place entirely in one compartment
        - Membrane transport reactions are at high risk for loops
        - Compartment must be in chassis model (ie, c, p, and e for E. coli)
    5) Does not contain a metabolite with more carbons than `max_num_carbons`
        (lipids, etc.) that is not a currency metabolite (ATP, cofactors, etc)
    6) Does not include tRNA metabolites
    7) Is fully mass balanced
    8) not in host model
    9) Is in a compartment which exists in host
    10) Metabolite's elemental composition are consistent with host model
        - if formular for a metabolite is not in the model, query BIGG database and curate

    """
    
    
    # define parameter
    # gene placeholder for spontaneous reactions
    try:
        spont = h_model.genes.get_by_id("s0001")
    except:
        # no placeholder for spontaneous reaction in genes list
        spont = None
      
    # copied reaction
    copy_sign = "_copy"
      

    # get host model compartments
    if only_host_compartments:
        c_host = set(model.compartments.keys())
    else:
        c_host = set()
    
    
    rxns_to_add = []
    for rxn in h_model.reactions:
        
        # compartments exclusively in the host?
        if not(rxn.compartments.issubset(c_host)) and only_host_compartments:
            continue
         
        # is in host model?
        try:
            # reaction is already in host organism
            model.reactions.get_by_id(rxn.id)
            continue
        except:           
            pass    
                    
        # check if reaction is a copy
        if copy_sign in rxn.id:
            rxn_id = rxn.id[:rxn.id.index(copy_sign)]
        else:
            rxn_id = rxn.id
            
        
        # is in exclude list?
        if rxn_id in exclude_list:
            continue
        
        # is boundary reaction or ATP synthase?
        if rxn.boundary or rxn_id == "ATPS":
            continue
        
        
        # may be the biomass equation
        if "biomass" in rxn_id.lower():
            continue
        
        # is spontaneous or not related to a gene
        if not(rxn.genes) \
            or spont in rxn.genes \
                or "spontaneous" in rxn.name:
            continue
        
        # get carbon content of participating metabolites
        high_carbon = False
        cs = set()
        for met in rxn.metabolites:
            cs.add(met.id[-1])
            if met.id[:-2] in currency_mets:
                continue
            if met.elements.get('C', 0) > num_carbon_threshold:
                high_carbon = True
                break
        if high_carbon:
            continue
            
        # contains metabolites in different compartments?
        if len(cs) != 1:
            continue
        
        # is mass balance correct in other model
        mass_bal = rxn.check_mass_balance()
        if not(not(mass_bal)):
            continue
        
        
        # if metabolites are subset of host organism, check consistency
        is_inconsistent = False
        for met in rxn.metabolites:
            try:
                met_host = model.metabolites.get_by_id(met.id)
                # metabolite ID in host organism, check elements, if annotated
                if not(met_host.formula) or not(met.formula):
                    # formula of metabolite not provied
                    continue             
                
                for element, number in met.elements.items():
                    if element in list(met_host.elements.keys()):
                        if number != met_host.elements[element]:
                            # number of elements does not match                     
                            is_inconsistent = True                            
                    else:
                        # elements do not match
                        is_inconsistent = True
                    
                if is_inconsistent:
                    print(met.id + " in " + rxn_id + " is inconsistent")
                    break
                    
            except:
                # metabolite not in host model, continue
                pass
        if is_inconsistent:
            continue
        
        
        
        

        # is mass balance correct in host model?
        with model:
            # add heterologous reaction
            model.add_reactions([rxn.copy()])
            # check mass balance
            mass_balance_host = model.reactions.get_by_id(rxn.id).check_mass_balance()
            if not(mass_balance_host):
                # mass balance is OK
                pass               
            else:               
                # mass balance fails, check formula of metabolites
                # print(model.reactions.get_by_id(rxn.id).metabolites.items())
                for met in model.reactions.get_by_id(rxn.id).metabolites.keys():
                    met_id = met.id

                    if not(met.formula):
                        # get formula from BIGG
                        
                        bigg_query = bigg.get(met_id[:-2], "metabolites", "universal")
                        if not(bigg_query):
                            # no entry in BIGG database found
                            print(met_id[:-2] + " not found in BIGG database!")
                        else:
                            # add formula and charge into host and other model
                            if not(bigg_query["formulae"]):
                                # no formula provided
                                continue
                            else:
                                met.formula = bigg_query["formulae"][0]
                                h_model.metabolites.get_by_id(met_id).formula = bigg_query["formulae"][0]
                                
                            if not(bigg_query["charges"]):
                                # no charges provided
                                pass
                            else:
                                met.charge = bigg_query["charges"][0]
                                h_model.metabolites.get_by_id(met_id).charge = bigg_query["charges"][0]
            
                            print("Formula of " + met_id + " added: " + bigg_query["formulae"][0])
                        
                # check mass balance again
                mass_balance_host = model.reactions.get_by_id(rxn.id).check_mass_balance()
                if not(not(mass_balance_host)):  
                    print("Mass balance fails in host organism for: " + rxn.id)
                    continue


        
        # rxns_to_add[rxn_id] = h_model.reactions.get_by_id(rxn.id).copy()

        rxns_to_add.append(rxn.id)
        # change reaction ID if truncated
        # rxns_to_add[rxn_id].id = rxn_id
    
        
        
    return rxns_to_add
        

def evaluate_reaction_directionalities(model, rxns_to_add_in, config=config_default):
    """
    evaluate reaction directions based on thermodynamics, flux variabilities, and annotations

    Parameters
    ----------
    model : cobra.Model
        COBRA model.
    rxns_to_add_in : list
        DESCRIPTION.
    config : TYPE, optional
        DESCRIPTION. The default is config_default.

    Returns
    -------
    rxns_to_add : TYPE
        DESCRIPTION.
    directions : TYPE
        DESCRIPTION.

    """

    # import thermodynamic evaluation modules
    # do it here to avoid loading Equilibrator if this method here is not used
    from growth_coupling_suite.reaction_directionality \
        import reaction_thermodynamics_evaluation as rt
    
    # check config file
    check_config(config, config_default)
    
    # copy rxn_to_add_in
    rxns_to_add = rxns_to_add_in.copy()
    
    # save direction in data frame
    directions = pd.DataFrame()
    
    # thermodynamic and annotation assessment of reaction directions
    print("Assess directions of heterologous reactions...")
    rxns_to_add_iter = rxns_to_add.copy()
    for rxn_id, rxn in rxns_to_add_iter.items():
        print(rxn_id)
        # get direction by thermodynamic assessment   
       
        try:
             direction_thermodynamic, eq_rxn, reversibility_measure \
                = rt.get_reaction_direction(rxn, config=config)           
        except: #(TypeError) as e:
            # print(e.args[0])
            # thermodynamic evaluation failed
            print("Determination of reaction direction failed")
            direction_thermodynamic = -3
          
        # get annotated direction
        direction_annotated = get_annotated_reaction_direction(rxn)
        
        print(str(direction_thermodynamic) + " " + str(direction_annotated))
        
        if direction_thermodynamic==-3 and \
            (config.consider_thermodynamic_valid_reactions_only or direction_annotated==-3):
                # unable to define a reaction direction or thermodynamically invalid not allowed
                del(rxns_to_add[rxn_id])
                continue
        else:
            directions = pd.concat([directions, pd.DataFrame({
                "thermodynamic": direction_thermodynamic,
                "annotated": direction_annotated},
                index=[rxn_id])])
            
            # set reaction direction
            if abs(direction_thermodynamic) == 1:
                # thermodynamically irreversible reaction
                set_reaction_direction(rxns_to_add[rxn_id], direction_thermodynamic)
            elif direction_thermodynamic == 0:
                # thermodynamically reversible reaction
                if abs(direction_annotated) == 1:
                    # annotation suggests irreversibility
                    set_reaction_direction(rxns_to_add[rxn_id], direction_annotated)
                else:
                    # reversible reaction
                    set_reaction_direction(rxns_to_add[rxn_id], 0)
            else:
                # inconclusive thermodynamic assessment
                set_reaction_direction(rxns_to_add[rxn_id], direction_annotated)
                    


    return rxns_to_add, directions


def consolidate_directionality_results(directions, directions_fva):
    """
    consolidate reaction directions suggestions from thermodynamic and FVA results
    
    Heuristic (restrictive, favors irreversible reactions):
        - reaction is blocked for any condition according to FVA remove reaction
        - thermodynamic analysis suggests irreversibility
            - if reaction needs to operate in the opposite direction for one or more
              condiions accordsing to FVA, make reaction reversible
             - otherwise set reaction to irreversible
        - thermodynamic analysis suggest reversible reaction
            - if FVA shows uniform direction under any condition, make reaction irreversible
            - if reaction is annotated as irrversible, make reaction irreversible
        - thermodynamic resuls are inconclusive
            - if FVA suggest irreversible reaction, set reaction as irreversible
        
             
    
    
    """
    
    directions_dict_out = {}
    
    # consolidate reaction IDs
    rxn_list = list(set(directions.index) & set(directions_fva.index))
    
    for rxn_id in rxn_list:
        # load reaction directionalities
        rxn_dir_fva = directions_fva.loc[rxn_id]
        rxn_dir_thermodynamic = directions.loc[rxn_id, "thermodynamic"]
        rxn_dir_annotated = directions.loc[rxn_id, "annotated"]
        
        # is reaction blocked?
        if np.all(rxn_dir_fva == -2) or np.all(rxn_dir_fva == -3):
            # disregard reaction
            directions_dict_out[rxn_id] = -2
        
        # irreversible according to thermodynamic analysis?
        elif abs(rxn_dir_thermodynamic) == 1:
            # FVA suggests operation in reverse direction?
            if np.any(rxn_dir_fva == -rxn_dir_thermodynamic):
                # make reaction reversible
                directions_dict_out[rxn_id] = 0
            else:
                # make reaction irreversible
                directions_dict_out[rxn_id] = rxn_dir_thermodynamic
                
        # reversible according to thermodynamics?
        elif rxn_dir_thermodynamic == 0:
            # FVA suggest irrversibility?
            if np.all(rxn_dir_fva == 1):
                directions_dict_out[rxn_id] = 1           
            elif np.all(rxn_dir_fva == -1):
                directions_dict_out[rxn_id] = -1 
            # does annotation suggest an irreversible reaction?
            elif abs(rxn_dir_annotated) == 1:
                directions_dict_out[rxn_id] = rxn_dir_annotated
            else:
                directions_dict_out[rxn_id] = 0
         
        # no thermodynamic results        
        else:
            # FVA suggest irrversibility?
            if np.all(rxn_dir_fva == 1):
                directions_dict_out[rxn_id] = 1           
            elif np.all(rxn_dir_fva == -1):
                directions_dict_out[rxn_id] = -1 
            else:
                directions_dict_out[rxn_id] = rxn_dir_annotated
                
        
    
    
    return directions_dict_out



def get_annotated_reaction_direction(rxn):
    """
    Assess reaction direction from reaction bounds

    Parameters
    ----------
    rxn : cobra.Reaction
        reaction in COBRA format.

    Returns
    -------
    direction : int
        depicts reaction direction
            0: reversible
            1: irreversible in forward direction
            -1: irreversible in backward direction
            -2: blocked.
            -3: assessment failed

    """
    # analyze reaction direction based on flux bounds
    if rxn.lower_bound < 0 and rxn.upper_bound > 0:
        direction = 0
    elif rxn.lower_bound >= 0 and rxn.upper_bound > 0:
        direction = 1
    elif rxn.lower_bound < 0 and rxn.upper_bound <= 0:
        direction = -1
    elif rxn.lower_bound == 0 and rxn.upper_bound == 0:
        direction = -2
    else:
        print("Warning: Unknown reaction direction for " + rxn.id)
        direction = -3
        
    return direction


def set_reaction_direction(rxn, direction):
    """
    Set reaction bounds according to specified direction

    Parameters
    ----------
    rxn : cobra.Reaction
         reaction in COBRA format.
    direction : int
        depicts direction of reaction.

    Returns
    -------
    None.

    """
    # set boundaries of reaction according to specified direction
    if direction == 1:
        rxn.lower_bound = 0
        rxn.upper_bound = 1000
    elif direction == 0:
        rxn.lower_bound = -1000
        rxn.upper_bound = 1000
    elif direction == -1:
        rxn.lower_bound = -1000
        rxn.upper_bound = 0
    else:
        print("Warning: Unknown reaction direction: " + str(direction))
        
        
def load_heterologous_database_model(model_name, model_path=hrd_dir, reaction_direction_assessed=True):
    """
    Load a heterologous reaction database model

    Parameters
    ----------
    model_name : str
        Name of the heterologous reaction database model.
    model_path : str, optional
        Path to the heterologous reaction database model. The default is hrd_dir.
    reaction_direction_assessed : boolean, optional
        indicates if directions of reaction in the model have beena assessed. The default is True.

    Returns
    -------
    hrd_model : cobra.Model
        Heterologous reaction database model in cobra format.

    """
    
    # load (existing) heterologous reaction database model and load heterologous models

    # check if heterologous reactions database exist for model, directions not assessed
    if reaction_direction_assessed:
        hrd_model_feature = "_dir_assessed"
    else:
        hrd_model_feature = ""
    
    model_dir = Path(model_path).joinpath(model_name + hrd_model_feature + ".json")
    if model_dir.is_file():
        print("Load existing heterologous reaction database model...\n" +
              "\t" + model_name + hrd_model_feature)

        hrd_model = cobra.io.load_json_model(str(model_dir))
        
        
    else:
        print("No heterologous reaction database model found for",
              model_name + hrd_model_feature)
        hrd_model = []
        
    return hrd_model

        
def save_heterologous_database_model(model, model_name, model_path=hrd_dir, reaction_direction_assessed=True):   
    """
    Save an heterologous reaction database model

    Parameters
    ----------
    model : cobra.Model
        Heterologous reaction database model in cobra format.
    model_name : str
        Model name.
    model_path : str, optional
        Path to the heterologous reaction database model. The default is hrd_dir.
    reaction_direction_assessed : TYPE, optional
        Indicates if directions of reaction in the model have beena assessed. The default is True.

    Returns
    -------
    None.

    """
    
    # check if heterologous reactions database exist for model, directions not assessed
    if reaction_direction_assessed:
        hrd_model_feature = "_dir_assessed"
    else:
        hrd_model_feature = ""

    model_dir = Path(model_path).joinpath(
        model_name + hrd_model_feature
        )
        
    # save in .json format
    cobra.io.save_json_model(model, str(model_dir)+".json")

