"""
Process COBRA models for strain design optimization purposes
"""

from growth_coupling_suite.model_processing import config_model_processing_default as config_default

from growth_coupling_suite.heterologous_reactions_processing \
    import heterologous_reactions_processing as hrp

from growth_coupling_suite.util.util import check_config

import cobra



def add_heterologous_reactions(model, hr_database_model=None, 
                               config=config_default):
    """
    Extract, process, and add heterologous reactions to a (host) model
    Heterologous reactions will be marked with a "__hr" in the ID

    Parameters
    ----------
    model : cobra.model
        COBRA model.
    hr_database_model : cobra.model
        Heterologou reaction database model. The default is None   
    config : module, optional
        Optional parameters. The default is config_default.

    Returns
    -------
    model : cobra.model
        COBRA model with integrated heterologous reactions.
    heterologous_reactions : list
        IDs of all added heterologous reactions.

    """
    
    # parameter
    rxn_id_format = "{0}__hr"
    rxn_name_format = "{0} heterologous"
    
    
    # was a heterolgous reaction database model provided?
    if not(hr_database_model):
        # get heterologous reaction from database 
        with model:
            hr_database_model_out, hr_database_model_origin_out = hrp.get_heterologous_reactions(model,
                                                          config=config,
                                                          reprocess=False)
    
        # choose assessed or unassessed heterologous database
        if hasattr(config, "directionality_assessment"):
            if config.directionality_assessment:
                # directionality assessed
                hr_database_model = hr_database_model_out
            else:
                # directionality not assessed
                hr_database_model = hr_database_model_origin_out
                
        else:
            # directionality is assessed by default
            hr_database_model = hr_database_model_out

    # add and rename heterologous reactions to the original model
    heterologous_reactions_dict = {}
    heterologous_reactions = []
    rxns_to_add = []
    for rxn in hr_database_model.reactions:
        
        # reaction explicitly excluded?
        if rxn.id in config.addin_exclude_list:
            continue
        
        # add reaction to model, create from scratch to avoid name assignment problems
        rxn_to_add = cobra.Reaction(rxn_id_format.format(rxn.id),
                                    name=rxn_name_format.format(rxn.name),
                                    lower_bound=rxn.lower_bound,
                                    upper_bound=rxn.upper_bound)

        # add metabolites
        mets_to_add ={}
        for met, coeff in rxn.metabolites.items():
            mets_to_add[met] = coeff
        rxn_to_add.add_metabolites(mets_to_add)
        
        rxns_to_add.append(rxn_to_add)
        
        # save heterologous reaction ID
        heterologous_reactions.append(rxn_to_add.id)
      
    # add heterologous reaction to the model
    model.add_reactions(rxns_to_add)
    
    for rxn_id in heterologous_reactions:
        heterologous_reactions_dict[rxn_id] = model.reactions.get_by_id(rxn_id)

    return model, heterologous_reactions, hr_database_model



def get_deletion_target_space(model, config=config_default):
    """
    Finds the reactions in the chassis model that can be knocked out.

    Must meet following criteria:
    1) Must have a GPR and not spontaneous (s0001 in GPR)
    2) Not have a subsystem in a user-defined exclusion list
    3) Not be in a user-defined reaction exclusion list
    4) Only have metabolites in one compartment
        Knocking out transport reactions can be problematic
    5) Have all metabolites with fewer than 35 carbons
        Unless the metabolite is a currency metabolite
      

    Parameters
    ----------
    model : cobra.Model
        Metabolic model in COBRA format.
    config : module, optional
        optimization and model parameters. The default is config_default.

    Returns
    -------
    deletion_targets : list
        IDs of reaction deletion targets.

    """

    
    # load config
    check_config(config, config_default)

    
       
    deletion_targets = []   
    
    try:
        spont = model.genes.get_by_id('s0001')
    except:
        spont = None
            
    for rxn in model.reactions:
        # is reaction specified or excluded by user as target?
        if rxn.id in config.deletion_targets:
            # accept reaction right away
            deletion_targets.append(rxn.id)
            continue
        elif rxn.id in config.deletion_exclude_list:
            # exclude reaction right away
            continue
        
        # is boundary reaction?
        if rxn.boundary:
            continue
        
        # gene related?
        if not(rxn.genes):
            continue
        
        # exclude by subsystem?
        if rxn.subsystem in config.subsystem_exclude_list:
            continue
        
        # is reaction spontaneous?
        if spont in rxn.genes or "spontaneous" in rxn.name:
            continue
        
        # critical exchange reaction?
        if rxn.id in config.exchanges_not_to_knockout:
            continue
        
        # excluded by the user
        if rxn.id in config.user_exclude_list:
            continue
        
        # reaction acts on high carbon metabolites?
        if not(not(config.num_carbons_threshold)):
            high_carbon = False
            cs = set()
            for met in rxn.metabolites:
                cs.add(met.id[-1])
                if met.id[:-2] in config.currency_mets:
                    continue
                if met.elements.get('C', 0) > config.num_carbons_threshold:
                    high_carbon = True
                    break
        if high_carbon:
            continue
            
        # contains metabolites in different compartments?
        if not(config.allow_transporter_deletion) and (len(cs) != 1):
            continue
        
        # accept reaction
        deletion_targets.append(rxn.id)
      
    return deletion_targets



def set_source_targets(model,
                       allow_carbon_source_switch=False,
                       allow_non_carbon_sources=True,
                       allow_cofeed_targets=True,
                       config=config_default):
    """
    introduce source reactions for carbon substrates  or cofeeds
    Do NOT use existent exchange reactions as sources
      - the target metabolite may be secreted in designs where it is not a substrate
      - split exchange reaction into forward and backward reaction
    exclude certain uptake reactions, e.g., for CO2

    Parameters
    ----------
    model : cobra.Model
        Metabolic model in COBRA format.
    allow_carbon_source_switch : bool, optional
        Enable changing of carbon sources as target variable. The default is False.
    allow_non_carbon_sources : bool, optional
        Enable changing of non-carbon sources as target variable. The default is True.
    allow_cofeed_targets : bool, optional
        Enable the addition of a cofeed as target variable. The default is True.
    config : module, optional
        Parameters for setting source or cofeed targets. The default is config_default.

    Returns
    -------
    model : cobra.Model
        Metabolic model including source or cofeed target reactions.
    source_targets : list
        IDs of source reactions in the model
    cofeed_targets : list
        IDs of cofeed reactions in the model

    """

    
    # define functions
    def create_source_reaction(met, lb, id):
        # create reaction
        rxn_source = cobra.Reaction(id=id,
                                    lower_bound=lb,
                                    upper_bound=0)
        # add metabolite
        rxn_source.add_metabolites({met: -1})
        # set bounds
        # rxn_source.bounds = (lb, 0)
        return rxn_source 
        
    
    
    # check config
    check_config(config, config_default)    
    
    # merge "to keep" reaction lists
    exchanges_to_keep = config.exchanges_to_keep
    exchanges_to_keep.extend(config.cofeed_targets)
    exchanges_to_keep.extend(config.source_targets)
    
    # exclusively take user defined source targets if specified by source exclude list
    if config.source_exclude_list==-1 and len(config.source_targets)>0:
         source_exclude_list = [ex.id for ex in model.boundary if ex.id not in config.source_targets]
    else:
         source_exclude_list = config.source_exclude_list
        
    # exclusively take user defined cofeed targets if specified by cofeed exclude list
    if config.cofeed_exclude_list==-1 and len(config.cofeed_targets)>0:
         cofeed_exclude_list = [ex.id for ex in model.boundary if ex.id not in config.cofeed_targets]
    else:
         cofeed_exclude_list = config.cofeed_exclude_list   
                
    # scan exchange reactions
    cofeed_targets = []
    source_targets = []
    for ex in model.boundary:
        
        mets = ex.metabolites
        # get metabolite
        met = list(mets.keys())[0] 
        # is exchange reaction explicitly included?
        if ex.id not in exchanges_to_keep:           
            # exclude exchange reaction?
            if ex.id in config.exchanges_not_to_add:
                # exchange reaction will not be deleted and replaced by a variable source reaction
                continue
            # generally excluded
            elif ex.id in config.user_exclude_list:
                continue
            # is a boundary reaction (only one metabolite participating)?
            elif len(mets) !=1:
                continue
                       
            # too high carbon content?
            # metabolite's element list
            if not(list(met.elements.keys())):
                print("Elemental composition not specified for exchange metabolite " + met.id)
                print('\t Exchange reaction', ex.id, 'will not be considered as source/cofeed target')
                continue
            
            # calculate carbon content
            if "C" in list(met.elements.keys()):
                if met.elements["C"] > config.num_exchange_carbons_threshold:
                    continue
                else:
                    c_content = met.elements["C"] 
            elif allow_non_carbon_sources:
                c_content = 0
            else:
                continue
            
            # calculate nitrogen content
            if "N" in list(met.elements.keys()):         
                n_content = met.elements["N"]
            else:
                n_content = 0
            
            # is excretion enforced?
            if ex.lower_bound > 0:
                continue
            
            # is exchange already an uptake reaction or is flux enforced?
            elif ex.lower_bound < 0:
                if c_content == 0 or not(allow_carbon_source_switch):
                    # not a carbon substrate or carbon source switch not allowed
                    continue
                else:
                    # enable switch of carbon uptake flux, disable original exchange     
                    ex.lower_bound = 0
                    
        else:
            # calculate carbon content
            if "C" in list(met.elements.keys()):         
                c_content = met.elements["C"]
            else:
                c_content = 0
                
            # calculate nitrogen content
            if "N" in list(met.elements.keys()):         
                n_content = met.elements["N"]
            else:
                n_content = 0
                
            # is exchange already an uptake reaction or is flux enforced?
            if ex.lower_bound < 0 and allow_carbon_source_switch:
                # enable switch of carbon uptake flux, disable original exchange     
                ex.lower_bound = 0    
                  
               

        # define uptake bound
        uptake_bound_c = 1000
        uptake_bound_n = 1000
        
        if c_content > 0:
            # carbon substrate
            uptake_bound_c = config.maximum_source_carbon_flux/c_content
            
        if (n_content > 0) and config.maximum_source_nitrogen_flux:
            uptake_bound_n = config.maximum_source_nitrogen_flux / n_content
        
        if (c_content == 0) and ((n_content == 0) or (not(config.maximum_source_nitrogen_flux))):
            # no carbon or nitrogen in metabolite, apply absolute general flux bound
            molar_mass = met.formula_weight
            if not molar_mass:
                # no molar mass provided             
                uptake_bound = config.maximum_source_non_carbon_flux
            else:
                uptake_bound = config.maximum_source_non_carbon_mass_flux/(molar_mass/1000)
                
        else:
            # take most restrictive uptake bound
            if uptake_bound_c >= uptake_bound_n:
                uptake_bound = uptake_bound_n
            else:
                uptake_bound = uptake_bound_c
         
        # create new exchange (uptake) reaction for metabolite target
        
        # create source metabolite that the resembles the exchange metabolite
        met_source = cobra.Metabolite(id=met.id+"_src", name=met.name,
                                      charge=met.charge, formula=met.formula,
                                      compartment="src")
        
        # create source to exchange (diffusion) reaction
        source_to_exchange = cobra.Reaction(id="EX2SRC_"+met.id,
                                            lower_bound=-uptake_bound,
                                            upper_bound=0)  
        source_to_exchange.add_metabolites({
            met_source: 1,
            met: -1})
        model.add_reaction(source_to_exchange)

        
        
        # add a cofeed target reaction?
        if allow_cofeed_targets:
            # is metabolite on exclude list target or explicitly enforced?
            if ((ex.id in config.cofeed_targets) or (ex.id not in cofeed_exclude_list))\
                and (ex.id not in config.source_targets):         
                # create and add source reaction for metabolite
                cofeed_reaction = create_source_reaction(met_source, -uptake_bound*config.cofeed_rate_factor, "cofeed_" + met.id)
                model.add_reaction(cofeed_reaction)
                # print(model.reactions.get_by_id(source_reaction.id))
                # save target
                cofeed_targets.append(cofeed_reaction.id)
                # cofeed_database_model.add_reactions([source_reaction])
            
        # add a carbon source target reaction?
        if allow_carbon_source_switch:
            # is metabolite on exclude list target or explicitly enforced?
            if ((ex.id in config.source_targets) or (ex.id not in source_exclude_list))\
                and (ex.id not in config.cofeed_targets): 

                # create and add source reaction for metabolite
                source_reaction = create_source_reaction(met_source, -uptake_bound, "source_" + met.id)
                model.add_reaction(source_reaction)
                # print(model.reactions.get_by_id(source_reaction.id))
                # save target
                source_targets.append(source_reaction.id)
                # cofeed_database_model.add_reactions([source_reaction])
    
   
    return model, source_targets, cofeed_targets
      


def get_objective_reaction(model):
    """
    return reactions embedded in the objective function

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.

    Returns
    -------
    objective_reactions : list
        IDs of reactions embedded in the objective function.

    """
    # get reaction of the objective function, normally the biomass formation reaction
    
    # protect model
    with model:
        objective_reactions = []
        for rxn in model.reactions:
            # check objective coefficient of reaction
            if rxn.objective_coefficient != 0:
                objective_reactions.append(rxn.id)
                
    return objective_reactions
                
            
            