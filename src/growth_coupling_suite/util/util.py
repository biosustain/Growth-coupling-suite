""" Utility functions """
import importlib
import numpy as np
from warnings import warn

def check_config(config, config_default=None):
    """
    check config module, if attributes are missing use default values
    If no default config is provided check against all default config files available

    Parameters
    ----------
    config : module
        includes config attributes.
    config_default : module
        default config module against which config is checked.

    Returns
    -------
    None.

    """

    if not(config_default):
        # compare with all default config files of the GC suite
        
        # import default config modules
        module_base = "growth_coupling_suite.{0}"
        modules_name = [
            "gcOpt_algorithm.config_gcOpt_default",
            "heterologous_reactions_processing.config_heterologous_reactions_default",
            "model_processing.config_model_processing_default",
            "reaction_directionality.config_reaction_directionality_default"
            ]
        modules = {}
        for m in modules_name:
            modules[m] = importlib.import_module(module_base.format(m))
            
        # compare attributes with default config
        for m, default_config_module in modules.items():
            config_attr = vars(default_config_module)
            for attr, value in config_attr.items():
                if not(hasattr(config, attr)):
                    setattr(config, attr, value) 
                
    else:
       
        # compare attributes with default config
        config_attr = vars(config_default)
        for attr, value in config_attr.items():
            if not(hasattr(config, attr)):
                setattr(config, attr, value)             
  
                       
def is_reaction_in_model(model, rxn_id):
    """
    Check if reaction is in COBRA model

    Parameters
    ----------
    model : cobra.model
        COBRA model.
    rxn_id : string
        reaction identifier.

    Returns
    -------
    in_model : boolean
        Depicts if reaction 'rxn_id' is in model.

    """
    # check if reaction is in model
    try:
        model.reactions.get_by_any(rxn_id)
        in_model = True
    except:
        in_model = False   
        
    return in_model


def GPR_linked_reaction_knockouts(model, rxn_ids: list, eval_gpr: bool = False,
                                      gpr_expressions: dict = {},
                                      gene_stat_dict: dict = {}) -> list:
    """
    Determine all reactions which are affected by a gene-based knockout of a set of reactions
    
    Parameters
    ----------
    model : cobra.core.model
        metabolic model in COBRA format.
    rxn_ids : list
        Reaction IDs.

    Returns
    -------
    dependent_rxns : List of reactions which may be affected by gene-based
                     knockout of rxn_id. rxn_id is NOT included itself.

    """
    
    dependent_rxns = []
    if not(eval_gpr):
        # define dependent reaction knockouts by any gene-based link
        # do not consider GPR logical expressions

        for rxn_id in rxn_ids:
            # check existence of reaction
            if rxn_id not in model.reactions:
                warn(rxn_id + " not in model. Could not find GPR-related reaction knockouts", UserWarning)
                continue

            
            # determine genes and GPR of reaction
            genes = model.reactions.get_by_id(rxn_id).genes
            
            # determine all reactions dependent on extracted genes
            if len(genes) > 0:
                for g in genes:
                    for r_g in g.reactions:
                        if (r_g.id not in rxn_ids) and (r_g.id not in dependent_rxns):
                            dependent_rxns.append(r_g.id)
                            
    else:
        # find dependent reaction knockouts according to GPR logical expressions
        
        # create GPR expressions if not provided
        if not(gpr_expressions):
            gpr_expressions, gene_stat_dict = create_GPR_expressions(model)
        # create gene status dict if not provided
        if not(gene_stat_dict):
            gene_stat_dict = {g.id: True for g in model.genes}
            
        # determines genes which need to be knocked out to disable reactions 
        for rxn_id in rxn_ids:
            # check existence of reaction
            if rxn_id not in model.reactions:
                warn(rxn_id + " not in model. Could not find GPR-related reaction knockouts", UserWarning)
                continue
            
            for g in model.reactions.get_by_id(rxn_id).genes:
                # assign gene knockout in gene status dict
                gene_stat_dict[g.id] = False
                
        # evaluate GPR expressions and determine reaction knockouts
        for r, expr in gpr_expressions.items():
            if not(eval(expr)) and (r not in rxn_ids):
                # reaction is disabled due to gene knockouts
                dependent_rxns.append(r)
            
                
             
                                        
    return dependent_rxns


def create_GPR_expressions(model):
    """
    Transform gene (protein) reaction rules (GPR) into machine readable expressions

    Parameters
    ----------
    model : TYPE
        DESCRIPTION.

    Returns
    -------
    gpr_expressions : dict
        machine executable GPR rules. Keys are reaction IDs
    gene_stat_dict : dict
        Status of genes. True: active, False: knockout

    """
    
    # helper functions
    def transform_gpr(gpr, gene_ids):
        # transform GPR into a machine readable form
        # set logical operators
        gpr_expr = gpr.replace("and", "&").replace("or", "|")
        # replace gene ids with calls to gene status array
        for g in gene_ids:
            gpr_expr = gpr_expr.replace(g, "gene_stat_dict['" + g + "']")

        return gpr_expr        
    
    # allocate gene status array
    gene_stat_dict = {g.id: True for g in model.genes}

    
    gpr_expressions = {}  
    for r in model.reactions:
        gpr = r.gene_reaction_rule
        if len(gpr) == 0: continue # no GPR available for reaction
        
        # determine IDs of genes involved in GPR
        gene_ids = [g.id for g in r.genes]
        # transform GPR
        gpr_expressions[r.id] = transform_gpr(gpr, gene_ids)
    
    
    
    return gpr_expressions, gene_stat_dict




    


