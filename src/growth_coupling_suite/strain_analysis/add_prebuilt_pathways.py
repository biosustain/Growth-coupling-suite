# add custom made and pre-built pathways

from cobra import Reaction, Metabolite



# %% add tetrahydrobiopterin regeneration
def tetrahydrobiopterin(model):
    
    
    # add Tetrahydrobiopterin-4a-carbinolamine dehydratase
    if "THBPT4ACAMDASE" not in model.reactions:
       
        THBPT4ACAMDASE = Reaction(id="THBPT4ACAMDASE",
                             lower_bound=0,
                             upper_bound=1000)
        try:
            dhbpt_c = model.metabolites.get_by_id("dhbpt_c")
        except:
            dhbpt_c = Metabolite(id="dhbpt_c",
                                       compartment="c",
                                       formula="C9H13N5O3",
                                       charge=0)
            
        try: 
            thbpt4acam_c = model.metabolites.get_by_id(" thbpt4acam_c")
        except:
            thbpt4acam_c = Metabolite(id="thbpt4acam_c",
                                       compartment="c",
                                       formula="C9H15N5O4",
                                       charge=0)  
            
        THBPT4ACAMDASE.add_metabolites({
            thbpt4acam_c: -1,
            model.metabolites.get_by_id("h2o_c"): 1,
            dhbpt_c: 1
            })
        
        model.add_reactions([THBPT4ACAMDASE])
    
    
    # add 6,7-dihydropteridine reductase
    if "DHPR" not in model.reactions and "DHPR2" not in model.reactions:
        
        DHPR = Reaction(id="DHPR",
                        lower_bound=0,
                        upper_bound=1000)
        
        try:
            thbpt_c = model.metabolites.get_by_id("thbpt_c")
        except:
            thbpt_c = Metabolite(id="thbpt_c",
                                    compartment="c",
                                    formula="C9H15N5O3",
                                    charge=0)
            
        DHPR.add_metabolites({
            model.metabolites.get_by_id("h_c"): -1,
            model.metabolites.get_by_id("nadh_c"): -1,
            model.metabolites.get_by_id("dhbpt_c"): -1,
            model.metabolites.get_by_id("nad_c"): 1,
            thbpt_c: 1,
            })
        
        model.add_reactions([DHPR])
      
    
    
    return model