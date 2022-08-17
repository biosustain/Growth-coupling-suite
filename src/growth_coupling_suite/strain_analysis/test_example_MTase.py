from strain_design_analysis import StrainDesignAnalyzer
from growth_coupling_suite import models
from growth_coupling_suite.heterologous_reactions_processing import heterologous_reactions_database
import cobra
import pickle
import json

from os.path import abspath, dirname

#%% load model
model_dir = dirname(abspath(models.__file__))
model_name = "iML1515"
model = cobra.io.load_json_model(f"{model_dir}\{model_name}.json")


# %% add MTase reaction
rxn = cobra.Reaction("MTase_generic")
rxn.lower_bound = 0
rxn.upper_bound = 1000
# get metabolites from model
amet = model.metabolites.get_by_id("amet_c")
ahcys = model.metabolites.get_by_id("ahcys_c")
# create new metabolites
R = cobra.Metabolite("R", formula = "", name="generic compound", compartment="c", charge=0)
R_meth = cobra.Metabolite("R_meth", formula = "CH3", name="generic compound methylated", compartment="c", charge=1)
# add metabolites
rxn.add_metabolites({
    amet: -1,
    ahcys: 1,
    R: -1,
    R_meth: 1
})


# create sinks and sources for generic metabolites
R_source = cobra.Reaction("R_source")
R_source.lower_bound = -1000
R_source.upper_bound = 0
R_source.add_metabolites({R: -1})
R_meth_sink = cobra.Reaction("R_meth_sink")
R_meth_sink.lower_bound = 0
R_meth_sink.upper_bound = 1000
R_meth_sink.add_metabolites({R_meth: -1})


# add reaction to model
model.add_reactions([rxn, R_source, R_meth_sink])
# add sinks and sources for generic metabolites

# check model
with model:
    model.objective = "MTase_generic"
    sol = model.optimize()
    print(sol.objective_value)
    
    
    
# %% literature strain design
# Luo et al. 2019
solutions = {}
sol_lit = {"interventions": {"cysE": {"ID": "SERAT",
                                     "type": "deletion",
                                     "lower_bound" : 0,
                                     "upper_bound": 0},
                            "CYSTL": {"ID": "CYSTL",
                                     "type": "addin",
                                     "lower_bound": -1000,
                                     "upper_bound": 1000},
                            "CYSTGL": {"ID": "CYSTGL",
                                      "type": "addin"},
                            "EX_met__L_e": {"ID": "met__L_e",
                                   "type": "cofeed",
                                   "lower_bound": -10,
                                   "upper_bound": 0}}}
solutions["sol_lit"] = sol_lit

# %% Computed strain designs
res_dir = "C:/Users/Tobi/Bull/NAS Drive/Postdoc CFB/Projekte/01 Growth coupling selection systems/05_Growth coupling suite/MTase example/results/g01_i6_d5_a5_c1"
with open(res_dir + "/callback_solutions_dict_MTase_generic_g01_i6_d5_a5_c1.json", "r") as f:
    solutions_comp = json.load(f)
    
solutions = {**solutions, **solutions_comp}


# %% load heterologous reactions
hrd_dir = abspath(dirname(heterologous_reactions_database.__file__))
with open(hrd_dir+"/iML1515_hr_database_dir_assessed.pickle", "rb") as f:
    hrd = pickle.load(f)

# %% load strain design analyzer
sd = StrainDesignAnalyzer(model,
                            design_solutions_dict=solutions,
                            heterologous_reaction_database_model=hrd,
                            save_model=False)



# apply design
target_rxn = "MTase_generic"
# prepare model
sd.model.exchanges.EX_glc__D_e.lower_bound = -10
sd.model.exchanges.EX_lac__D_e.lower_bound = 0
# apply solution
sd.solutions.sol_lit.load()

# %% investigate design
with sd.model as m:
    # rescue growth-coupling
    m.reactions.MTase_generic.bounds = (0, 0)
    # m.reactions.EX_lac__D_e.lower_bound = -0.0546
    # m.reactions.EX_co2_e.lower_bound = 0
    # m.reactions.source_met__L_e.lower_bound = -0.1
    # cofeed metabolites
    # m.exchanges.EX_cys__L_e.lower_bound = 0
    # create source
    source = cobra.Reaction("source_rxn",
                             lower_bound=-10,
                             upper_bound=0)
    met = m.metabolites.get_by_id("cpe180_c")
    source.add_metabolites({met: -1})
    # m.add_reaction(source)
    
    sol = m.slim_optimize()
    # calculate flux space
    sd.flux_space_projection(target_rxn)
    # calculate auxotrophies
    # aux_mets = sd.determine_auxotrophies()
    # calculate precursor availability
    precursor_fluxes, target_rxn_rescue_flux = sd.precursor_availability()
    # determine precursor bottleneck
    target_flux_after_removal = sd.precursor_bottleneck()


# coupling strength
# print(sd.coupling_strength(target_rxn))
# sd.flux_space_projection(target_rxn)
# a = sd.flux_bound_limitation_changes(target_rxn,
                                      # objective_direction="min",
                                      # grid_size=20,)


