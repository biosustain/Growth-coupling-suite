"""
gcOpt config file
For more optional attributes refer to the default config files in the growth coupling suite directories
"""

import os


# %% mandatory gcOpt configurations
output_dir = os.getcwd() + "/results_parallel" # directory in which results are saved
output_file_id = "succinate_gc"

biomass_id = "BIOMASS_Ec_iJO1366_core_53p95M" # BIOMASS_Ec_iML1515_core_75p37M,  BIOMASS_Ec_iJO1366_core_53p95M
growth_rate_fix = 0.3   # fix biomass value


# %% main gcOpt configurations

# solver specific
time_limit = 60 # runtime limit of (gurobi) solver [sec]
processes = 4 # number of threads gurobi is using
node_memory_limit = None # memory limit (in gb) for each node from which it will be saved on the harddrive


# specific number of interventions
num_total_interventions = 5 # number of total interventions

num_deletions = 4 # maximum number of reaction deletions
num_addins = 4 # maximum number of heterologous reaction additions
num_cofeeds = 1 # maximum number of additional metabolite cofeeds
num_carbonsources = 1 # maximum number of carbon sources

# evaluate GPR relations for each design solution
eval_gpr = True

# evaluate growth phenotype for each solution
# Only solutions allowing for a growth rate greater than eval_max_growth are valid
eval_max_growth = 0.5 

# %% extensive functionalities
# explicit targets or exclusions

# enforce target space
deletion_targets = []
addin_targets = [] # has no effect

cofeed_targets = ["EX_tyr__L_e", "EX_trp__L_e", "EX_sucr_e", "EX_tre_e", "EX_adn_e",
                     "EX_3hpppn_e", "EX_3hcinnm_e", "EX_cytd_e", "EX_dad_2_e",
                     "EX_dcyt_e", "EX_dgsn_e", "EX_din_e", "EX_duri_e", "EX_gsn_e",
                     "EX_3gmp_e", "EX_ins_e", "EX_phe__L_e", "EX_lcts_e",
                     "EX_malthx_e", "EX_maltpt_e", "EX_malt_e", "EX_maltttr_e",
                     "EX_malttr_e", "EX_melib_e", "EX_pppn_e", "EX_gthrd_e",
                     "EX_sucr_e", "EX_thm_e", "EX_thymd_e", "EX_tre_e", "EX_xtsn_e",
                     "EX_uri_e", "EX_xtsn_e"] # all reactions specified here are automatically excluded as source targets

source_targets = ["EX_glc__D_e", "EX_glyc_e", "EX_lac__D_e", "EX_ac_e",
                       "EX_glyclt_e", 
                       "EX_gthrd_e", "EX_acald_e", "EX_akg_e", "EX_mal__L_e", 
                       "EX_cit_e", "EX_etoh_e", "EX_fum_e", 
                       "EX_succ_e"] # all reactions specified here are automatically excluded as cofeed targets
mediareduction_targets = []

# explicitly exclude from target space
deletion_exclude_list = []
addin_exclude_list = []
cofeed_exclude_list = ["EX_glc__D_e", "EX_glyc_e", "EX_lac__D_e", "EX_ac_e",
                       "EX_f6p_e", "EX_g6p_e", "EX_glyc3p_e", "EX_glyclt_e", 
                       "EX_gthrd_e", "EX_acald_e", "EX_akg_e", "EX_mal__L_e", 
                       "EX_cit_e", "EX_etoh_e", "EX_fum_e", "EX_for_e", 
                       "EX_gam6p_e", "EX_succ_e"] # exclude explicit (carbon) cofeed targets
source_exclude_list = -1
mediareduction_exclude_list = []


# exclusion lists
subsystem_exclude_list = ['Cell Envelope Biosynthesis',
                          'Exchange',
                          'Inorganic Ion Transport and Metabolism',
                          'Lipopolysaccharide Biosynthesis / Recycling',
                          'Murein Biosynthesis',
                          'Murein Recycling',
                          'Transport, Inner Membrane',
                          'Transport, Outer Membrane',
                          'Transport, Outer Membrane Porin',
                          'tRNA Charging',
                          'TRANSPORT, EXTRACELLULAR',
                          'Transport, extracellular',
                          'Secondary transporters',
                          'Transport, Extracellular'# from recon
                          ]



exchanges_not_to_knockout = ['EX_co2_e', 'EX_h2o_e', 'EX_h_e', 'EX_fe2_e', 'EX_fe3_c',
                     'EX_na1_e', 'EX_cl_e', 'EX_k_e', 'EX_o2_e', 'EX_glc__D_e', 'EX_nh4_e']


# disregard exchange reactions as cofeed and carbon source target
exchanges_not_to_add = ["EX_co2_e", "EX_h2_e", "EX_h2s_e", "EX_o2_e", "EX_ch4_e", "EX_o2s_e",
                         "EX_h2o2_e", "EX_h2o_e", "EX_h_e", "EX_n2o_e", "EX_no_e", "EX_cynt_e",
                        "EX_meoh_e", "EX_cyan_e", "EX_mso3_e", "EX_mepn_e", "EX_tcynt_e"
                         ]

# explicitly include exchange reactions as cofeed and carbon source target
# has no effect if cofeed/source_exclude_list is set to -1 
exchanges_to_keep = []
# maximum carbon content of cofeed or carbon source targets
# has no effect if cofeed/source_exclude_list is set to -1 
num_exchange_carbons_threshold = 8


# %% options for processing deletions targets
# Essentiality of reactions (not considered as deletions targets) is based on the original wild-type model
consider_wildtype_essentiality = True

# Control if reaction with metabolites from different compartments can be deleted
# If allowed, erase transport subsystem from subsystem_exclude_list
allow_transporter_deletion = False

# %% options for processing carbon source and cofeed choice
maximum_source_non_carbon_flux = 10 # [mmol/gDW/h]
maximum_source_non_carbon_mass_flux = 1.80 # [g/gDW/h]
maximum_source_carbon_flux = 60 # [c-mmol/gDW/h]

# specifically define constraints for nitrogen source or cofeed targets, disabled by default
maximum_source_nitrogen_flux = None # [N-mmol/gDW/h]

# cofeed uptake rate is the maximum source uptake rate multiplied by this factor
# a relatively small cofeed rate prevents the cofeed from becoming the sole carbon source
cofeed_rate_factor = 0.25 


user_exclude_list = []




# %% options for processing heterologous reactions (addins)


# assess directions of heterologous reaction by thermodynamic and flux variability analysis
directionality_assessment = True

# if thermodynamic assessment of reactions is done, only consider reactions for which
# Gibbs free energy of reaction is accessable
consider_thermodynamic_valid_reactions_only = True

# currency metabolites do not count as high carbon metabolites
currency_mets = {'h2o', 'co2', 'o2', 'h2o2', 'nh4', 'no2', 'no3', 'no', 'h2s',
                 'so3','so4','h','h2','pi','ppi','coa','accoa','ppcoa','aacoa',
                 'butcoa','succoa','atp','gtp','adp','gdp','amp','gmp','nad',
                 'nadp','nadh','nadph','fad','fadh','na1','ahcys','amet','thf',
                 'mlthf', 'q8h2','q8','mql8','mqn8','2dmmql8','2dmmq8'}


num_carbons_threshold = 35



# %% model processing options
# adapt model bounds and coefficients to avoid numerical issues when solving the MILP
improve_model_numerics = True



