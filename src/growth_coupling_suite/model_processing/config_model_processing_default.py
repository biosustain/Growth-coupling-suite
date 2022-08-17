"""
Default config file for the model processing routine
"""

# enforce target space
# assign "-1" to exclude list if explcitly specified targets should be used exclusively
deletion_targets = []
addin_targets = [] # has no effect
cofeed_targets = []
source_targets = []
mediareduction_targets = []

# explicitly exclude from target space
# assign "-1" to exclude list if explcitly specified targets should be used exclusively
deletion_exclude_list = []
addin_exclude_list = []
cofeed_exclude_list = ["EX_glc__D_e", "EX_glyc_e", "EX_lac__D_e", "EX_ac_e",
                       "EX_f6p_e", "EX_g6p_e", "EX_glyc3p_e", "EX_glyclt_e", 
                       "EX_gthrd_e", "EX_acald_e", "EX_akg_e", "EX_mal__L_e", 
                       "EX_cit_e", "EX_etoh_e", "EX_fum_e", "EX_for_e", 
                       "EX_gam6p_e", "EX_succ_e"]
source_exclude_list = []
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

currency_mets = {'h2o', 'co2', 'o2', 'h2o2', 'nh4', 'no2', 'no3', 'no', 'h2s',
                 'so3','so4','h','h2','pi','ppi','coa','accoa','ppcoa','aacoa',
                 'butcoa','succoa','atp','gtp','adp','gdp','amp','gmp','nad',
                 'nadp','nadh','nadph','fad','fadh','na1','ahcys','amet','thf',
                 'mlthf', 'q8h2','q8','mql8','mqn8','2dmmql8','2dmmq8'}

exchanges_not_to_knockout = ['EX_co2_e', 'EX_h2o_e', 'EX_h_e', 'EX_fe2_e', 'EX_fe3_c',
                     'EX_na1_e', 'EX_cl_e', 'EX_k_e', 'EX_o2_e', 'EX_glc__D_e', 'EX_nh4_e']


# # disregard metabolite as cofeed target
# cofeed_exclude_list = ["glc__D", "glyc", "lac__D", "ac", "f6p", "g6p", "glyc3p", 
#                        "glyclt", "gthrd", "acald", "akg", "mal__L", "cit",
#                        "etoh", "fum", "for", "gam6p", "succ"]

# # only use these metabolites as potential (carbon) source targets
# source_metabolites_set = ["glc__D", "glyc", "lac__D", "ac", "f6p", "g6p", "glyc3p", 
#                        "glyclt", "gthrd", "acald", "akg", "mal__L", "cit",
#                        "etoh", "fum", "for", "gam6p", "succ"]

# disregard exchange reactions as cofeed and carbon source target
exchanges_not_to_add = ["EX_co2_e", "EX_h2_e", "EX_h2s_e", "EX_o2_e", "EX_ch4_e", "EX_o2s_e",
                         "EX_h2o2_e", "EX_h2o_e", "EX_h_e", "EX_n2o_e", "EX_no_e", "EX_cynt_e",
                        "EX_meoh_e", "EX_cyan_e", "EX_mso3_e", "EX_mepn_e", "EX_tcynt_e",
                         "DM_amob_c", "DM_5drib_c", "DM_oxam_c", "DM_aacald_c", "DM_4crsol_c",
                         "DM_mththf_c"]
# explicitly include exchange reactions as cofeed and carbon source target
exchanges_to_keep = ["EX_tyr__L_e", "EX_trp__L_e", "EX_sucr_e", "EX_tre_e", "EX_adn_e",
                     "EX_3hpppn_e", "EX_3hcinnm_e", "EX_cytd_e", "EX_dad_2_e",
                     "EX_dcyt_e", "EX_dgsn_e", "EX_din_e", "EX_duri_e", "EX_gsn_e",
                     "EX_3gmp_e", "EX_ins_e", "EX_phe__L_e", "EX_lcts_e",
                     "EX_malthx_e", "EX_maltpt_e", "EX_malt_e", "EX_maltttr_e",
                     "EX_malttr_e", "EX_melib_e", "EX_pppn_e", "EX_gthrd_e",
                     "EX_sucr_e", "EX_thm_e", "EX_thymd_e", "EX_tre_e", "EX_xtsn_e",
                     "EX_uri_e", "EX_xtsn_e", "EX_for_e", "EX_urea_e"]
# maximum carbon content of cofeed or carbon source targets
num_exchange_carbons_threshold = 8

# define constraints for carbon/non-carbon source  targets
maximum_source_non_carbon_flux = 10 # [mmol/gDW/h]
maximum_source_non_carbon_mass_flux = 1.80 # [g/gDW/h]
maximum_source_carbon_flux = 60 # [C-mmol/gDW/h]

num_carbons_threshold = 35

# cofeed uptake rate is the maximum source uptake rate multiplied by this factor
# a relatively small cofeed rate prevents the cofeed from becoming the sole carbon source
cofeed_rate_factor = 0.25 


# specifically define constraints for nitrogen source or cofeed targets, disabled by default
maximum_source_nitrogen_flux = None # [N-mmol/gDW/h]

user_exclude_list = []

# deletion target space

# Control if reaction with metabolites from different compartments can be deleted
# If allowed, erase transport subsystem from subsystem_exclude_list
allow_transporter_deletion = False


