"""
Default config for the reaction direction evaluation
"""

from os.path import dirname
filedir = dirname(__file__) + "\Data"

# %% thermodynamic evaluation
# default values were taken from Noor et al. 2012 (10.1093/bioinformatics/bts317)

# Standadrd Equilibrator parameters
ionic_strength = "0.1 M"
temperature = "298.15 K"
p_h = "7.2"

# reaction direction evaluation
# "RI": reversibility index
# "dG": Gibbs free energy of reaction
reaction_direction_eval_type = "dG" 

# consider (experimental) physiological metabolite data
use_experimental_concentrations = True
metabolite_concentration_filename \
    = filedir + "\conc_Sauer_avg_cofactors.xlsx"
        
# consider worst case szenarios with regard to endogeneous metabolite concentration ranges
evaluate_metabolite_range = True
metabolite_concentration_range = [0.01, 20] # unit mM (Feist 2007)
# metabolite_concentration_range = [0.003, 3] # unit mM (Memote)

# metabolite concentration parameter
standard_physiological_concentration = 1  # unit mM


# reversibility index /dG prime calculation and comparison parameters
ri_threshold = 3 # threshold for the (absolute) reversibility index above which a reaction is irreversible
dG_threshold = 40 # threshold for the (absolute) Gibbs free energy of reaction [kJ/mol] above which a reaction is irreversible
# check_pH_range = False
# pH_limits = [4.8, 9.0]   # minimal and maximal pH in bacteria (Breeuwer et al. 1996 ())
consider_standard_deviation = True


# %% flux variability evaluation

# allowable fraction of deviation from optimal growth
fraction_of_optimum = 0.9

# number of processes to be used
processes = 6

# provide exchange reaction IDs of carbon substrates
filename_carbon_exchange_reactions = filedir + "\carbon_exchange_reaction_list.xlsx"  # use standard list of substrates from Feist et al. 2007