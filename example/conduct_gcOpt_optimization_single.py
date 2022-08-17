# %%
from growth_coupling_suite.gcOpt_algorithm import gcOpt
from growth_coupling_suite.strain_analysis.strain_design_analysis import StrainDesignAnalyzer
import cobra

# %% Check "gcOpt_config_file.py" for further options
import gcOpt_config_file_single as config

# %% load "mid-scale" metabolic model of E. coli (https://doi.org/10.1038/srep39647)
model_name = "ECC2.json"
model = cobra.io.load_json_model("Models/" + model_name)

# %% load heterolgous reaction database model
hr_database_model = cobra.io.load_json_model("Models/ECC2_hr_database_dir_assessed.json")

# %% prepare model
# set glucose uptake rate
model.exchanges.EX_glc__D_e.lower_bound = -10
# set oxygen uptake rate
model.exchanges.EX_o2_e.lower_bound = -20
# disable co2 uptake
model.exchanges.EX_co2_e.lower_bound = 0


# %% define target reaction -> succinate exchange
target_reaction = "EX_succ_e"


# %% preapre gcOpt config
config.exchanges_not_to_add.append(target_reaction)


# %% run a single gcOpt instance
# load gcOpt class
GCS = gcOpt.GCOpt(model, target_reaction, config=config, build_gcopt_problem=True)
# solve gcOpt MILP problem
GCS.optimize(init_DesignAnalyzer=False)


# %% analyze, save, and plot results
# load Strain Design Analyzer
sda = StrainDesignAnalyzer()
# load results files, all files in the specified folder will be loaded
sda.load_strain_design_files_from_dir(
    config.output_dir,
    eval_gpr=True
    )

# get a summary of computed growth-coupled designs 
sda.growth_coupling_summary(
    results_filename="gcOpt_summary_result",
    results_dir=config.output_dir,
    determine_significant_designs=True, # duplicate solutions are disregarded,
                                        # design objects (deletions, add-ins, etc.) that do not contribute to the coupling are stripped from solution
    save_results=True,
    save_flux_space_plots=True,
    eval_gpr=True
    )
    

