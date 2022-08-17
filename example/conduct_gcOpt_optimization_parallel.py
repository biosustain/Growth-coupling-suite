from growth_coupling_suite.gcOpt_algorithm import gcOpt
from growth_coupling_suite.strain_analysis.strain_design_analysis import StrainDesignAnalyzer
import cobra

# %% Check "gcOpt_config_file.py" for further options
import gcOpt_config_file_parallel as config


if __name__ == "__main__": # necessary when running GCS on parallel workers
    # %% load "mid-scale" metabolic model of E. coli (https://doi.org/10.1038/srep39647)
    model_name = "ECC2.json"
    model = cobra.io.load_json_model("Models/"+model_name)
    
    
    # %% prepare model
    # set glucose uptake rate
    model.exchanges.EX_glc__D_e.lower_bound = -10
    # set oxygen uptake rate
    model.exchanges.EX_o2_e.lower_bound = -20
    # disable co2 uptake
    model.exchanges.EX_co2_e.lower_bound = 0
    
    
    # %% define target reaction -> succinate exchange
    target_reaction = "EX_succ_e"
    
    # %% load heterolgous reaction database model
    hr_database_model = cobra.io.load_json_model("Models/ECC2_hr_database_dir_assessed.json")   
    
    # %% preapre gcOpt config
    config.exchanges_not_to_add.append(target_reaction)
    

    # %% run multiple gcOpt instances in parallel/sequentially
    # Can be useful since the Gurobi solver doesn't make use of all available 
    # threads after running the MILP problem for some time.
    
    # set up config parameters
    # Any config file attribute can be added to the parameter dictionaries 
    parameters_parallel = [
        {"num_total_interventions": 3, "num_deletions": 2, "num_addins":2, "output_file_id": "i3_d2_a2_cf1_cs1"},
        {"num_total_interventions": 6, "num_deletions": 5, "num_addins":5, "output_file_id": "i6_d5_a5_cf1_cs1"},
        {"num_total_interventions": 9, "num_deletions": 8, "num_addins":8, "output_file_id": "i9_d8_a8_cf1_cs1"},
    ]
    # "output_suffix" is added to the "output_file_id"
    parameters_sequential = [
        {"growth_rate_fix": 0.1, "output_suffix": "gr_01"},
        {"growth_rate_fix": 0.3, "output_suffix": "gr_03"},
        {"growth_rate_fix": 0.6, "output_suffix": "gr_06"},
    ]


    # initialize gcOpt class
    GCS = gcOpt.GCOpt(
        model,
        target_reaction,
        hr_database_model=hr_database_model,
        config=config,
        build_gcopt_problem=True
        )
    # run gcOpt in parallel
    GCS.optimize_series(
        parameters_parallel=parameters_parallel,
        parameters_sequential=parameters_sequential,
        max_workers = 3, # number of parallel workers, compare with the number of processes allocated in the config file to avoid overloading your machine
        init_DesignAnalyzer=False
    )

    
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
        determine_significant_designs=True, # duplicate solutions are disrearded,
                                            # design objects (deletions, add-ins, etc.) that do not contribute to the coupling are stripped from solution
        save_results=True,
        save_flux_space_plots=True,
        eval_gpr=True
        )
        

