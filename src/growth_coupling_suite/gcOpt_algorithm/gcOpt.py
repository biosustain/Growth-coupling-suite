"""
Initializes the setting up and solving of a gcOpt problem for determining
growth-coupled strain designs:
gcOpt principle: Maximize a minimal guaranteed target reaction rate at a fix growth rate
"""

from growth_coupling_suite.bilevel_optimization_algorithm.bilevel_optimization \
    import BilevelModel
from growth_coupling_suite.gcOpt_algorithm import config_gcOpt_default as config_default
from growth_coupling_suite.util.util import check_config, is_reaction_in_model
from growth_coupling_suite.model_processing import model_processing, reduce_model
from growth_coupling_suite.strain_analysis.strain_design_analysis import StrainDesignAnalyzer
from growth_coupling_suite.util.monitor import Monitor

import cobra

from types import SimpleNamespace
from copy import deepcopy
import pickle
from os import mkdir, listdir
from os.path import isdir, isfile
import shutil
from warnings import warn
from pathlib import Path
from time import sleep

from concurrent import futures
from multiprocessing import cpu_count


# Disable gurobi logging output
try:
    import gurobipy
    gurobipy.setParam("OutputFlag", 0)
except ImportError:
    pass


class GCOpt():
    
    
    # gcopt class for setting up and handling of growth-coupling_problems 
    def __init__(self, model, target_rxn,
                 hr_database_model=None,
                 config=config_default,
                 build_gcopt_problem=True):
        """
        Initializes a gcOpt object

        Parameters
        ----------
        model : cobra.Model
            Metabolic model in cobra format.
        target_rxn : str
            ID of the growth-coupling target reaction in the model.
        hr_database_model : cobra.Model
            Heterologous reaction database model in cobra format.    
        config : TYPE, optional
            configuration module including optimization and model parameters. The default is config_default.
        build_gcopt_problem : TYPE, optional
            Build the bilevel optimization problem from the model and parameters. The default is True.

        Returns
        -------
        None.

        """
        
        # save config as SimpleNamespace with attribute access
        self.config = config
        
        # monitor changes in config module
        self._monitor = Monitor()
        self._monitor.add_object("config", self.config)
                
        # load and save inputs
        self.model = model
        self._model_original = model.copy()
        self.target_rxn = target_rxn
        self.hr_database_model = hr_database_model
        
        # create necessary directories for saving results
        if not(isdir(config.output_dir)): mkdir(config.output_dir)
        self.callback_output_dir = config.output_dir + "/callback_solutions"
        if not(isdir(self.callback_output_dir)): mkdir(self.callback_output_dir)
        
        # initialize gcOpt problem  
        if build_gcopt_problem:
            self._build_gcopt_problem()
            self._rebuild_bilevel_model = False
        else:
            self._rebuild_bilevel_model = True
        
        # initialize callback solutions
        self._callback_solutions = {}
        self._status = []
        
        

    
    # define properties
    @property
    def model(self):
        return self._model
    
    @model.setter
    def model(self, model):           
        self._model = model
       
    @property
    def target_rxn(self):
        return self._target_rxn    
       
    @target_rxn.setter
    def target_rxn(self, value):
        if is_reaction_in_model(self.model, value):
            self._target_rxn = value
            # adapt bilevel model
            self._rebuild_bilevel_model = True
        else:
            warn("Target reaction " + value + " is not in model", UserWarning)
    
    @property
    def callback_output_dir(self):
        return self._callback_output_dir
    
    @callback_output_dir.setter
    def callback_output_dir(self, dirname):
        self._callback_output_dir = dirname
    
    @property
    def config(self):
        return self._config
  
    @config.setter
    def config(self, module):
        # check config module with all default config files in the GC suite
        check_config(module)
        # translate module into SimpleNamespace
        # TODO: generally exclude modules
        disregard_attr = ["os", "json"]
        attributes = dir(module)
        self._config = SimpleNamespace()
        for attribute in attributes:
            if not(attribute[0] == "_") and not(attribute[:1] == "__") \
                and not(attribute in disregard_attr):
                setattr(self._config, attribute, getattr(module, attribute))



       
    def optimize(self, save_solution=True, init_DesignAnalyzer=True):
        """
        Solve and optimize the bilevel model

        Parameters
        ----------
        save_solution : bool, optional
            save solution in seperate files. The default is True.
        init_DesignAnalyzer : TYPE, optional
            Create a StrainDesigner object and load solutions into it. The default is True.

        Returns
        -------
        None.

        """

        
        # check configuration has changed
        if self._rebuild_bilevel_model or self._monitor.is_changed("config", self.config):
            print("Config and parameters changed. Rebuild gcOpt problem...")
            # load original model
            self.reset_model()
            # self.model = self._model_original.copy()
            # config changes, first rebuild bilevel model
            self._build_gcopt_problem()
            # reset rebuild problem parameter
            self._rebuild_bilevel_model = False
         
        # move and save old solution files
        self._move_previous_solution_files()
         
        # optimize model 
        self._callback_solutions, self._status = self._bilevel_model.optimize(
            base_model=self._base_model.copy()
            )
        
        # save solution
        if save_solution:
            self.save_minimal_solution()
        
        # create Design Analyzer class
        if init_DesignAnalyzer:
            self._init_design_analysis()

            
         
    def optimize_series(self, parameters_parallel=[], parameters_sequential=[], max_workers=3,
                        save_solution=True, init_DesignAnalyzer=False):
        """
        Solving of multiple gcOpt problems in paralell and series.
        Sequential parameter sets are run in each parallel instance

        Parameters
        ----------
        parameters_parallel : list, optional
            list of dicts with parameters for each single gcOpt instance run in parallel. The default is [].
        parameters_sequential : list, optional
            list of dicts with parameters for each single gcOpt instance run in series. The default is [].
        max_workers : int, optional
            Number of workers allocated to each gcOpt instance run in parallel. The default is 3.
        save_solution : bool, optional
            Save solutions in separate files. The default is True.
        init_DesignAnalyzer : bool, optional
            Create a StrainDesigner object and load solutions into it. The default is False.

        Raises
        ------
        Exception
            Missing parameters.

        Returns
        -------
        None.

        """
        
        # check if individual file IDs are provided for parameter sets
        for parameters in parameters_parallel:
            if not any(key in parameters for key in ["output_file_id", "output_suffix"]):
                raise Exception("No unique file ID for each parallel parameter set provided!")
        
        for parameters in parameters_sequential:
            if not any(key in parameters for key in ["output_file_id", "output_suffix"]):
                raise Exception("No unique file ID for each sequential parameter set provided!")
        
        # check number of processes
        if max_workers*self.config.processes >= cpu_count():
            warn("Too many parallel processes allocated. CPU may be overloaded.", UserWarning)

        # Spwan a process for each parallel parameter set
        with futures.ProcessPoolExecutor(max_workers=max_workers) as process:
            for parameters in parameters_parallel:
                # submit process
                process.submit(self._call_optimize_series,
                               parameters,
                               parameters_sequential,
                               save_solution,
                               init_DesignAnalyzer)
                
                # delay next process submission (sec)
                sleep(20)


    def save_minimal_solution(self):
        """
        Save models and design solutions 

        Returns
        -------
        None.

        """
        # construct and save solution for later use
        filename = "gcOpt_solution_dict_{0}{1}.pickle"
        if len(self.config.output_file_id) > 0:
            filename = filename.format(self.target_rxn, "_" + self.config.output_file_id)
        else:
            filename = filename.format(self.target_rxn, "")
        
        # check solution
        sol = {}
        if not(self._callback_solutions):
            print("No solution available")
        else:   
            # interventions   
            sol["solutions_dict"] =   self._callback_solutions  
            # heterologous reaction database model
            sol["heterologous_reaction_database_model"] = self.hr_database_model
            # model ID
            sol["model_id"] = self.model.id
            # model
            sol["model"] = self._model_original
            # target reaction
            sol["target_reaction"] = self.target_rxn
            # config file
            sol["config"] = self.config
            # mutables
            sol["mutables"] = {"deletion_targets": self.deletion_targets,
                               "addin_targets": self.addin_targets,
                               "source_targets": self.source_targets,
                               "cofeed_targets": self.cofeed_targets,
                               "mediareduction_targets": self.mediareduction_targets}
        
        # save solution dictionary       
        with open(self.config.output_dir + "/" + filename, "wb") as f:
            pickle.dump(sol, f)
            
   
            
    # function for sequentially calling optimize
    def _call_optimize_series(self, parameters, parameters_sequential, 
                             save_solution=True, init_DesignAnalyzer=True):
        """
        Call optimize with a series of parameters. A gcOpt problem is built
        for each sequential parameter set
        

        Parameters
        ----------
        parameters : dict
            Parameter set for this specific gcOpt instance.
        parameters_sequential : list
            list of parameter sets for gcOpt instances run in series.
        save_solution : bool, optional
            Save solutions in separate files. The default is True.
        init_DesignAnalyzer : bool, optional
            Create a StrainDesigner object and load solutions into it. The default is True.

        Returns
        -------
        None.

        """
        # apply global parameters
        self._change_config_parameters(self.config, parameters)
    
        # save base output file id
        output_file_id_base = self.config.output_file_id
        # loop through sequential parameters if provided
        if len(parameters_sequential) == 0:
            # execute gcOpt
            self._rebuild_bilevel_model = True
            self.optimize(save_solution=save_solution,
                              init_DesignAnalyzer=init_DesignAnalyzer)
           
        else:
            for p_seq in parameters_sequential:
    
                # use base output file id
                self.config.output_file_id = output_file_id_base
                # change sequential parameters
                self._change_config_parameters(self.config, p_seq)
                                    
                # execute gcOpt
                self._rebuild_bilevel_model = True
                self.optimize(save_solution=save_solution,
                              init_DesignAnalyzer=init_DesignAnalyzer)

    # function for changing config parameters
    def _change_config_parameters(self, config, parameter_dict): 
        """
        Change config parameters

        Parameters
        ----------
        config : SimpleNamespace
            config dict.
        parameter_dict : dict
            parameters (keys) to be changed in config.

        Returns
        -------
        None.

        """
        
        
        suffix = ""
        for parameter, value in parameter_dict.items():
            if parameter == "output_suffix":
                suffix = value
            else:
                setattr(config, parameter, value)
        # set suffix
        if len(suffix) > 0:
            config.output_file_id += "_" + suffix


    def _build_gcopt_problem(self):
        """
        Build a bilevel gcOpt model from the metabolic and the optimization parameters

        Returns
        -------
        None.

        """
        
        # process target reactions
        self._process_mutable_reactions()
        
        # create and save base model
        self._create_base_model()
        
        # reduce model and target space
        self._reduce_model_complexity()
        
        # prepare model
        self._prepare_model()
        
        # create bilevel model
        self._create_bilevel_model()
        
        # reset config monitor
        self._monitor.add_object("config", self.config) 






    def _init_design_analysis(self):
        """
        create StrainDesignAnalyzer object and parse solutions

        Returns
        -------
        None.

        """
        
        self.DesignAnalyzer = StrainDesignAnalyzer(self._model_original,
                                design_solutions_dict=self._callback_solutions,
                                heterologous_reaction_database_model=self.hr_database_model,
                                  )   

    def _move_previous_solution_files(self):
        """
        move files of previously computed solutions to a separate folder

        Returns
        -------
        None.

        """

        previous_files = []
        # check if old callback results exist
        cb_dir = Path(self.callback_output_dir)
        files = listdir(cb_dir) 
        for file in files:
            if isfile(cb_dir.joinpath(file)) \
                and (self.target_rxn in file) \
                    and (self.config.output_file_id in file)\
                        and file.startswith("callback_"):
                previous_files.append(cb_dir.joinpath(file))
        
        
        # check if old gcOpt results files exist
        op_dir = Path(self.config.output_dir)
        files = listdir(op_dir)
        for file in files:
            if isfile(op_dir.joinpath(file)) \
                and (self.target_rxn in file) \
                    and (self.config.output_file_id in file)\
                        and file.startswith("gcOpt_solution_"):
                previous_files.append(op_dir.joinpath(file))
                                            
        # move previous files
        if len(previous_files) > 0:
            # store old solution in extra folder
            previous_results_path_format = self.config.output_dir + "/previous_results_{}/"
            for i in range(1000):
                if not isdir(previous_results_path_format.format(i+1)):
                    previous_results_path = Path(previous_results_path_format.format(i+1))
                    mkdir(previous_results_path)
                    break
                    
            # move files
            for file in previous_files:
                shutil.move(file, previous_results_path.joinpath(file.name))
            


        
    
        
    def _prepare_model(self):
        """
        Prepare COBRA model for gcOpt optimization
        - check validity of the passed biomass equation
        - improve the model numerics by
            - restricting flux bounds
            - removal of biomass components with low stoichiometric coefficients
            - factorizing stoichiometric coefficients of biomass components
    
        
        Returns
        -------
        None.
    
        """
        
        # load parameter
        growth_rate_fix = self.config.growth_rate_fix
        
        
        # check biomass formation reaction id
        if not(is_reaction_in_model(self.model, self.config.biomass_id)):
            # extract biomass equation from objective function
            objective_reactions = model_processing.get_objective_reaction(self.model)
            if len(objective_reactions) == 1:
                self.config.biomass_id = objective_reactions[0]
            else:
                raise Exception("No or invalid Biomass equation ID provided")
        

        # save medium composition
        self._medium_specs = {}
        for ex in self.model.reactions:
            if "EX_" in ex.id and ex.lower_bound < 0:
                # uptake reaction
                self._medium_specs[ex.id] = {"exchange_reaction_id": ex.id,
                                            "lower_bound": ex.lower_bound}
      
        
        # adapt bounds and coefficients to improve numerics
        if self.config.improve_model_numerics:
            # restrict lower/upper bounds to -100/100
            for r in self.model.reactions:
                if r.upper_bound > 100: r.upper_bound = 100
                if r.lower_bound < -100: r.lower_bound = -100
                
            # delete biomass precursors with very low coeffcients
            biomass_rxn = self.model.reactions.get_by_id(self.config.biomass_id)
            remove_met = []
            print("Remove metabolites from biomass equation:")
            for m, c in biomass_rxn.metabolites.items():
                if abs(c) < 1e-04:
                    remove_met.append(m)
                    print("\t", m.id, " coefficient: ", str(c))
                    
            if len(remove_met) > 0:
                for m in remove_met:
                    biomass_rxn.add_metabolites({m: 0}, combine=False)
            else:
                print("\t None")
                
            # factorize biomass equation coefficients and bounds
            factor = 10
            for m, c in biomass_rxn.metabolites.items():
                biomass_rxn.add_metabolites({m: c * factor}, combine=False)
                
            growth_rate_fix = growth_rate_fix / factor
                
        # fix biomass formation rate
        bm_rxn = self.model.reactions.get_by_id(self.config.biomass_id)
        bm_rxn.lower_bound = growth_rate_fix
        bm_rxn.upper_bound = growth_rate_fix
        # add biomass enforcement constraint
        bio_fix = cobra.Metabolite("Biomass_fix")
        bm_rxn.add_metabolites({bio_fix: 1})
        bio_fix.constraint.ub = growth_rate_fix
        bio_fix.constraint.lb = growth_rate_fix
        
        
    def _process_mutable_reactions(self):
        """
        Define reactions which can be de-/activated for optimizing growth-coupling
        
        The resulting model contains all heterologous reactions and has original
        exchange reactions for substrate uptake replaced by source reactions
        for the respective metabolite
        

        Returns
        -------
        None.

        """
        # get deletions, addin targets etc
        
        # load config 
        config = deepcopy(self.config)
        
        # analyze number of mutable reactions
        if (config.num_deletions + config.num_addins + config.num_mediareductions\
            + config.num_cofeeds + config.num_carbonsources) <= 0:
            # at least consider deletion targets
            self.config.num_deletions = self.config_default.num_deletions
        
        # put target reaction on exclude list
        config.user_exclude_list.append(self.target_rxn)
               
        
        # get deletion targets
        if config.num_deletions > 0:
            self.deletion_targets \
                = model_processing.get_deletion_target_space(self.model,
                                                             config=config)
                
        else:
            self.deletion_targets = []
        print('\t', len(self.deletion_targets), ' deletion targets')
            
        # get addin targets and integrate in the model (heterologous reactions)
        if config.num_addins > 0:
            self.model, self.addin_targets, self.hr_database_model \
                = model_processing.add_heterologous_reactions(self.model,
                                                              hr_database_model=self.hr_database_model,
                                                              config=config)
        else:
            self.addin_targets = []
        print('\t', len(self.addin_targets), ' addin targets')
                    
        # get mediareduction targets
        if config.num_mediareductions > 0:
            # get mediareduction targets
            pass
        else:
            self.mediareduction_targets = []
            
            
        # get cofeed and carbon source targets (source targets)
        if config.num_cofeeds > 0 or config.num_carbonsources > 0:
            # get and setcofeed targets
            # include all exchange reactions which are not uptake reactions atm
            self.model, self.source_targets, self.cofeed_targets\
                = model_processing.set_source_targets(self.model,
                                                      allow_carbon_source_switch=config.num_carbonsources > 0,
                                                      allow_non_carbon_sources=config.num_cofeeds > 0,
                                                      allow_cofeed_targets=config.num_cofeeds > 0,
                                                      config=config)
        else:
            self.cofeed_targets = [] 
            self.source_targets = []
            
        print('\t', len(self.cofeed_targets), ' cofeed targets')
        print('\t', len(self.source_targets), ' source targets')
            

   

    def _reduce_model_complexity(self):
        """
        Minimize the target space by removing variables (target and model reactions)

        Returns
        -------
        None.

        """
        # reduce model variables and target space
        
        # remove blocked reactions
        print("Remove blocked reactions...")
        self.model, blocked_reactions \
            = reduce_model.blocked_reactions(self.model, remove_blocked_reactions=True)
        print("\t", len(blocked_reactions), ' reactions blocked')

        # remove blocked reactions as target variables
        self.deletion_targets \
            = list(set(self.deletion_targets).difference(set(blocked_reactions)))
        self.addin_targets \
            = list(set(self.addin_targets).difference(set(blocked_reactions)))
        self.source_targets \
            = list(set(self.source_targets).difference(set(blocked_reactions)))
        self.cofeed_targets \
            = list(set(self.cofeed_targets).difference(set(blocked_reactions)))
        self.mediareduction_targets \
            = list(set(self.mediareduction_targets).difference(set(blocked_reactions)))
        self.blocked_reactions = blocked_reactions

        # identify essential reactions not be knocked out
        print("Identify essential reactions...")
        if self.config.consider_wildtype_essentiality:
            # use original wild-type model without heterologous reactions, alternative carbon sources etc
            essential_reactions = reduce_model.essential_reactions(self._model_original,
                                               reaction_list=self.deletion_targets,
                                               objective_cutoff=self.config.growth_rate_fix)            
        else:
            # use model including heterologous reactions
            essential_reactions = reduce_model.essential_reactions(self.model,
                                               reaction_list=self.deletion_targets,
                                               objective_cutoff=self.config.growth_rate_fix)
        self.deletion_targets \
            = list(set(self.deletion_targets).difference(set(essential_reactions)))
        self.essential_reactions = essential_reactions
        print("\t", len(essential_reactions), ' reactions essential')

        
    def _create_bilevel_model(self):
        """
        Create a gcOpt bilevel model

        Returns
        -------
        None.

        """
        
        # load config
        config = deepcopy(self.config)
        
        # setup bilevel options
        bilevel_options = {"time_limit": config.time_limit,
               "output_dir": self.callback_output_dir,
               "output_file_id": config.output_file_id,
               "medium": self._medium_specs,
               "eval_gpr": self.config.eval_gpr,
               "eval_max_growth": config.eval_max_growth,
               "n_total_interventions": config.num_total_interventions,
               "n_deletion": config.num_deletions,
               "n_addin": config.num_addins,
               "n_cofeed": config.num_cofeeds,
               "n_source": config.num_carbonsources,
               "n_mediareduction": config.num_mediareductions,
               "n_threads": config.processes,
               "node_memory_limit": config.node_memory_limit}
        
        # set up gcOpt problem model
        self._bilevel_model = BilevelModel(self.model.copy(),
                                  self.target_rxn, "min",
                                  self.target_rxn, "max",
                                  deletion_targets=self.deletion_targets,
                                  addin_targets=self.addin_targets,
                                  source_targets=self.source_targets,
                                  cofeed_targets=self.cofeed_targets,
                                  mediareduction_targets=self.mediareduction_targets,
                                  options=bilevel_options)
        
        
    def _create_base_model(self):
        """
        Create a base model from the current model
        
        The base model is primarily used to assess the feasibility of a strain design
        in the original COBRA format of the model and to compute design properties
        like the maximum growth rate
        

        Returns
        -------
        None.

        """
        
        
        
        # copy from processed model
        self._base_model = self.model.copy()
        
        # if available, block source and cofeed reactions
        # only source/cofeed exchanges from specified strain designs will be allowed 
        for r in self._base_model.reactions:
            if r.id.startswith('source_') or r.id.startswith('cofeed_'):
                r.lower_bound = 0
                
        # if available, block heterologous reactions
        # heterologous reactions will only be activated if a design demands it
        for r in self._base_model.reactions:
            if r.id.endswith('__hr'):
                r.lower_bound = 0
                r.upper_bound = 0
        
        


    def reset_model(self):
        """
        reload originally provided model and replace model

        Returns
        -------
        None.

        """
        
        self.model = self._model_original.copy()