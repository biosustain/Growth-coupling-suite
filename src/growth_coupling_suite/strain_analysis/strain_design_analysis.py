# Analyze strain designs
# inherit from StrainAnalyzer class

from growth_coupling_suite.strain_analysis.strain_analysis import StrainAnalyzer
from growth_coupling_suite.strain_analysis.strain_design_solution import Solution
from growth_coupling_suite.util import util
from growth_coupling_suite.strain_analysis.design_variable_significance import design_significance
from growth_coupling_suite.strain_analysis import add_prebuilt_pathways

import cobra

import pandas as pd

from types import SimpleNamespace
from os.path import isfile
from os import listdir, getcwd, mkdir
import pickle
from itertools import combinations
from warnings import warn
from pathlib import Path
from copy import deepcopy
from matplotlib import pyplot as plt
import numpy as np

from typing import TYPE_CHECKING, Optional, Dict, List


class StrainDesignAnalyzer(StrainAnalyzer):
    def __init__(self, model=None,
                 filename="",
                 design={},
                 design_solutions_dict={},
                 config=SimpleNamespace(),
                 target_reaction=None,
                 heterologous_reaction_database_model=None,
                 copy_model=True, save_model=True):
        
        # initialize parameters
        self.design_key = "design_{0}"
        self._design = {}
        self._applied_design = None
        self._reverse_design = {}
        self.solutions = SimpleNamespace()
        self.target_reaction = target_reaction
        
    
        # load/initialize model
        if model:
            # call constructor of StrainAnalyzer class
            self._init_model(model, copy_model, save_model, target_reaction)
            
            # load provided designs
            self.heterologous_reaction_database_model = heterologous_reaction_database_model
            self.load_strain_design(design)
            # initialize solutions dictionary       
            self._parse_solutions(design_solutions_dict, config=config)
            # save wild-type medium
            self.save_current_medium()
        
               
        # save other inputs
        self._copy_model = copy_model
        self._save_model = save_model
        
        
        # load strain design from file
        # load model if not provided before       
        self.load_strain_design_file(filename,
                                    only_load_design=False,
                                    verbose=False,
                                    load_model=not(hasattr(self, "model")),
                                    target_reaction=target_reaction,
                                    )
        
        


        
    def save_strain_design_file(self, filename):
        # save current strain design solutions
        
        
        # design solutions available?
        if len(self.solutions.__dict__) == 0:
            warn("No strain design solutions available for saving", UserWarning)
            return
        
        design_dict = {}
        
        # reverse currently loaded strain design
        self.revert_strain_design(verbose=False)
        # save model
        design_dict["model"] = self.model
        design_dict["model_id"] = self.model.id
        
        # save heterologous reaction database
        if hasattr(self, "heterologous_reaction_database_model"):
            if not(self.heterologous_reaction_database_model):
                design_dict["heterologous_reaction_database_model"] = []
            else:
                design_dict["heterologous_reaction_database_model"] = self.heterologous_reaction_database_model
                
        # save target reaction ID
        design_dict["target_reaction"] = self.target_reaction
        
        # save solution dict and config
        design_dict["solutions_dict"] = {}
        design_dict["config"] = {}
        for key, sol in self.solutions.__dict__.items():
            design_dict["solutions_dict"][sol._solution_key] = sol._design_solution
            design_dict["config"][sol._solution_key] = sol._config
            

        
        # save design dict
        with open(filename, "wb") as f:
            pickle.dump(design_dict, f)
            
        

            
    def load_strain_design_file(self, filename, only_load_design=False, reduce_design=False,
                                eval_gpr=False, default_design={}, verbose=True,**kwargs):
                                
                                # target_reaction=None, heterologous_reaction_database_model=None,
                                # load_model=True, copy_model=True, save_model=True, verbose=True):
        
        loading_complete = True                           
        
        # initialize path to file
        self.filename = Path(filename)     
        if isfile(self.filename):
            # load results file
            with open(self.filename, "rb") as f:
                res = pickle.load(f)
                
        else:
            if verbose:
                print("No results or strain desgin file provided")
            return
        
        print("Load file: " + filename.name)

        if not(only_load_design):

            # load target reaction
            if "target_reaction" in res:
                self.target_reaction = res["target_reaction"]    
            elif verbose:
                warn("No target reaction provided. Abort loading of strain design")

            
            # load and initialize model
            if "model" in res:
                self._init_model(res["model"],
                                 copy_model=self._copy_model,
                                 save_model=self._save_model,
                                 target_reaction=self.target_reaction)
            else:
                loading_complete = False
                if verbose: warn("No model provided in strain design file", UserWarning)
        
            # load heterologous reaction database model
            if "heterologous_reaction_database_model" in res:
                self.heterologous_reaction_database_model \
                    = res["heterologous_reaction_database_model"]
                    
            elif verbose:
                warn("No heterologous reaction database model provided", UserWarning)
                            
        
        
        # load user parameters
        for key, value in kwargs.items():
            setattr(self, key, value)


        # load strain design solutions and config
        if "solutions_dict" in res:
            self._parse_solutions(res["solutions_dict"],
                                  config=res["config"],
                                  default_design=default_design,
                                  filename=self.filename,
                                  check_duplicate=True,
                                  eval_gpr=eval_gpr,
                                  reduce_design=reduce_design
                                  )
            # save wild-type medium
            self.save_current_medium()
            
        else:
            loading_complete = False
            warn("Provide strain design solutions in solutions_dict dictionary")
            
        return loading_complete
        
        
        
         
    def load_strain_design_files_from_dir(self, dirname, common_name=".pickle",
                                          default_design={},
                                          only_load_designs=False, reduce_design=False,
                                          eval_gpr=False):
        """
        Load all available strain design results files in the specified directory
        If model, target_reaction, heterologous reaction database model not provided load from first file

        Parameters
        ----------
        dirname : str
            path to results directory.

        Returns
        -------
        None.

        """

     
        dirname = Path(dirname)
        
        # parameters
        file_ending = ".pickle" # filetype must be .pickle
        filename_generic = 'gcOpt_solution_dict' # typical strain design file 
        
        # check if directory exists
        if not(dirname.is_dir()):
            warn("Path to directory does not exist: " + str(dirname))
            return
        
        # load all design files in directory
        filenames = listdir(dirname)
        first_file = True
        for f in filenames:
            # correct file type?
            if not(f.endswith(file_ending)) or not(f.startswith(filename_generic)):
                continue
            # filename includes common name?
            if common_name not in f:
                continue
            
            # create path to file
            path2file = dirname.joinpath(f)
            
            if first_file:
                # load model, target reaction etc and other information from file
                loading_complete = self.load_strain_design_file(path2file, only_load_design=False, 
                                                                default_design=default_design,
                                                                reduce_design=reduce_design, verbose=True,
                                                                eval_gpr=eval_gpr)
                if loading_complete: first_file=False
                
            else:
                # load strain design
                loading_complete = self.load_strain_design_file(path2file, only_load_design=True,
                                                                default_design=default_design,
                                                                reduce_design=reduce_design, verbose=False,
                                                                eval_gpr=eval_gpr)
            
        
            
    def load_strain_design(self, design={}, name="", medium={}, verbose=True) -> None:
    
        """
        load a strain design from a dict of interventions and medium composition
        
        :param dict design: interventions of strain design
        :param str name: name of the strain design
        :param dict medium: bounds of exchange reactions defining the extracellular medium
        :param bool verbose: allow printing of messages
        
        :return None
        
        """

        
        # helper functions
        def _create_reaction(rxn_id, rxn_name, lower_bound, upper_bound,
                                 mets_dict):                
                # add reaction to model, create from scratch to avoid name assignment problems
                rxn_to_add = cobra.Reaction(rxn_id,
                                            name=rxn_name,
                                            lower_bound=lower_bound,
                                            upper_bound=upper_bound)

                # add metabolites
                rxn_to_add.add_metabolites(mets_dict)
                return rxn_to_add
        

        # reverse current strain design
        self.revert_strain_design(verbose=verbose)
                    
            
        if verbose: print("Apply parsed strain design (" + name + ")...")    
            
        # apply medium
        old_medium = {}
        if len(medium) > 0:
            if verbose: print("\tSet medium composition")
            # reset current medium
            
            # save old medium bounds
            for ex in self.model.boundary:
                if "EX_" in ex.id and ex.lower_bound < 0:
                    old_medium[ex.id] = {"exchange_reaction_id": ex.id,
                                          "lower_bound": ex.lower_bound}
                    ex.lower_bound = 0
            # set new medium composition
            for key, medium_component in medium.items():
                # get reaction
                ex_rxn = self.model.reactions.get_by_id(medium_component["exchange_reaction_id"])
                # save old bound
                old_medium[key] = {
                    "exchange_reaction_id": medium_component["exchange_reaction_id"],
                    "lower_bound": ex_rxn.lower_bound}
                # set new bound
                ex_rxn.lower_bound = medium_component["lower_bound"]
                
        else:
            if verbose: print("\tNo medium composition provided")
        
        # save medium
        self._reverse_design["medium"] = old_medium   
            
            
            
            
        # load new strain design
        if len(design) > 0:
            if verbose: print("\tSet design interventions")  

            for pt_name, pt in design.items():
                # process deletion, cofeed, mediareduction targets
                if pt["type"] == "deletion" \
                            or pt["type"] == "mediareduction":
                    # get reaction
                    pt_rxn = self.model.reactions.get_by_id(pt["ID"])
                    # save current reaction bounds
                    self._reverse_design[pt_name] = {"ID": pt["ID"],
                                                     "type": pt["type"],
                                                     "lower_bound": pt_rxn.lower_bound,
                                                     "upper_bound": pt_rxn.upper_bound}
                    # apply perturbation
                    if "lower_bound" in pt:
                        pt_rxn.lower_bound = pt["lower_bound"]
                    else:
                        pt_rxn.lower_bound = 0
                    if "upper_bound" in pt:
                        pt_rxn.upper_bound = pt["upper_bound"]
                    else:
                        pt_rxn.upper_bound = 0   
    
                    
                elif pt["type"] == "source" or pt["type"] == "cofeed":
                    # process cofeed target, add respective source reaction
                    
                    # extract cofed metabolite
                    if "source_" in pt["ID"]:
                        met_id = pt["ID"].replace("source_", "")
                    elif "cofeed_" in pt["ID"]:
                        met_id = pt["ID"].replace("cofeed_", "")
                    else:
                        met_id = pt["ID"]
                        
                    # delete source compartment id
                    met_id = met_id.replace("_src", "")
                        
                    # find metabolite in the model
                    try:
                        met = self.model.metabolites.get_by_id(met_id)
                    except:
                        # warn("Metabolite not found in model: " + met_id, UserWarning)
                        raise Exception("Metabolite not found in model: " + met_id)
                    
                    
                    # add source reaction to model
                    # pt_id = "source_" + met_id 
                    pt_id = pt["ID"]
                    source_reaction = cobra.Reaction(pt_id)
                    source_reaction.add_metabolites({met: -1})
                    source_reaction.bounds = (pt["lower_bound"], pt["upper_bound"])
                    self.model.add_reaction(source_reaction)
                            
                    # save to reverse cofeed
                    self._reverse_design[pt_name] = {"ID": pt_id,
                                                     "type": pt["type"]}
                  
                # process addin targets    
                elif pt["type"] == "addin":
                
                    # check for heterologous reaction tag
                    if pt["ID"][-4:] == "__hr":
                        pt_id = pt["ID"][:-4]
                    elif pt["ID"][-3:] == "_hr":
                        pt_id = pt["ID"][:-3]
                    else:
                        pt_id = pt["ID"]
                       
                 
                    # search for heterologous reaction in database
                    if not(self.heterologous_reaction_database_model):
                        # no heterologous reaction database model available, create empty model
                        warn("No heterologous reaction database model provided", UserWarning)
                        self.heterologous_reaction_database_model = cobra.Model()
                        
                    if pt_id not in [rxn.id for rxn in self.heterologous_reaction_database_model.reactions]:
                        # reaction already exists? if yes copy
                        if pt_id in [rxn.id for rxn in self.model.reactions]:
                            print(pt_id + " already in native model. Copy reaction with new bounds")
                            rxn_to_add = _create_reaction(pt["ID"]+"_COPY",
                                             pt["ID"],
                                             pt["lower_bound"],
                                             pt["upper_bound"],
                                             self.model.reactions.get_by_id(pt_id).metabolites)

                        else:                               
                            raise Exception("Heterologous reaction " + pt_id + " not found in database model" )
                    else:                            
                        # extract reaction from database model
                        rxn_from_database = self.heterologous_reaction_database_model.reactions.get_by_id(pt_id)
                        # if not provided use bounds from heterologous reaction database
                        if "lower_bound" not in pt:
                            lb = rxn_from_database.lower_bound
                        else:
                            lb = pt["lower_bound"]
                            
                        if "upper_bound" not in pt:
                            ub = rxn_from_database.upper_bound
                        else:
                            ub = pt["upper_bound"]
                            
                        # create reaction
                        rxn_to_add = _create_reaction(pt["ID"],
                                             rxn_from_database.name,
                                             lb,
                                             ub,
                                             rxn_from_database.metabolites)
                            
                        
                        
                    # add reaction to model
                    self.model.add_reaction(rxn_to_add)

                    # save to reverse addin
                    self._reverse_design[pt_name] = {"ID": pt["ID"],
                                                     "type": pt["type"]}
                    
                else:
                    raise Exception("Unknown design or perturbation type: " + pt.type)
                    
            
        else:
            # no interventions provided in strain design solution
            if verbose: print("\tNo design interventions provided")
            pass
        
        # save current design
        self._design = design
        # save name of design
        self._applied_design = name
            
      
            
    def save_current_medium(self):
        """
        Save current medium composition, i.e., save the lower bounds of exchange reactions

        Returns
        -------
        None.

        """
        
        # save old medium bounds
        old_medium = self.medium_from_model(self.model)
        # old_medium = {}
        # for ex in self.model.boundary:
        #     if "EX_" in ex.id and ex.lower_bound < 0:
        #         old_medium[ex.id] = {"exchange_reaction_id": ex.id,
        #                              "lower_bound": ex.lower_bound}
                
        # save medium
        self._reverse_design["medium"] = old_medium
        
          
    def revert_strain_design(self, verbose=True):
        # if current design exist, revert
        
        if not(self._reverse_design):
            # no design currently applied
            pass        
        else:
            if verbose: print("Reverse previous strain design...")
            # revert interventions
            for pt_name, pt in self._reverse_design.items():
                if pt_name == "medium":
                    continue
                # get reaction
                pt_rxn = self.model.reactions.get_by_id(pt["ID"])
                # process deletion, cofeed, mediareduction targets
                if pt["type"] == "deletion" \
                        or pt["type"] == "mediareduction":
                    # restore original (unperturbed) reaction bounds
                    pt_rxn.bounds = (pt["lower_bound"], pt["upper_bound"])
                 
                # process addin targets
                elif pt["type"] == "addin" \
                    or pt["type"] == "source" \
                        or pt["type"] == "cofeed":
                    # remove reaction
                    self.model.remove_reactions([pt_rxn])

                else:
                    raise Exception("Unknown design or perturbation type: " + pt["type"])
                    
            # revert medium definition
            if "medium" in list(self._reverse_design.keys()):
                medium = self._reverse_design["medium"].copy()
                for ex in self.model.boundary:
                    if "EX_" in ex.id and ex.lower_bound < 0:
                        ex.lower_bound = 0
                # set previous medium composition
                for key, medium_component in medium.items():
                    # get reaction
                    ex_rxn = self.model.reactions.get_by_id(medium_component["exchange_reaction_id"])
                    # set new bound
                    ex_rxn.lower_bound = medium_component["lower_bound"]
                    
            self._reverse_design = {}
            self._reverse_design["medium"] = medium
            self._design = {}
            self._applied_design = None
     
        
    def mutable_flux_at_maximum_objective(self):
        # calculate fluxes of all interventions of the design
        if not(self.design):
            print("No strain design loaded. Provide a strain design first.")
            return None
        
        # mutable_ids = [mutable["ID"] for mutable in self.design.values()]
        # check IDs
        mutable_ids = []
        for mutable in self.design.values():
            if "__hr" in mutable["ID"]:
                mutable_ids.append(mutable["ID"][:-4])
            elif "_hr" in mutable["ID"]:
                mutable_ids.append(mutable["ID"][:-3])
            else:
                mutable_ids.append(mutable["ID"])
                
        return self.flux_at_maximum_objective(mutable_ids)
    
      
    def reaction_activity_changes(self, original_model=None, medium={}):
        """
        Determine reactions that are activated/disactivated in strain design at maxmum growth

        Returns
        -------
        None.

        """
        
        
        
        
        if not(not(original_model)):
            model_wt = original_model
                
        elif not(not(self._model_original)):
            model_wt = self._model_original
            
        else:
            # no original wildtype model available
            print("No wildtype model provided for calculating the reference state")
            
            return None
        
        # parameter
        flux_thres = 1e-05
        
        
        # protect wildtype model
        with model_wt:
            # prepare wildtype model
            # set medium
             # set new medium composition
            for key, medium_component in medium.items():
                # get reaction
                ex_rxn = self.model_wt.reactions.get_by_id(medium_component["exchange_reaction_id"])
                # set new bound
                ex_rxn.lower_bound = medium_component["lower_bound"]
                       
            # calculate wildtype flux distribution
            sol = model_wt.optimize()
            fluxes_wt = sol.fluxes
            
        # calculate mutant flux distribution
        # protect model
        with self.model as model:
            sol = model.optimize()
            fluxes_mutant = sol.fluxes
            
            
        # compare solutions and extract reaction activity changes
        rxn_activated = {}
        rxn_deactivated = {}
        for rxn_id, flux_wt in fluxes_wt.items():
            if abs(flux_wt) < flux_thres and abs(fluxes_mutant[rxn_id]) >= flux_thres:
                rxn_activated[rxn_id] = fluxes_mutant[rxn_id]
                
            elif abs(flux_wt) >= flux_thres and abs(fluxes_mutant[rxn_id]) < flux_thres:
                rxn_deactivated[rxn_id] = flux_wt
                
        
        return pd.Series(rxn_activated), pd.Series(rxn_deactivated)  

       

    def build_design_dict(self, reaction_ids, intervention_types, bounds,
                          medium={},
                          wildtype_bounds=[],
                          integrate_design=True, design_name=None):
        """
        Build a design dictionary

        Parameters
        ----------
        reaction_ids : TYPE
            list of target reaction IDs.
        intervention_type : TYPE
            list of corresponding intervention types. Available: deletion, mediareduction, cofeed, source, addin
        bounds : list
            DESCRIPTION. 
        medium : dict
            lower bounds of boundary reactions resembling medium components (cf. medium_from_model)
        wildtype_bounds : list, optional
            bounds of of target in a non-perturbed, wildtype state
        integrate_design: boolean, optional
            Is new design directly integrated as a solution in strain design model?

        Returns
        -------
        None.

        """
             
        if isinstance(reaction_ids, str):
            reaction_ids = [reaction_ids]
        
        if isinstance(intervention_types, str):
            # use type for all provided reactions
            intervention_types = [intervention_types for i in range(len(reaction_ids))]
            
        if isinstance(intervention_types, tuple):
            # use bounds for all provided reactions
            bounds = [bounds for i in range(len(reaction_ids))]
        
        # check input consistency
        if not(len(reaction_ids) == len(intervention_types) \
            and len(reaction_ids) == len(bounds) \
            and (len(reaction_ids) == len(wildtype_bounds) or len(wildtype_bounds) == 0)):
            raise Exception("Number of provided reaction IDs, intervention types or bounds do not match!")
            
        # build design dictionary from inputs
        interventions_dict = {}
        for i in range(len(reaction_ids)):
            interventions_dict[reaction_ids[i]] = {
                    "ID": reaction_ids[i],
                    "type": intervention_types[i],
                    "lower_bound": bounds[i][0],
                    "upper_bound": bounds[i][1]
                    }
            if not(len(wildtype_bounds) == 0):
                interventions_dict[reaction_ids[i]]["wildtype_bounds"] = wildtype_bounds[i]
        
        # build medium dictionary
        # load current medium
        if not(medium):
            medium = self.medium_from_model(self.model)
        
        design_dict = {"interventions": interventions_dict,
                       "medium": medium
                       }
                
                
        if integrate_design:
            if not(design_name):
                design_name = "_".join(list(interventions_dict.keys()))
                keep_design_key = False
            else:
                keep_design_key = True
            # parse design solution   
            self._parse_solutions({design_name: design_dict}, keep_design_key=keep_design_key)
            
                
        return design_dict
    
    
    
    
        
    def determine_significant_designs_from_solutions(self, target_reaction=None,
                                                     design_objective=None,
                                                     minimum_significance=0.1,
                                                     eval_gpr=False,
                                                     adopt_significant_designs=False):
        """
        scan all parsed design solutions for significant subsets of interventions

        Parameters
        ----------
        target_reaction : TYPE, optional
            DESCRIPTION. The default is None.
        design_objective : str, optional
            design objective by which strain designs were calculated. The default is None.
                "growth_coupling": Minimization of target reaction flux at a fix growth rate
                None: currently active objective function in model
        minimum_significance : float, optional
            minimum allowable fraction of objective function value for reduced design               
        adopt_significant_designs : bool, optional
            parse significant designs and delete the original set. The default is False.

        Returns
        -------
        designs_significant_dict : dict
            dictionary of significant designs.
            
        designs_significant_frame : DataFrame
            dataframe of significant designs

        """
        
        
        # check target reaction
        if not(target_reaction):
            if not(hasattr(self, "target_reaction")):
                raise Exception("No target reaction provided.")
            else:
                target_reaction = self.target_reaction
        
        # check design objective
        if design_objective not in ["growth_coupling", "target_at_max_growth"]:
            print("Unknown design objective", design_objective)
            print("\tSet design objective to growth_coupling")
            design_objective = "growth_coupling"
        
        
        
        # get current solutions and reset solutions
        solutions = self.solutions.__dict__.copy()
        # if adopt_significant_designs:
        #     # save old solution space
        #     self._solutions_deprecated = deepcopy(self.solutions)
        #     # initiate new solution space
        #     self.solutions = SimpleNamespace()
            
        # get significant designs from currently available solutions      
        count = 0
        designs_significant_dict = {}
        design_significant_list = []
        designs_config = {}
        designs_filename = {}
        for key_sol, sol in solutions.items():
            
            # load design solution
            sol.load()
            # determine significant subsets of currently loaded design
            design_subsets, smallest_design = self.determine_significant_designs(
                design_objective,
                target_reaction,
                eval_gpr=eval_gpr,
                minimum_significance=minimum_significance
                )
            
            # # calculate significant subset of current design
            # # use design objective specific parameters
            # if design_objective == "growth_coupling":  
            #     # check available config parameter
            #     parameter_needed = ["growth_rate_fix", "biomass_id"]
            #     design_kwargs = {}
            #     if not(all([hasattr(sol._config, parameter) for parameter in parameter_needed])):
            #         warn("config parameter for design objective " + design_objective \
            #              + " not available", UserWarning)
            #     else:                                
            #         # set gcOpt objective function
            #         gr_tol = 1e-5 # allow for a small tolerance on the fix growth rate, complement gcOpt solver tolerance
            #         design_kwargs["objective"] = target_reaction
            #         design_kwargs["objective_direction"] = "min"
            #         design_kwargs["additional_constraints"] \
            #             = {sol._config.biomass_id: 
            #                (sol._config.growth_rate_fix-gr_tol, sol._config.growth_rate_fix+gr_tol)}
                            
            # elif design_objective == "target_at_max_growth": 
            #     # design objective is to maximize the minimal target reaction rate at maximum growth (RobustKnock)
            #     # check available config parameter
                
            #     parameter_needed = ["biomass_id"]
            #     design_kwargs = {}
            #     if not(all([hasattr(sol._config, parameter) for parameter in parameter_needed])):
            #         warn("config parameter for design objective " + design_objective \
            #              + " not available", UserWarning)
            #     else:                                
            #         # set gcOpt objective function
            #         design_kwargs["objective"] = sol._config.biomass_id
            #         design_kwargs["objective_direction"] = "max"
                            
                            
            # # calculate design objective for all combinations of full design
            # # design_subsets, smallest_design = self.design_objects_significance(
            # #     target_reaction=target_reaction,
            # #     eval_gpr=eval_gpr,
            # #     **design_kwargs
            # #     )
            
            # design_subsets, smallest_design = design_significance(
            #     self.model,
            #     self._design,
            #     self._reverse_design,
            #     target_reaction=target_reaction,
            #     eval_gpr=eval_gpr,
            #     **design_kwargs
            #     )
            
            

            
            # check significance results
            if design_subsets.empty:
                # even the full design is infeasible
                warn("Design is infeasible: " + str(list(sol._design_solution["interventions"].keys())), UserWarning)
                continue
            
            # extract significant design subsets
            design_significant_dataframe = design_subsets.loc[design_subsets.loc[:, "significance"]>minimum_significance, :]
            # design_significant_list.append(design_significant_dataframe)
            
            # check if any significant solution exist
            if design_significant_dataframe.empty:
                warn("No significant design solutions: " + str(list(sol._design_solution["interventions"].keys())),
                     UserWarning)
                continue

            
            # get smallest design for similar significance
            bin_num = (design_significant_dataframe["significance"].max()\
                       -design_significant_dataframe["significance"].min())/0.01
            if bin_num < 1:
                bin_num = 1
            else:
                bin_num = round(bin_num)
                
            design_groups = design_significant_dataframe.groupby(
                pd.cut(design_significant_dataframe["significance"], bin_num)
                )
            
            
            sign_idx = [] # indices of significant, smallest designs
            for dt in design_groups:
                df = dt[1] # get dataframe
                if not df.empty:
                    # get smallest design in group
                    design_len = []
                    design_idx = []
                    for i, d in df.iterrows():
                        design_len.append(len(d.loc["interventions"]))
                        design_idx.append(i)
                    sign_idx.append(design_idx[design_len.index(min(design_len))])
                        
            # save significant design
            design_significant_list.append(design_significant_dataframe.loc[sign_idx])
            
            for i in sign_idx:
                # construct design solution
                design_i_subset = design_significant_dataframe.loc[i, "interventions"]
                design_sign = deepcopy(sol._design_solution)
                design_i = deepcopy(design_sign["interventions"])
                design_sign["interventions_excluded"] = design_significant_dataframe.loc[i, "excluded"]
                for key, obj in design_i.items():
                    if obj["ID"] not in design_i_subset:
                        del(design_sign["interventions"][key])

                # renew objective function value
                design_sign["objective_value"] = design_significant_dataframe.loc[i, "objective_value"]
                        
                # save design as dict
                designs_significant_dict[self.design_key.format(count)] = design_sign
                # save additional data for parsing
                designs_filename[self.design_key.format(count)] = Path(sol.filename)
                designs_config[self.design_key.format(count)] = sol._config
                count += 1
                
        
        
            
        # check if any significant designs were generated
        if not(design_significant_list):
            # no significant design encountered
            warn("Could not determine significant designs", UserWarning)
            return designs_significant_dict, None
 
        # merge significant design in datafrme
        designs_significant_frame = pd.concat(design_significant_list, ignore_index=True)    
        # reload design solutions   
        if adopt_significant_designs:
            # save old solution space
            self._solutions_deprecated = deepcopy(self.solutions)
            # initiate new solution space
            self.solutions = SimpleNamespace()
            # parse new, significant design solutions
            self._parse_solutions(
                designs_significant_dict,
                config=designs_config,
                filename=designs_filename,    
                check_duplicate=True,
                reduce_design=False,
                eval_gpr=False
                )
            
                
        
            
        
        return designs_significant_dict, designs_significant_frame
    
    
    def determine_significant_designs(
            self,
            design_objective: str,
            target_reaction: Optional[str] = None,
            eval_gpr: bool = False,
            minimum_significance: float = 0.1,
            parameter_objective: dict = {}
            ):
        """Determine significant design variable subsets from currently loaded design
        """
        
        # determine target reaction (objective)
        if not(target_reaction):
            if not(self.target_reaction):
                warn('No target reaction (objective) identified')
                return
        
        
        # get config of design solution
        if self._applied_design in self.solutions.__dict__:
            sol_config = self.solutions.__dict__[self._applied_design]._config
        else:
            warn('Design solution not found')
            return
        
        # calculate significant subset of current design
        # use design objective specific parameters
        if design_objective == "growth_coupling":  
            # check available config parameter
            parameter_needed = ["growth_rate_fix", "biomass_id"]
            design_kwargs = {}
            if not(all([hasattr(sol_config, parameter) for parameter in parameter_needed])):
                warn("config parameter for design objective " + design_objective \
                     + " not available", UserWarning)
            else:                                
                # set gcOpt objective function
                gr_tol = 1e-5 # allow for a small tolerance on the fix growth rate, complement gcOpt solver tolerance
                design_kwargs["objective"] = target_reaction
                design_kwargs["objective_direction"] = "min"
                design_kwargs["additional_constraints"] \
                    = {sol_config.biomass_id: 
                       (sol_config.growth_rate_fix-gr_tol, sol_config.growth_rate_fix+gr_tol)}
                        
        elif design_objective == "target_at_max_growth": 
            # design objective is to maximize the minimal target reaction rate at maximum growth (RobustKnock)
            # check available config parameter
            
            parameter_needed = ["biomass_id"]
            design_kwargs = {}
            if not(all([hasattr(sol_config, parameter) for parameter in parameter_needed])):
                warn("config parameter for design objective " + design_objective \
                     + " not available", UserWarning)
            else:                                
                # set gcOpt objective function
                design_kwargs["objective"] = sol_config.biomass_id
                design_kwargs["objective_direction"] = "max"
                        
                        
        # calculate design objective for all combinations of full design
        # design_subsets, smallest_design = self.design_objects_significance(
        #     target_reaction=target_reaction,
        #     eval_gpr=eval_gpr,
        #     **design_kwargs
        #     )
        
        design_subsets, smallest_design = design_significance(
            self.model,
            self._design,
            self._reverse_design,
            target_reaction=target_reaction,
            eval_gpr=eval_gpr,
            minimum_relative_significance=minimum_significance,
            **design_kwargs
            )
        
        return (design_subsets, smallest_design)
        
    

    
    

            
    def determine_escape_metabolites(self, target_reaction=None, biomass_reaction=None,
                                     fixed_growth_rate=0.2, number_metabolites=None,
                                     save=False, filename=None,
                                     ) -> pd.DataFrame:
        """Determine metabolites which diminish the design objective when being spiked
        into the metabolic model
        The measure of an escape potential is compute based on
            (1) shadow prices at optimality
            (2) objective value changes when metabolite is added
        For both cases, maximization of growth and minimization of the target reaction
        flux at a fixed growth rate are applied as objective functions
        """
        
        # get growth-coupling target reaction 
        if not target_reaction:
            if not self.target_reaction:
                warn("No target reaction specified!", UserWarning)
                return pd.DataFrame(), {}
            else:
                target_reaction = self.target_reaction
                
        if not(biomass_reaction):
                # identifiy biomass formation reaction objective function
                objective_rxn = []
                for rxn in self.model.reactions:
                    if rxn.objective_coefficient != 0:
                        objective_rxn.append(rxn.id)
                if len(objective_rxn) != 1:
                    raise Exception("Objective function does not include a single reaction")
                else:
                    biomass_reaction = objective_rxn[0]                    
        else:
            # reference reaction provided
            if biomass_reaction not in [rxn.id for rxn in self.model.reactions]:
                raise Exception("Reference reaction " + biomass_reaction + " is not in model")
                
                
        # calculate shadow prices
        
        # growth coupling objective
        with self.model as model:
            # change objective to target enzyme
            model.objective = target_reaction
            # minimize objective
            model.objective_direction = "min"
            # fix biomass
            model.reactions.get_by_id(biomass_reaction).lower_bound = fixed_growth_rate
            model.reactions.get_by_id(biomass_reaction).upper_bound = fixed_growth_rate
            # solve model
            sol_sp_t = model.optimize()
            # get and sort shadow prices
            # positive shadow price means objective value decreases when metabolite is provided
            sp_target = sol_sp_t.shadow_prices.sort_values()
        
        # growth objective
        with self.model as model:
            # change objective to biomass equation
            model.objective = biomass_reaction
            # change direction, maximize
            model.objective_direction = "max"
            # solve model
            sol_sp_g = model.optimize()
            # get and sort shadow prices
            # negative shadow price means objective value increases when metabolite is provided
            sp_growth = sol_sp_g.shadow_prices.sort_values()
        
        # calculate pFBA solution
        sol_fba = cobra.flux_analysis.pfba(self.model)    
        
        
        # Determine and save incremental escape potantial of metabolites
        # escape metabolites have a positive target reaction shadow price and 
        # a negative growth shadow price
        incr_escape_all_mets_dict = {}
        for met in self.model.metabolites:
            if "C" in met.elements and sp_target[met.id] > 0 and sp_growth[met.id] < 0:
                # calculate sum of fluxes for metabolite
                flux_sum = np.sum([abs(sol_fba.fluxes[r.id])/met.elements["C"] for r in met.reactions])
                
                incr_escape_all_mets_dict[met.id] = [
                    met.name,
                    met.elements["C"],
                    flux_sum,
                    sp_target[met.id]/met.elements["C"],
                    sp_growth[met.id]/met.elements["C"],
                    -sp_target[met.id]*sp_growth[met.id]/(met.elements["C"]**2) # incremental escape potential
                    ]
        # make a frame
        incr_escape_all_mets_frame = pd.DataFrame(
            data=incr_escape_all_mets_dict.values(),
            index=incr_escape_all_mets_dict.keys(),
            columns=["name", "carbon_atoms", "sum_fluxes", "sp_target", "sp_growth",
                    "sp_product_target_growth"]
            )
            
            
        # compute macroscopic escape potential
        incr_escape_all_mets_frame = incr_escape_all_mets_frame.sort_values(by="sp_product_target_growth", ascending=False)
        rxn_to_check = [r.id for r in self.model.reactions if abs(sol_fba.fluxes[r.id]) < 1e-9]
        C_mol_rate = 10
        flux_cutoff = 1e-3
        
        if not number_metabolites:
            number_metabolites = len(incr_escape_all_mets_frame)
        
        macro_escape_mets_dict = {}
        escape_mets_frame = pd.DataFrame(
            columns=list(incr_escape_all_mets_frame.columns)+["sp_target_macro", "sp_growth_macro",
                                                          "sp_product_target_growth_macro",
                                                          "maximum_growth_rate",
                                                          "minimum_target_rate",
                                                          "activated_reactions"])

        
        for i in range(number_metabolites):
           met = incr_escape_all_mets_frame.iloc[i].name
           row = incr_escape_all_mets_frame.iloc[i, :]
           
           with self.model as model:
               # constrain demand flux
               lb_met = -C_mol_rate/model.metabolites.get_by_id(met).elements["C"]
               # create demand reaction
               model.add_boundary(model.metabolites.get_by_id(met), type="sink",
                                  reaction_id="DM_"+met, lb=lb_met)
               
               
               # optimize for maximum growth
               sol_g = cobra.flux_analysis.pfba(model)
               sol_g_obj_val = sol_g.fluxes[biomass_reaction]
       
               # optimize gcOpt objective
               model.objective = target_reaction
               model.objective_direction = "min"
               model.reactions.get_by_id(biomass_reaction).lower_bound = fixed_growth_rate
               model.reactions.get_by_id(biomass_reaction).upper_bound = fixed_growth_rate
               sol_t = cobra.flux_analysis.pfba(model)
               sol_t_obj_val = sol_t.fluxes[target_reaction]
               
               # get activated reactions
               rxn_activated = []
               for rxn_id in rxn_to_check:
                   if (abs(sol_g.fluxes[rxn_id]) > flux_cutoff or abs(sol_t.fluxes[rxn_id]) > flux_cutoff):
                       rxn_activated.append(rxn_id)
                       
                       
           # save in dictionary
           macro_escape_mets_dict[met] = [
               (sol_t_obj_val-sol_sp_t.objective_value)/C_mol_rate, # macroscopic target shadow price 
               (sol_g_obj_val-sol_sp_g.objective_value)/C_mol_rate, # macroscopic growth shadow price 
               (sol_t_obj_val-sol_sp_t.objective_value)*(sol_g_obj_val-sol_sp_g.objective_value)/(C_mol_rate**2), # macroscopic escape potential
               sol_g_obj_val, # maximum growth rate
               sol_t_obj_val, # minimum target flux rate
               ";".join(rxn_activated) # activated reactions
               ]
           
           # merge with incremental results
           for idx in row.index:
               escape_mets_frame.loc[met, idx] = row[idx]
           
           escape_mets_frame.loc[met, "sp_target_macro"] = (sol_t_obj_val-sol_sp_t.objective_value)/C_mol_rate
           escape_mets_frame.loc[met, "sp_growth_macro"] = (sol_g_obj_val-sol_sp_g.objective_value)/C_mol_rate
           escape_mets_frame.loc[met, "sp_product_target_growth_macro"] = -(sol_t_obj_val-sol_sp_t.objective_value)*(sol_g_obj_val-sol_sp_g.objective_value)/(C_mol_rate**2)
           escape_mets_frame.loc[met, "maximum_growth_rate"] = sol_g_obj_val
           escape_mets_frame.loc[met, "minimum_target_rate"] = sol_t_obj_val
           escape_mets_frame.loc[met, "activated_reactions"] = ";".join(rxn_activated)
         
           
        escape_mets_frame = escape_mets_frame.sort_values(by="sp_product_target_growth_macro", ascending=False)
        
        
        # save frame
        if save:
            if not filename:
                filename = "escape_metabolites_results"
            escape_mets_frame.to_excel(filename+".xlsx", engine="openpyxl",
                                       sheet_name="escape_metabolites")
         
        return escape_mets_frame

    
    
    def add_cofactor_regeneration(self, cofactor):
        
        

        # choose cofactor regeneration to be added
        if cofactor == "tetrahydrobiopterin":
            # add tetrahydrobiopterin regeneration
            self.model = add_prebuilt_pathways.tetrahydrobiopterin(self.model)
            
        else:
            print("Unknown cofactor identifier: ", cofactor)



    def growth_coupling_summary(self, results_dir=None, results_filename="gcOpt_solution_summary",
                                determine_significant_designs=False, 
                                minimum_significance=0.1,
                                save_results=False, save_flux_space_plots=False,
                                eval_gpr=False, design_objective="growth_coupling"
                                ):
        """
        analyse all solutions regarding growth-coupling and summarize findings
        Load all consistent results files in dedicated directory and analyze solutions accordingly

        Parameters
        ----------
        results_dir : str, optional
            Path to storing the results file. The default is "".
            
        results_filename : str, optional
            filename for saving growth coupling summary result. The default is "gcOpt_solution_summary".
            
        save_results : boolean, optional
            save summary results to file. The default is False
                            
        save_flux_space_plots : boolean, optional
            Plot and save 2D projections of flux spaces. The default is True.
            
        determine_significant_designs : boolean, optional
            determine subset of significant interventions from design solutions. The default is False
            
        eval_gpr : boolean, optional
            Take GPR relations into account when determining significant designs. The default is False

        Returns
        -------
        summary : DataFrame
            Summary data of growth-coupling solutions.

        """
        
        self.revert_strain_design()

        # check if design solutions exist
        if not(self.solutions.__dict__):
            warn("No design solutions available!", UserWarning)
            return

        # set up results directory
        if not(results_dir):
            results_dir = getcwd()

        # check directory
        results_dir = Path(results_dir)
        if not(results_dir.is_dir()):
            warn("Not a valid path to a directory: " + results_dir, UserWarning)
            results_dir = Path(getcwd())
                       
        # set up directory for saving flux space plots
        if save_flux_space_plots:
            plots_dir = results_dir.joinpath("flux_space_plots")
            if not(plots_dir.is_dir()): mkdir(plots_dir)
            

        # get significant designs
        if determine_significant_designs:
            self.determine_significant_designs_from_solutions(
                target_reaction=self.target_reaction,
                design_objective=design_objective,
                adopt_significant_designs=True,
                minimum_significance=minimum_significance,
                eval_gpr=eval_gpr
                )

        # save strain design solutions
        if save_results:
            self.save_strain_design_file(
                results_dir.joinpath(results_filename + "_dict.pickle")
                )
            

        
        solution_index = list(self.solutions.__dict__)
        columns = ["filename",
                   "key_in_file",
                    "objective_value",
                     "number_interventions",
                     "number_genetic_interventions",
                     "interventions",  
                     "interventions_excluded",
                     "genetic_interventions",
                     "unique_genetic_interventions",
                     "carbon_uptake_bounds",
                     "mutables_flux_at_maximum_growth",
                     "coupling_strength",
                     "max_growth_rate",
                     "target_flux_at_maximum_growth",
                     "score",
                    "biomass_precursor_auxotrophy"]
        
        # init data frame
        summary = pd.DataFrame(columns=columns, index=solution_index)
        
        # analyse each solution
        for sol_name, sol_class in self.solutions.__dict__.items():
            # load design
            sol_class.load()
            design = sol_class._design_solution
            
            # save filename
            if hasattr(sol_class, "filename"):
                summary.loc[sol_name, "filename"] = getattr(sol_class, "filename")
                
            # save key of solution in file
            if hasattr(sol_class, "_solution_key_in_file"):
                summary.loc[sol_name, "key_in_file"] = getattr(sol_class, "_solution_key_in_file")
            
            # get objective value
            if "objective_value" in design:
                 summary.loc[sol_name, "objective_value"] = design["objective_value"]
            
            # get interventions
            interventions_list = list(design["interventions"])
            summary.loc[sol_name, "interventions"] = ";".join(interventions_list)
            summary.loc[sol_name, "number_interventions"] = len(interventions_list)
            
            # get interventions excluded from original design
            if "interventions_excluded" in design:
                summary.loc[sol_name, "interventions_excluded"] = ";".join(design["interventions_excluded"])
            else:
                summary.loc[sol_name, "interventions_excluded"] = ""
            
            # translate reaction interventions to genetic interventions
            genetic_interventions_list = self.reactions_to_genes_design(interventions_list)         
            summary.loc[sol_name, "genetic_interventions"] = ";".join([",".join(gi) for gi in genetic_interventions_list])
            
            # count number of unique genetic interventions
            unq_genes = []
            for gi in genetic_interventions_list:
                for g in gi:
                    if g not in unq_genes:
                        unq_genes.append(g)
                        
            summary.loc[sol_name, "number_genetic_interventions"] = len(unq_genes)       
            summary.loc[sol_name, "unique_genetic_interventions"] = ";".join(unq_genes)       
            
            # extract carbon source(s)
            carbon_uptake_bounds = {}
            for ex_rxn in self.model.boundary:
                # is an uptake reaction?
                # Also consider source and cofeed reactions of the design
                if ex_rxn.lower_bound < 0:
                    # is a carbon source?
                    met = list(ex_rxn.metabolites.keys())[0]
                    if "C" in met.elements:
                        carbon_uptake_bounds[ex_rxn.id] = ex_rxn.lower_bound
            # summary.loc[sol_name, "carbon_uptake_bounds"] = [carbon_uptake_bounds]
            summary.loc[sol_name, "carbon_uptake_bounds"] = ";".join([rxn + ":" + str(flux) for rxn, flux in carbon_uptake_bounds.items()])
            
            # calculate growth coupling strength
            summary.loc[sol_name, "coupling_strength"] = self.coupling_strength(self.target_reaction)
            
            # calculate target reaction flux at maximum growth rate
            summary.loc[sol_name, "target_flux_at_maximum_growth"] \
                = self.flux_at_maximum_objective([self.target_reaction])[self.target_reaction]
                
            # mutable reaction fluxes at maximum growth rate
            summary.loc[sol_name, "mutables_flux_at_maximum_growth"] \
                = ",".join([rxn + ":" + str(flux) for rxn, flux in self.flux_at_maximum_objective(interventions_list).items()])
            
            # maximum growth rate
            summary.loc[sol_name, "max_growth_rate"] = self.model.slim_optimize()
            
        
            # calculate whoch biomass precursors cannot be synthesized without active target reaction
            precursor_fluxes, biomass_flux_rescue = self.precursor_availability(disable_reactions=[self.target_reaction])
            precursor_auxotrophy = []
            for met, biomass_flux in biomass_flux_rescue.items():
                if biomass_flux > 0:
                    # precursor cofeed enables growth
                    precursor_auxotrophy.append(met)
            
            summary.loc[sol_name, "biomass_precursor_auxotrophy"] = ";".join(precursor_auxotrophy)
            
            # plot flux space projection
            if save_flux_space_plots:
                save_name = plots_dir.joinpath(sol_name)
                fp = self.flux_space_projection(self.target_reaction, plot=True, save_name=save_name)
                plt.close("all")
            
        # assign a score to designs
        for i, row in summary.iterrows():
            if not(row.loc["coupling_strength"]):
                # infeasible growth coupling design solution
                row.loc["score"] = 0
                continue
            
            row.loc["score"] = (
                row.loc["coupling_strength"]
                + row.loc["max_growth_rate"]
                ) / np.log10(row.loc["number_interventions"])
            

        # save summary to excel file
        if save_results:
            
            # helper function 
            def add_excel_table(xl_sheet, dataframe):
                # create an excel table
                (max_row, max_col) = dataframe.shape
                # determine headers for the table
                column_settings = [{"header": dataframe.index.name}]
                for header in dataframe.columns:
                    column_settings.append({'header': header})
                # add excel table
                xl_sheet.add_table(0, 0, max_row, max_col, {'columns': column_settings})
                xl_sheet.set_column(0, max_col, 12)
            
            
            
            file_name = "{0}.xlsx".format(results_filename)
            file_path = results_dir.joinpath(file_name)
            
            
            
            sheet_name = "Summary"
            if Path.is_file(file_path):
                # delete old summary
                with pd.ExcelWriter(
                        file_path,
                        engine='openpyxl', 
                        mode='a',
                        if_sheet_exists='replace',
                ) as writer:
                    pd.DataFrame().to_excel(writer, sheet_name=sheet_name)
                    writer.save()
                        
                # save new summary 
                with pd.ExcelWriter(
                        file_path,
                        engine='xlsxwriter', 
                ) as writer:
                    # save summary to excel
                    summary.to_excel(writer, sheet_name=sheet_name)
                    # create excel table
                    add_excel_table(writer.sheets[sheet_name], summary)
              

                
            else:
                with pd.ExcelWriter(
                        file_path,
                        engine='xlsxwriter',
                ) as writer: 
                    # save summary to excel
                    summary.to_excel(writer, sheet_name="Summary")
                    # create excel table
                    add_excel_table(writer.sheets[sheet_name], summary)
                    
    

           
            
        return summary
    
    
    
    def reactions_to_genes_design(self, interventions_list: list) -> list:
        """
        determine necessary genetic interventions from a list of reaction interventions

        Parameters
        ----------
        interventions_list : list
            list of reaction identifier

        Returns
        -------
        genetic_interventions_list : list
            list of gene identifiers

        """
        
        genetic_interventions_list = []
        genetic_interventions_sets = []
        
        # protect model
        # with self.model as model:  
        #     with self.heterologous_reaction_database_model as model_hr:
        for ri in interventions_list:
            if (ri in self.model.reactions) and (not(ri.endswith("__hr"))):
                # get corresponding genes for reaction deletion
                gi = set([g.id for g in self.model.reactions.get_by_id(ri).genes])                        
                    
            elif ri.endswith("__hr"):
                # heterologous reaction encountered                       
                if ri.replace("__hr", "") in self.heterologous_reaction_database_model.reactions:
                    # get corresponding genes for reaction
                    gi = set([g.id+"__hr" for g in self.heterologous_reaction_database_model.reactions.get_by_id(ri.replace("__hr", "")).genes])

            # check if genes set was already encountered
            if (len(gi) > 0) and (gi not in genetic_interventions_sets):
                genetic_interventions_sets.append(gi)
            
        # construct genetic interventions list
        genetic_interventions_list = [list(gi) for gi in genetic_interventions_sets]
        
        return genetic_interventions_list
    
    
    
    def medium_from_model(self, model):
        """Extract medium composition from model
            medium is saved as lower bounds of boundary reactions
        """

        medium = {}
        for ex in model.boundary:
            if "EX_" in ex.id and ex.lower_bound < 0:
                medium[ex.id] = {"exchange_reaction_id": ex.id,
                                 "lower_bound": ex.lower_bound}
        
        return medium

            
    # def _reduce_to_smallest_significant_subset(self, design_solution):
    #     """
    #     analyze strain design object for the significant subsets of interventions
    #     adopt signifcant subset

    #     Parameters
    #     ----------
    #     design_solution : dict
    #         strain design object.

    #     Returns
    #     -------
    #     design_solution : TYPE
    #         strain design object of only significant subset of the original design.

    #     """
    #     # strip unnecessary interventions from design 
        
    #     if not(design_solution["interventions"]):
    #         warn("No interventions in design solution. Design reduction is not possible.", UserWarning)
    #         return design_solution
              
    #     # load design
    #     self.load_strain_design(
    #         design=design_solution["interventions"],
    #         medium=design_solution["medium"],
    #         verbose=False
    #         )
        
        
    #     # get smallest significant subset of the design
    #     subset_sign, design_reduced = self.design_objects_significance()
               
    #     # revert strain design
    #     self.revert_strain_design(verbose=False)
        
    #     if design_reduced:
    #         # return reduced design
    #         print("Reduced significant subset of design found")
    #         design_solution["interventions"] = design_reduced
            
    #     return design_solution
        
        
            
    def _parse_solutions(self, design_solutions, config=SimpleNamespace(), default_design={},
                         filename=Path(), check_duplicate=True,
                         reduce_design=False, keep_design_key=False,
                         eval_gpr=False):
        """
        parse and save design solution

        Parameters
        ----------
        design_solutions : dict
            strain design object.
        config : SimpleNamespace, optional
            config parameters of optimization algorithm. The default is SimpleNamespace().
        default_design : dict, optional
            set of interventions already applied in the current model. The default is {}.
        filename : str, optional
            filename of stored strain design solution. The default is Path().
        check_duplicate : bool, optional
            exclude duplicate solutions. The default is True.
        reduce_design : TYPE, optional DEPRECATED
            determine and parse only significant subsets of the design. The default is False.
        keep_design_key : Boolean, optional
            use keys in desing_solutions dict for keys in SimpleNamespace, The default is False.

        Returns
        -------
        None.

        """
        # 
        # check provided config and filenames
        if isinstance(config, dict):
            # check length of dictionary
            if len(design_solutions) != len(config):
                warn("Numbers of design solutions and config do not match", UserWarning)
                config=SimpleNamespace()
                
        if isinstance(filename, dict):
            # check length of dictionary
            if len(design_solutions) != len(filename):
                warn("Numbers of design solutions and filenames do not match", UserWarning)
                filename=Path()          
    
        
        if len(design_solutions) > 0:
            for sol_key, sol in design_solutions.items():
                # disregard empty solutions
                if not(sol["interventions"]) and not(sol["medium"]):
                    # empty solution
                    print("\tEmpty design solution encountered")
                    continue
                
                # add default design interventions
                if len(default_design) > 0:
                    sol = self._integrate_default_design(sol, default_design)
                
                # check for linked reaction deletions due to GPR dependencies and add them to design
                if eval_gpr:
                    reaction_id_to_add = []
                    for i_key, i_dict in sol["interventions"].items():
                        # only evaluate deletions
                        if not(i_dict["type"] == "deletion"): continue
                        # check for linked reactions
                        rxns_linked = util.GPR_linked_reaction_knockouts(
                                            self.model,
                                            [i_dict["ID"]],
                                            eval_gpr=True
                                            )
                        # add linked reactions to design
                        for r in rxns_linked:
                            # make sure linked reaction is not already in design
                            if (r not in [del_obj["ID"] for del_obj in sol["interventions"].values()]) \
                                and (r not in reaction_id_to_add):
                                    
                                reaction_id_to_add.append(r)
                     
                    # add new design objects to design
                    if len(reaction_id_to_add) > 0:       
                        design_objects_to_add = self.build_design_dict(
                                                    reaction_id_to_add,
                                                    ["deletion" for t in range(len(reaction_id_to_add))],
                                                    [[0, 0] for t in range(len(reaction_id_to_add))],
                                                    integrate_design=False,
                                                    )
                        sol["interventions"] = {**sol["interventions"], **design_objects_to_add["interventions"]}
                        sol["GPR_linked_interventions_added"].extend(list(design_objects_to_add["interventions"].keys()))
                    
                
                # # get smallest significant subset of design
                # if reduce_design:
                #     sol = self._reduce_to_smallest_significant_subset(sol)
                
                # check for solution duplicates
                if check_duplicate:
                    if self._is_design_duplicate(sol):
                        # duplicate design solution
                        print("\tDuplicate design solutions encountered")
                        continue
                
                # specify key for new solution
                if keep_design_key:
                    key = sol_key
                else:
                    key = self.design_key.format(str(len(self.solutions.__dict__)+1))
                
                # construct solution
                if isinstance(config, dict):
                    config_parse = config[sol_key]
                else:
                    config_parse = config
                if isinstance(filename, dict):
                    filename_parse = filename[sol_key].name
                else:
                    filename_parse = filename.name  
                    
                self.solutions.__setattr__(
                    key, Solution(
                        self,
                        sol,
                        config=config_parse,
                        filename=filename_parse,
                        solution_key=key,
                        key_in_file=sol_key
                        )
                    )

        else:
            print("No solutions in provided dictionary")
        
     
        
    def _integrate_default_design(self, design_solution, default_design):
        """
        integrate default design into the strain design object
        a default design must be activated in the model

        Parameters
        ----------
        design_solution : dict
            strain design object.
        default_design : dict 
            set of interventions already applied in the current model.
            
        Raises
        ------
        Exception
            default intervention reaction is not in the model or heterologous reaction database.

        Returns
        -------
        design_solution : dict
            strain design object with integrate default design.

        """
        
        
        default_design_copy = deepcopy(default_design)
        
        # parameter
        addin_suffix = "__hr"
        
  
        
        # check if heterologous reaction database exists
        if hasattr(self, "heterologous_reaction_database_model"):
            if not(self.heterologous_reaction_database_model):
                # create empty database model
                self.heterologous_reaction_database_model = cobra.Model()
        else:
            # create empty database model
            self.heterologous_reaction_database_model = cobra.Model()
        # use model shortcut
        hrd_model = self.heterologous_reaction_database_model
            
            
        for key, intervent in default_design_copy.items():
                          
            # distinguish type of intervention
            if intervent["type"] in ["addin"]:
                # addin targets
                
                # mark addin target
                if not(intervent["ID"].endswith("_hr")) and not(intervent["ID"].endswith("__hr")):
                    target_id = intervent["ID"] + addin_suffix
                else:
                    target_id = intervent["ID"] 
   
                if intervent["ID"] not in hrd_model.reactions \
                    and intervent["ID"] not in self.model.reactions:
                        raise Exception("Default addin target " + intervent["ID"] + " not found in model or database!")
   
                # if target is in model, delete and add to heterologous database model
                hrd_rxn = None
                if intervent["ID"] in self.model.reactions:
                    rxn = self.model.reactions.get_by_id(intervent["ID"])
                    hrd_rxn = self._create_reaction(intervent["ID"],
                                                    intervent["ID"],
                                                    rxn.lower_bound,
                                                    rxn.upper_bound,
                                                    rxn.metabolites)
                    
                    rxn.delete()
                    
                    # copy reaction to heterologous database model
                    if intervent["ID"] in hrd_model.reactions:
                        # delete reaction
                        hrd_model.reactions.get_by_id(intervent["ID"]).delete()
                        
                    hrd_model.add_reaction(hrd_rxn)
                 
                # if target is a deletion target in the design solution, delete
                design_iter = deepcopy(design_solution)
                for key, obj in design_iter["interventions"].items():
                    if obj["ID"] == intervent["ID"] and obj["type"] == "deletion":
                        del design_solution["interventions"][key]
                 
                # change ID in intervention dict
                intervent["ID"] = target_id
                
    
        
            elif intervent["type"] in ["deletion"]:
                # deletion targets
                
                target_id = intervent["ID"]
                # target must be found in model
                if intervent["ID"] not in self.model.reactions:
                    raise Exception("Default design target " + intervent["ID"] + " not found in model!")
                    
                if "wildtype_bounds" not in intervent:
                    raise Exception("Specify <wildtype_bounds> for default deletion targets")
                    
                # activate reaction
                rxn = self.model.reactions.get_by_id(intervent["ID"])
                rxn.lower_bound = intervent["wildtype_bounds"][0]
                rxn.upper_bound = intervent["wildtype_bounds"][1]
                       
            
            else:
                warn("Intervention type " + intervent["type"] + " unkown", UserWarning)
                continue
        
            # add intervention to design
            design_solution["interventions"] = {**design_solution["interventions"],
                                                target_id: intervent}
            
            
        return design_solution
            
    
    
    def _is_design_duplicate(self, design_solution):
        """
        check if design solution already exists in solution dict

        Parameters
        ----------
        design_solution : dict
            new design solution dict.

        Returns
        -------
        boolean
            if True, design_solution already exists.

        """
        
        
        def create_parameter_list(sol_dict, i_keys):
            # create list of parameter list of interventions
            par_list = []
            for key in i_keys:
                i = sol_dict[key]
                p_keys = set(list(i))
                par_list.extend([i[p] for p in p_keys])
            return par_list
        
        # get interventions and medium of input design solution
        i_new = design_solution["interventions"]
        m_new = design_solution["medium"]
        i_new_keys = list(i_new)
        m_new_keys = list(m_new)
        # create parameter list of input design solution
        i_new_par_list = create_parameter_list(i_new, i_new_keys)
        m_new_par_list = create_parameter_list(m_new, m_new_keys)
        
        is_intervent_duplicate = False
        is_medium_duplicate = False
        is_duplicate = False
        for key_sol, sol in self.solutions.__dict__.items():
            # get interventions and medium
            i_dict = sol._design_solution["interventions"]
            m_dict = sol._design_solution["medium"]
            
            # check for same interventions
            if set(i_new_keys) == set(list(i_dict)):
                # same target reactions, check type and bounds
                i_dict_par_list = create_parameter_list(i_dict, i_new_keys)
                # check if solutions share the same intervention parameters
                if i_dict_par_list==i_new_par_list:
                    is_intervent_duplicate = True
                    
            # check for same medium
            if set(m_new_keys) == set(list(m_dict)):
                # same target reactions, check type and bounds
                m_dict_par_list = create_parameter_list(m_dict, m_new_keys)
                # check if solutions share the same intervention parameters
                if m_dict_par_list==m_new_par_list:
                    is_medium_duplicate = True
            
            is_duplicate = is_medium_duplicate and is_intervent_duplicate        
            if is_duplicate:
                break

        return is_duplicate
                
                

    def _create_reaction(self, rxn_id, rxn_name, lower_bound, upper_bound,
                             mets_dict):      
        """
        savely copy/create reaction from set of reaction parameter

        Parameters
        ----------
        rxn_id : str
            reaction ID.
        rxn_name : str
            DESCRIPreaction nameTION.
        lower_bound : float
            lower flux bound of reaction.
        upper_bound : float
            upper flux bound of reaction.
        mets_dict : dict
            dictionary of metabolites as returned by <model.reactions.XXX.metabolites>.

        Returns
        -------
        rxn_to_add : cobra.Reaction
            reaction object.

        """

        rxn_to_add = cobra.Reaction(rxn_id,
                                    name=rxn_name,
                                    lower_bound=lower_bound,
                                    upper_bound=upper_bound)
    
        # add metabolites
        rxn_to_add.add_metabolites(mets_dict)
        return rxn_to_add       
        
            
            
    def _init_model(self, model, copy_model, save_model, target_reaction):
        super().__init__(model,
                 copy_model=copy_model,
                 save_model=save_model,
                 target_reaction=target_reaction)
    

     
    @property
    def design(self):
        return self._design
    
    @design.setter
    def design(self, design):
        self.load_strain_design(design)
