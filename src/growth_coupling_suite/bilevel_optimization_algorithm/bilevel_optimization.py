"""Class for a bilevel model with an inner and outer objective function transformed
to a single-level optimization problem exploiting the dual theorem
"""

from __future__ import print_function, absolute_import

from growth_coupling_suite.bilevel_optimization_algorithm import lp_setup, lp_solving
import growth_coupling_suite as gcs

from os import mkdir
from os.path import isdir


class BilevelModel():
    
    def __init__(self, model, inner_objective, inner_objective_direction,
                 outer_objective, outer_objective_direction, deletion_targets=[], 
                 addin_targets=[], source_targets=[], cofeed_targets=[],
                 mediareduction_targets=[],
                 options={}):
        """
        

        Parameters
        ----------
        model : cobra.Model
            Metabolic model.
        inner_objective : str
            Reaction ID of the optimization target of the inner problem.
        inner_objective_direction : str
            Direction of the inner objective function (max: maximization, min: minimization).
        outer_objective : str
            Reaction ID of the optimization target of the outer problem..
        outer_objective_direction : str
            Direction of the outer objective function (max: maximization, min: minimization)..
        deletion_targets : list, optional
            Reaction IDs of deletion targets. The default is [].
        addin_targets : list, optional
            Reaction IDs of insertion targets. The default is [].
        source_targets : list, optional
            Reaction IDs of source targets. The default is [].    
        cofeed_targets : list, optional
            Reaction IDs of cofeed targets. The default is [].
        mediareduction_targets : lis, optional
            Reaction IDs of media deletion targets. The default is [].
        options : dict
            Optional inputs. The default is {}.
                n_total_interventions: Maximal number of total interventions
                n_deletion:         Maximal number of deletions
                n_addin:            Maximal number of reaction insertions
                n_cofeed:           Maximal number of cofeed metabolites
                n_mediareduction:   Maximal number of medium reductions
                n_carbonsource:     Maximal number of carbon sources
                
                time_limit:         time limit for solver [sec]
                n_threads:          number of workers for solving
                
                output_dir:         directory for saving outputs and results


        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
  
        # load parameters
        # mandatory input
        self.original_model = model.copy()
        self.inner_objective = inner_objective # list of variable identifiers
        self.outer_objective = outer_objective # list of variable identifiers
        # maximization: "max" or 1, minimization: "min" or -1
        self.inner_objective_direction = inner_objective_direction 
        self.outer_objective_direction = outer_objective_direction
        # optional inputs
        self.target_rxns = {"deletion": deletion_targets,
                            "addin": addin_targets,
                            "mediareduction": mediareduction_targets,
                            "cofeed": cofeed_targets,
                            "source": source_targets}
        
        # load options
        self.options = {}

        # intervention numbers
        # also consider intervention types in "lp_solving.py" and "lp_setup.py"
        available_interventions = ["n_total_interventions", # number of total interventions
                             "n_deletion", # number of deletions
                             "n_addin", # number of reaction addins
                             "n_cofeed", # number of cofed metabolites
                             "n_source", # number of source metabolites
                             "n_mediareduction"] # number of media reductions]
        
        intervention_numbers = {}
        for kind in available_interventions:
            if kind in options: intervention_numbers[kind] = options[kind]
        # total number of interventions must be specified
        if "n_total_interventions" not in intervention_numbers:
            raise Exception("No total intervention set size specified!")
            
        self.options["intervention_numbers"] = intervention_numbers
        
    
        # output directory and files
        if "output_dir" not in options.keys():
            self.options["output_dir"] = gcs.__path__[0] + "/results"
        else:
            self.options["output_dir"] = options["output_dir"]
            
        if not(isdir(self.options["output_dir"])):
            mkdir(self.options["output_dir"])
            
        print("Output directory: " + self.options["output_dir"])
         
        if "output_file_id" not in options.keys():
            self.options["output_file_id"] = ""
        else:
            self.options["output_file_id"] = options["output_file_id"]
        
        
        
        # model specifics
        # medium specifics of original model
        if "medium" in options:
            self.options["medium"] = options["medium"]
        else:
            self.options["medium"] = {}
            
        # consider GPR relations 
        if "eval_gpr" in options:
            self.options["eval_gpr"] = options["eval_gpr"]
        else:
            self.options["eval_gpr"] = False
            
        # evaluationof the growth phenotype of every solution
        if "eval_max_growth" in options:
            self.options["eval_max_growth"] = options["eval_max_growth"]
        else:
            self.options["eval_max_growth"] = 0    
            
        # solver options
        # solver time limit
        if "time_limit" in options:
            self.options.update(time_limit=options["time_limit"])
        else:
            self.options["time_limit"] = 400000
        # number of threads  
        if "n_threads" in options:
            self.options.update(n_threads=options["n_threads"])
        else:
            self.options["n_threads"] = 4
        # node memory limit for gurobi solver (gb)
        if "node_memory_limit" in options:
            self.options["node_memory_limit"] = options["node_memory_limit"]
        else:
            self.options["node_memory_limit"] = None
    
 
          
        # initialize bilevel model
        self.create_bilevel_model()
                
                
     
    # define properties
    @property
    def output_directory(self):
        return self.options["output_dir"]
    
    @output_directory.setter
    def output_directory(self, value):
        if isinstance(value, str):
            self.options["output_dir"] = value
            if not(isdir(self.options["output_dir"])):
                mkdir(self.options["output_dir"])
        
     
        
    def create_bilevel_model(self):
        """
        Create a bilevel model including booleans for target reaction decision
        variables and transform it into a mixed integer linear program (MILP)

        Returns
        -------
        None.

        """
        
        #  prepare metabolic model
        # load original metabolic model
        self.bilevel_model = self.original_model.copy()
    
        
        
        # create dual problem and combine with primal
        print("Dualize primal problem ...")
        dual_problem \
            = lp_setup.combine_dual_and_primal_problem(self.bilevel_model,
                                                 self.inner_objective,
                                                 self.inner_objective_direction)
        
        # add integer variables
        print("Add decision variables ...")
        y_vars = {}
        y_vars_kind = {}
        constrained_dual_vars = []
        for kind, rxn_list in self.target_rxns.items():
            for rxn_id in rxn_list:
                    y_var, constrained_vars \
                        = lp_setup.add_decision_variable(self.bilevel_model,
                                                         rxn_id, kind=kind)
                    y_vars[y_var] = rxn_id
                    y_vars_kind[y_var] = kind
                    constrained_dual_vars.extend(constrained_vars)
                    
        # link decision variables according to GPR relations
        if self.options["eval_gpr"]:
            print("\tLink deletion decision variables ...")
            
            gpr_linked_y_vars = {} # dict of decision variables to be linked
            y_var_count = {} # decision variables limited by nnumber of interventions constraint
            y_vars_kind_count = {} # kind of countable decision variables
            for y_var, kind in y_vars_kind.items():
                # only consider deletion targets
                if kind != "deletion":
                    # save reaction specific decision variable
                    y_var_count[y_var] = y_vars[y_var]
                    y_vars_kind_count[y_var] = y_vars_kind[y_var]
                    continue
                
                # determine gene-reaction rule
                gpr_rule = self.original_model.reactions.get_by_id(y_vars[y_var]).gene_reaction_rule
                # if GPR is not available continue
                if len(gpr_rule) == 0:
                    # save reaction specific decision variable
                    y_var_count[y_var] = y_vars[y_var]
                    y_vars_kind_count[y_var] = y_vars_kind[y_var]
                    continue
                
                # get set of genes
                gene_set = set(self.original_model.reactions.get_by_id(y_vars[y_var]).genes)
                
                gpr_rule_exist = None
                for gpr_rule_item, vars_genes in gpr_linked_y_vars.items():
                    if gene_set == vars_genes["gene_set"]:
                        # reserve addition decision variable to reaction cluster
                        gpr_rule_exist = gpr_rule_item
                        break
                    
                        
                # add decision variable to reaction cluster
                if gpr_rule_exist:
           
                    gpr_linked_y_vars[gpr_rule_exist]["y_vars"].append(y_var)
                    # print(gpr_rule_exist)
                    # print(gpr_linked_y_vars[gpr_rule_exist]["y_vars"])
                else:
         
                    # create new entry
                    gpr_linked_y_vars[gpr_rule] = {
                        "y_vars": [y_var],
                        "gene_set": gene_set
                        }
                           
                    
            # link dependent reaction cluster with an additional decision variable
            for y_var_cluster in gpr_linked_y_vars.values():
                # exclude clusters with only one reaction
                if len(y_var_cluster["y_vars"]) == 1:
                    y_var_count[y_var_cluster["y_vars"][0]] = y_vars[y_var_cluster["y_vars"][0]]
                    y_vars_kind_count[y_var_cluster["y_vars"][0]] = y_vars_kind[y_var_cluster["y_vars"][0]]
                    
                else:               
                    # link decision variables
                    yl_var = lp_setup.link_decision_variables(self.bilevel_model, y_var_cluster["y_vars"])
                    y_var_count[yl_var] = yl_var.name 
                    y_vars_kind_count[yl_var] = "deletion"
            
            
        else:
            # use reaction-specific decision variables as countable interventions
            y_var_count = y_vars
            y_vars_kind_count = y_vars_kind
                               
        
        # add inner optimality constraint (strong duality theory)
        # -> at optimality the dual and primal objective function values are equal
        print("Add inner optimality constraint (strong duality) ...")
        lp_setup.add_strong_duality_constraint(self.bilevel_model, dual_problem,
                                         constrained_dual_vars)
        
        # add constraints for the number of interventions
        print("Add intervention number constraints ...")
        lp_setup.add_intervention_number_constraints(self.bilevel_model,
                                               self.options["intervention_numbers"],
                                               y_var_count, y_vars_kind_count)
        
        # add outer objective
        print("Add outer objective function ...")
        lp_setup.add_objective(self.bilevel_model, self.outer_objective,
                         self.outer_objective_direction)
        
        
    
    
    
    def optimize(self, base_model=None):
        """Optimize the MILP problem of the transformed bilevel model
        

        Parameters
        ----------
        eval_gpr : bool, optional
            Evaluate GPR relations. The default is False.
        base_model : cobra.model, optional
            Base COBRA model on which the problem was built. The default is None.
        callback_parameters : dict, optional
            Parameters needed for gurobi solution callback . The default is {}.

        Returns
        -------
        callback_solutions : dict
            Solutions found for the MIPL problem.
            keys: objective function value at solution
            values: lists reaction that are deleted, added etc in solution
        solution_status : str
            Final status of the solver.

        """

        
        # MILP options
        opt = self.options
        
        # callback parameters
        callback_parameters = {
            "medium": opt["medium"],
            "eval_gpr": opt["eval_gpr"],
            "eval_max_growth": opt["eval_max_growth"]
            
            }
        
        # optimize MILP of the bilevel model
        self.callback_solutions, self.solution_status \
            = lp_solving.run_optimization(self.bilevel_model, self.outer_objective,
                     opt["output_dir"],
                     output_file_id=opt["output_file_id"],
                     time_limit=opt["time_limit"],
                     callback=True,
                     callback_parameters=callback_parameters,
                     num_threads=opt["n_threads"],
                     node_memory_limit=opt["node_memory_limit"],
                     base_model=base_model
                     )
    
        return self.callback_solutions, self.solution_status
    