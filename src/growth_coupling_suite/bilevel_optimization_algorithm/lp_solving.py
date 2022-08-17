"""
Set of methods for solving bilevel linear programs derived from 
stoichiometric metabolic models in a COBRA format
"""
from __future__ import print_function, division, absolute_import

from growth_coupling_suite.util.util import GPR_linked_reaction_knockouts

import cobra
import json
from os import mkdir
from os.path import exists, isfile
from shutil import rmtree
from pathlib import Path
# import pandas as pd
import numpy as np

try:
    from gurobipy import GRB
    import gurobipy as gp
except:
    print("Gurobi not installed")

# %%

def run_optimization(milp_model, target_rxn, output_dir, output_file_id="",
                     time_limit=400000, callback=True, callback_parameters={},
                     num_threads=4, node_memory_limit=None, base_model=None):
    """
    

    Parameters
    ----------
    milp_model : cobra.core.model
        single-level representation of bilevel problem.
    target_rxn : str
        engineering objective 
    output_dir : str
        directory for saving outputs.
    time_limit : int, optional
        time limit for solver in seconds. The default is 400000.
    callback : boolean, optional
        denotes if callback solutions from gurobi are saved. The default is True.
    num_threads : int, optional
        number of processes. The default is 4.
    node_memory_limit : int, optional
        memory limit of a single node in GB. The default is 1.
    base_model : cobra.Model
        Original metabolic model used to set up the MILP problem. The default is None.

    Returns
    -------
    solutions : dict
        returns all callback solutions from solver.
            keys: objective function value at solution
            values: lists reaction that are deleted, added etc in solution
    status : TYPE
        DESCRIPTION.

    """
    
    
    global output_file_format_for_callback
    output_file_format_for_callback = "callback_{0}_{1}_{2}.{3}"

    lp =  milp_model.solver.problem
    milp_model.solver.problem.setParam("FeasibilityTol", 1e-9)
    # Tolerate the smallest possible violation of integer bounds
    milp_model.solver.problem.setParam("IntFeasTol", 1e-9)
    milp_model.solver.problem.setParam("OptimalityTol", 1e-9)
    # Markowitz tolerance, default: 0.0078125
    # limits numerical issues the larger the value
    # milp_model.solver.problem.setParam("MarkowitzTol", 0.999)

    # enable lazy constraints
    if callback_parameters["eval_gpr"] or callback_parameters["eval_max_growth"]:
        milp_model.solver.problem.setParam("LazyConstraints", 1)
    else:
        milp_model.solver.problem.setParam("LazyConstraints", 0)

    # make output more verbose and allow progress to be tracked
    milp_model.solver.problem.setParam("OutputFlag", 1)
    if time_limit:
        milp_model.solver.problem.setParam("TimeLimit", time_limit)

    # objective function cutoff
    # milp_model.solver.problem.setParam("Cutoff", 0.001)
    milp_model.solver.problem.setParam("Presolve", 2)
    
    # save solver output in log files
    if not exists("Gurobi log files/"):
        mkdir("Gurobi log files/")
    milp_model.solver.problem.setParam("LogFile", "Gurobi log files/{}_log".format(target_rxn))

    if num_threads:
        milp_model.solver.problem.setParam("Threads", num_threads)
    if node_memory_limit:
        # Save memory by saving nodefiles to disk (number in gigB)
        milp_model.solver.problem.setParam("NodefileStart", node_memory_limit)
        nodeFiles_path = output_dir + "/" + "NodeFiles_" + target_rxn + "_" + output_file_id
        if not exists(nodeFiles_path):
            mkdir(nodeFiles_path)
        milp_model.solver.problem.setParam("NodefileDir", nodeFiles_path)

    
    # create callback parameters
    if "medium" not in callback_parameters:       
        callback_parameters["medium"] = {}
    if "eval_gpr" not in callback_parameters:
        callback_parameters["eval_gpr"] = False
    if "eval_max_growth" not in callback_parameters:
        callback_parameters["eval_max_growth"] = False
    

    print("Start solving MILP ...")
    solutions, status = solve_milp(lp,
                                   milp_model,
                                   target_rxn,
                                   output_dir,
                                   output_file_id,
                                   callback=callback,
                                   callback_parameters=callback_parameters,
                                   base_model=base_model)
    
    # delete nodeFiles from disk
    if node_memory_limit:
        if exists(nodeFiles_path):
            print("Remove NodeFiles created by the solver...")
            try:
                rmtree(nodeFiles_path)
                print("\tCompleted")
            except:
                print("\tFailed")
                
    
    return solutions, status
    


def evaluate_callback_solution(lp, cobra_model, target_rxn, parameters,
                               base_model):
    """
    Evaluate a new solution from the solver

    Parameters
    ----------
    lp : TYPE
        gurobi model of the MILP problem.
    cobra_model : cobra.core.model
        metabolic model in cobra format.
    target_rxn : str
        target reaction of the objective function.
    parameters : dict
        parameters and options for the callback.

    Returns
    -------
    solution_is_valid : boolean
        Denotes if a returned solution is numerically and metabolically valid.
    x_dict : dict
        design solution in dictionary.
    out_rxns : list
        reaction IDs targeted by solution.
    out_solution : dict
        details of the new solution.

    """
        
    
    # get solution details from solver
    x_dict = {r.VarName: value for r, value in zip(lp._vars, lp.cbGetSolution(lp._vars))}
    objective = lp.cbGet(GRB.Callback.MIPSOL_OBJ)
    time = lp.cbGet(GRB.Callback.RUNTIME)
    
    # only consider solutions above an objective value threshold
    solution_is_valid = True
    out_rxns = []
    all_added_grp_linked_reactions = []
    interventions = {} # interventions including GPR dependencies
    interventions_solver = {} # interventions returned by the solver
    if abs(float(objective)) > 1e-5:
        
        print("\n> Evaluate callback solution ...")
        # save uptake rates of cofeeds and sources
        cofeed_fluxes = {}
        source_fluxes = {}
        
        # evaluate feasibility of design solution according to GRP relations
        for variable, value in x_dict.items():
            # determine interger variables which enforce a knockout
            if variable.startswith("y_") and value > .99:
                r_id = variable[2:] # determine ID of deleted reaction
                
                
                if r_id is not None:
                    # save reaction ID
                    out_rxns.append(r_id)
                    
                    # get and save ID and target type, and fluxes of cofeed/source targets
                    available_types = ["_deletion", "_addin", "_source", "_cofeed",
                                       "_mediareduction"]
                    for kind in available_types:
                        if r_id.endswith(kind):
                            
                            if kind in ["_source", "_cofeed"]:
                                # save flux rates
                                # get reverse flux
                                for v, val in x_dict.items():
                                    if v.startswith(r_id.replace(kind, "_reverse")):
                                        r_id_reverse = v
                                        break
                                    
                                if kind == "_cofeed":
                                    cofeed_fluxes[r_id_reverse] = x_dict[r_id_reverse]
                                elif kind == "_source":
                                    source_fluxes[r_id_reverse] = x_dict[r_id_reverse]
                                    
                                # print(r_id, x_dict[r_id.replace(kind, "")], x_dict[r_id_reverse])
                            
                            intervention_dict = {"ID": r_id[0:r_id.rindex(kind)],
                                                 "var_name": variable,
                                                 "type": kind[1:]
                                                 }
                            # save reaction bounds
                            if kind in ["_deletion", "_mediareduction"]:
                                intervention_dict["lower_bound"] = 0
                                intervention_dict["upper_bound"] = 0
                            else:
                                rxn = cobra_model.reactions.get_by_id(intervention_dict["ID"])
                                intervention_dict["lower_bound"] = rxn.lower_bound
                                intervention_dict["upper_bound"] = rxn.upper_bound
                            
                            break
                    # save reaction information
                    interventions[intervention_dict["ID"]] = intervention_dict
                    interventions_solver[intervention_dict["ID"]] = intervention_dict

                    
                    
                    # determine GPR-linked reactions for deletion targets 
                    # GPR-linked reactions are only added to the set of interventions
                    # evaluation of extended design solution is done after optimization is finished  
                    if parameters["eval_gpr"]:
                    
                        
                        
                        
                        kind = "_deletion"
                        if r_id.endswith(kind):
                            
                            # get reaction ID in cobra mdoel
                            r_id_model = r_id[0:r_id.rindex(kind)]
                            
                            # check if targeted reaction has already been evaluated 
                            # and if lazy constraints has been added
                            if r_id in parameters["gpr_evaluated_objects"]:
                                # print("\tReload GPR-linked reactions for gene-based knockout of ", r_id_model)
                                # load specifics and continue with next intervention
                                interventions = {**interventions,
                                    **parameters["gpr_evaluated_objects"][r_id]["interventions"]
                                    }
                                all_added_grp_linked_reactions.extend(parameters["gpr_evaluated_objects"][r_id]["all_added_grp_linked_reactions"])
                                solution_is_valid = solution_is_valid & parameters["gpr_evaluated_objects"][r_id]["solution_is_valid"]
                                continue
                                
                            
                            # determine reactions also affected by gene-based knockout of r_id
                            dep_rxns = GPR_linked_reaction_knockouts(cobra_model, [r_id_model], eval_gpr=True)
                
                            # create additional lazy constraint to enforce simultaneous deletion of reactions
                            # add additional dependent reaction deletions
                            if len(dep_rxns) > 0:
                                print("\tFound GPR-linked reactions for gene-based knockout of", r_id_model)
                                
                                # allocate dict to store linked reaction if design object occurs in future solutions
                                parameters["gpr_evaluated_objects"][r_id] = {
                                    "interventions": {},
                                    "all_added_grp_linked_reactions": [],
                                    "all_grp_linked_reactions!": [],
                                    "solution_is_valid": True
                                    }
                                
                                # allocate lazy constraint
                                lazy_cnstr_lhs_list = []
                                gpr_linked_rxns = [r_id_model]
                                for var in lp._vars:
                                    if var.VarName == variable:
                                        lazy_cnstr_lhs_list.append(var)
                                        break
                                    
                                # get decision variables for dependent reactions from solver
                                # add reaction knockout to set of interventions
                                for var in lp._vars:
                                    var_name_model = var.VarName.replace("y_", "").replace(kind, "")
                                    if var.VarName.startswith("y_") and (var_name_model in dep_rxns):
                                        # save (decision) variable
                                        lazy_cnstr_lhs_list.append(var)
                                        gpr_linked_rxns.append(var_name_model)
                                        all_added_grp_linked_reactions.append(var_name_model)
                                                                                
                                        # add reaction to set of interventions
                                        interventions[var_name_model] = {
                                            "ID": var_name_model,
                                            "var_name": var.VarName,
                                            "type": "deletion",
                                            "lower_bound": 0,
                                            "upper_bound": 0
                                            }
                                        
                                        # save if design object occurs in future solutions
                                        parameters["gpr_evaluated_objects"][r_id]["interventions"][var_name_model] = interventions[var_name_model]
                                        parameters["gpr_evaluated_objects"][r_id]["all_added_grp_linked_reactions"].append(var_name_model)
                                
                                parameters["gpr_evaluated_objects"][r_id]["all_grp_linked_reactions"] = gpr_linked_rxns
                                        
                                # check if lazy constraints for the same set of dependent reactions have already been integrated
                                # current_linked_rxns_set = parameters["gpr_evaluated_objects"][r_id]["all_added_grp_linked_reactions"].copy()
                                # current_linked_rxns_set.append(r_id)
                                current_linked_rxns_set = set(gpr_linked_rxns)
                                    
                                
                                lazy_cnstr_duplicate = False
                                for curr_r_id, eval_obj in parameters["gpr_evaluated_objects"].copy().items():
                                    # skip current reaction seet
                                    if curr_r_id == r_id: continue
                                
                                    # saved_linked_rxns_set = eval_obj["all_added_grp_linked_reactions"].copy()
                                    # saved_linked_rxns_set.append(curr_r_id)
                                    saved_linked_rxns_set = set(eval_obj["all_grp_linked_reactions"])
                                        
                                    if saved_linked_rxns_set == current_linked_rxns_set:
                                        # lazy constraints for reaction set already integrated
                                        lazy_cnstr_duplicate = True
                                        solution_is_valid = eval_obj["solution_is_valid"]
                                        parameters["gpr_evaluated_objects"][curr_r_id]["solution_is_valid"] = solution_is_valid
                                        break
                                        
                                           
                                
                                # if there is no decision variable for any dependent reaction dismiss reaction target       
                                if lazy_cnstr_duplicate:
                                    # do not add duplicate lazy constraint
                                    print("\t\tLazy constraints for linked reaction set already existent")
                                                                    
                                
                                elif len(lazy_cnstr_lhs_list) != (len(dep_rxns)+1):
                                    # disable reaction knockout, knockout would target a protected reaction
                                    print("\t\tGPR-linked reactions are protected:", ", ".join(list(set(dep_rxns)-set(gpr_linked_rxns))))
                                    lp.cbLazy(lazy_cnstr_lhs_list[0] <= .01)
                                    solution_is_valid = False
                                    parameters["gpr_evaluated_objects"][r_id]["solution_is_valid"] = False
                                else:
                                    # set lazy constraint, either all or none of the dependent reactions should be knocked out
                                    # decision variable values for all dependent reactions must be the same
                                    # set pairwise equalities
                                    print("\t\tLink", ", ".join(gpr_linked_rxns), "with lazy constraints")
                                    for var in lazy_cnstr_lhs_list[1:]:                                  
                                        lp.cbLazy(lazy_cnstr_lhs_list[0] - var <= .02)
                                        lp.cbLazy(lazy_cnstr_lhs_list[0] - var >= -.02)
                            
        
                 
        
        
        # compute objective functions with original model -> phenotype
        if parameters["eval_max_growth"] > 0:
            with base_model as model:
                # apply interventions to model (including GPR dependencies)
                # design_reactions = [r for r in interventions.keys()]
                for rxn_id, intervention in interventions.items():
                    rxn = model.reactions.get_by_id(rxn_id)
                    # set bounds
                    rxn.lower_bound = intervention['lower_bound']
                    rxn.upper_bound = intervention['upper_bound']
                    
                # compute optimal objective function value
                sol_objective = model.slim_optimize()
                
            # Is strain design feasible in the original COBRA model?
            if sol_objective < parameters["eval_max_growth"]:
                # strain design is not feasible, set lazy constraints to prohibit solution
                solution_is_valid = False
                # get variables of interventions (only interventions from the solver)
                i_var_names = [i['var_name'] for rxn_id, i in interventions_solver.items()]
                lazy_cnstr_lhs_list = [
                    var
                    for var in lp._vars
                    if var.VarName in i_var_names
                    ]
                
                # to prohibit solution, sum of intervention decision variable values
                # has to be below the number of interventions in the solution
                lp.cbLazy(gp.quicksum(var for var in lazy_cnstr_lhs_list) <= len(lazy_cnstr_lhs_list)*.99)
                
                print('\n\t Solution violates phenotypic constraints')
                print("\t\tLink", ", ".join([r for r in interventions_solver.keys()]), "with lazy constraints")
                     
        
        
    else:
        # invalid solution
        print('> Solution dismissed. Minimum objective function constraint violated.')
        solution_is_valid = False
        
        
    # summarize solution output   
    out_solution = {"interventions": interventions,
                    "GPR_linked_interventions_added": list(np.unique(all_added_grp_linked_reactions)),
                    "objective_value": objective,
                    "time": time,
                    "node_count_explored": lp.cbGet(GRB.Callback.MIPSOL_NODCNT),
                    "best_bound": lp.cbGet(GRB.Callback.MIPSOL_OBJBND),
                    "medium": parameters["medium"]
                    }
        
        
    return solution_is_valid, x_dict, out_rxns, out_solution
        
        


def save_callback_solution(lp, cobra_model, target_rxn,
                           out_solution, out_rxns, is_valid):
    """
    

    Parameters
    ----------
    lp : TYPE
        gurobi model of the MILP problem..
    cobra_model : cobra.Model
        metabolic model in cobra format..
    target_rxn : str
        target reaction of the objective function..
    out_solution : dict
        interventions of the solution.
    out_rxns : list
        IDs of reactions being part of the solution.
    is_valid : bool
        Denotes if a solution is valid.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    global return_callback_solutions
    global output_dir_for_callback
    global output_dir_for_dismissed_solutions
    global output_text_file_for_callback
    global output_file_for_callback


    
    print("... saving callback solution", '[{}]'.format(
        ','.join(rxn_id for rxn_id in out_solution['interventions'].keys())
        ))

    
    if is_valid:
        # save solution to valid solutions folder
        output_path_for_callback = str(Path(output_dir_for_callback+'/'+output_file_for_callback))
        output_text_path_for_callback = str(Path(output_dir_for_callback+'/'+output_text_file_for_callback))
        
    else:
        # save solution to dismissed solutions folder
        output_path_for_callback = str(Path(output_dir_for_dismissed_solutions+'/'+output_file_for_callback))
        output_text_path_for_callback = str(Path(output_dir_for_dismissed_solutions+'/'+output_text_file_for_callback))

    # load output file
    if isfile(output_path_for_callback):
        with open(output_path_for_callback, "r") as out_file:
                out_rxns_dict = json.load(out_file)
    else:
        out_rxns_dict = {}
                    
    out_rxns_dict["sol_" + str(len(out_rxns_dict)+1)] = out_solution

    # save and close file
    with open(output_path_for_callback, "w") as out_file:
        json.dump(out_rxns_dict, out_file)
        
    # only return solution if valid
    if is_valid:
        return_callback_solutions = out_rxns_dict  
    
    
    # save solution in text format
    out_array = [out_solution["objective_value"],
                 out_solution["time"],
                 "|".join(out_rxns)]
    

    if not isfile(output_text_path_for_callback):
        with open(output_text_path_for_callback, "a") as f:
            f.write(str(["objective", "time", "reactions"]) + "\n")
    with open(output_text_path_for_callback, "a") as f:
        f.write(str(out_array) + "\n")




# def return_milp_solution_df(reaction_changes_from_mip, save_loc):
#     # Output alternative solutions to csv
#     out_df = pd.DataFrame()
#     for iter_num in reaction_changes_from_mip:
#         i = 0
#         for label in reaction_changes_from_mip[iter_num]:
#             if type(label) is float:
#                 name = "Objective Value"
#             else:
#                 name = "Change_%i" % i
#             out_df.loc[iter_num, name] = label
#             i += 1
#     out_df.to_csv(save_loc)


def solve_milp(lp, milp_model, target_rxn, output_dir, output_file_id,
               callback=True, callback_parameters={}, base_model=None):
    """A performance tunable method for updating a model problem file
    """
    global target_for_callback
    global model_for_callback
    global output_dir_for_callback
    global output_dir_for_dismissed_solutions
    global output_file_for_callback
    global output_text_file_for_callback
    global return_callback_solutions
    global output_file_id_for_callback
    global parameters_for_callback
    global original_model_for_callback

    model_for_callback = milp_model
    target_for_callback = target_rxn
    output_dir_for_callback = output_dir
    return_callback_solutions = []
    output_file_id_for_callback = output_file_id
    parameters_for_callback = callback_parameters
    parameters_for_callback["gpr_evaluated_objects"] = {} # save design objects for which GPR-linked reactions have already been evaluated
    original_model_for_callback = base_model
    
    # save milp model
    cobra.io.save_json_model(milp_model, output_dir + "/" 
                + output_file_format_for_callback.format("model", target_rxn,
                                                         output_file_id, "json"))
    
    # create dictionary and text file for storing solutions
    output_file_for_callback = output_file_format_for_callback.format(
        "solutions_dict",
        target_rxn,
        output_file_id,
        "json")
        
    output_text_file_for_callback = output_file_format_for_callback.format(
        "solutions",
        target_rxn,
        output_file_id,
        "txt")    
    
    # create folder for dismissed solutions
    output_dir_for_dismissed_solutions = output_dir+"/dismissed_solutions/"
    if not exists(output_dir_for_dismissed_solutions):
        mkdir(output_dir_for_dismissed_solutions)
        

    lp.update()
    if callback:
        lp._vars = lp.getVars()
        lp.optimize(callback=callback_function)
    else:
        lp.optimize()
    status = milp_model.solver.status
    
    return return_callback_solutions, status




def callback_function(lp, where):
    """
    return and save an incumbent MILP solution when a new one is found

    Parameters
    ----------
    lp : TYPE
        gurobi model of the MILP problem...
    where : str
        Indicates from whre in the Gurobi solver the callback is called.

    Returns
    -------
    None.

    """

    if where == GRB.Callback.MIPSOL:
        
        # evaluate new solution
        solution_is_valid, x_dict, out_rxns, out_solution \
            = evaluate_callback_solution(lp,
                                           model_for_callback,
                                           target_for_callback,
                                           parameters_for_callback,
                                           original_model_for_callback)
        
        save_callback_solution(
            lp,
            model_for_callback,
            target_for_callback,
            out_solution,
            out_rxns,
            solution_is_valid
            )
  