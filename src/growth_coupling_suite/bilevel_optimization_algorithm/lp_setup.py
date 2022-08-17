"""
Set of methods for setting up bilevel linear programs derived from 
stoichiometric metabolic models in a COBRA format
"""
from __future__ import print_function, absolute_import
from optlang.duality import convert_linear_problem_to_dual
import optlang

    
def combine_dual_and_primal_problem(model, inner_objective, inner_objective_direction):
    """
    

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.
    inner_objective : str
        Reaction ID of the optimization target of the inner problem.
    inner_objective_direction : str
        Direction of the inner objective function (max: maximization, min: minimization).

    Raises
    ------
    NameError
        DESCRIPTION.

    Returns
    -------
    dual_problem : cobra.Model
        Model combining the primal (metabolic) model and its dual.

    """
  
    # Get interface to use for class methods, etc
    interface = model.solver.interface
    
    # first set objective to minimization for inner optimality constraint
    model.objective = model.problem.Objective(0) # reset objective to default
    # determine objective coefficient
    if (inner_objective_direction == "max") or (inner_objective_direction == 1):
        objective_coefficient = 1  
        info_str = "\tInner objective function maximizes {0}"
    elif (inner_objective_direction == "min") or (inner_objective_direction == -1):
        objective_coefficient = -1  
        info_str = "\tInner objective function minimizes {0}"
    else:
        raise NameError("Unknown inner objective direction: " + str(inner_objective_direction))
        
    # set inner objective                 
    for r in model.reactions:
        if r.id == inner_objective:
            r.objective_coefficient = objective_coefficient
            print(info_str.format(inner_objective))
            
    # copy and test model
    copied_model = model.copy()
    copied_model.optimize()
    
    # get primal + dual and combine
    dual_problem = convert_linear_problem_to_dual(copied_model.solver)
    primal_problem = model.solver
    for var in dual_problem.variables:  # All variables in the dual are
        # copied to the primal
        var = interface.Variable.clone(var)
        primal_problem.add(var)
    for const in dual_problem.constraints:  # All constraints in the dual
        # are copied to the primal
        const = interface.Constraint.clone(const, model=primal_problem)
        primal_problem.add(const)

    dual_problem.optimize()
    
    return dual_problem


def add_decision_variable(model, rxn_id, kind="deletion"):
    """
    Create a boolean variable to control the target reaction and add it to the model

    Parameters
    ----------
    model : cobra.Model
        Metabolicmodel.
    rxn_id : str
        ID of the target reaction.
    kind : str, optional
        mutation type of the target. The default is "deletion".

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    y_var : cobra.Model.solver.interface.Variable
        boolean model variable.
    constrained_vars : list
        variables of the dual constraints.

    """

    # integrate binary variable
    interface = model.solver.interface
    reaction = model.reactions.get_by_id(rxn_id)
    y_var = interface.Variable("y_" + reaction.id + "_"+ kind, type="binary")
    y_var.decision_reaction_id = rxn_id
    
    # define variable expression
    if kind in ["deletion", "mediareduction"]:
        # for y_var = 1 reaction is deleted
        y_var_expr = 1 - y_var
        y_var_dual_exp = y_var
    elif kind in ["addin", "source", "cofeed"]:
        # for y_var = 1 reaction is activated/added
        y_var_expr = y_var
        y_var_dual_exp = 1 - y_var
    else:
        raise Exception("Kind {} not valid".format(kind))
        
    
    # add constraints
    # v <= ub * y  -->  v - ub * y <= 0
    model.solver.add(
        interface.Constraint(reaction.flux_expression -
                             reaction.upper_bound * y_var_expr, ub=0,
                             name="primal_y_const_" + reaction.id + "_ub"))
    # v >= lb * y  -->  v - lb * y >= 0
    model.solver.add(
        interface.Constraint(reaction.flux_expression -
                             reaction.lower_bound * y_var_expr, lb=0,
                             name="primal_y_const_" + reaction.id + "_lb"))    
    
    constrained_vars = []

    if reaction.upper_bound != 0:
        dual_forward_ub = model.solver.variables[
            "dual_" + reaction.forward_variable.name + "_ub"]
        model.solver.add(
            interface.Constraint(dual_forward_ub - 1000 * (y_var_dual_exp),
                                 ub=0, name="dual_y_const_" + reaction.forward_variable.name + "_ub"))
        constrained_vars.append(dual_forward_ub)
    if reaction.lower_bound != 0:
        dual_reverse_ub = model.solver.variables[
            "dual_" + reaction.reverse_variable.name + "_ub"]
        model.solver.add(
            interface.Constraint(dual_reverse_ub - 1000 * (y_var_dual_exp),
                                 ub=0, name="dual_y_const_" + reaction.forward_variable.name + "_lb"))
        constrained_vars.append(dual_reverse_ub)

    return y_var, constrained_vars


def link_decision_variables(model, y_vars):
    """
    link decision variables by pairwise equality constraints with new decision variable,
    e.g., to account for reactions catalyzed by the same enzyme

    Parameters
    ----------
    model : cobra.Model
        DESCRIPTION.
    y_vars : list
        boolean decision variables to be linked.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    yl_var : interface.Variable
        boolean decision variable linking y_vars.

    """
    
    # create new linking, binary decision variable
    # continue numbering
    yl_num = str(len([v.name for v in model.solver.variables if v.name.startswith("yl_")]))
    # create variable
    interface = model.solver.interface
    yl_var = interface.Variable("yl_" + yl_num, type="binary")
    
    
    # set pairwise equality constraints for the parsed decision variables
    count = 0
    for y_var in y_vars:
        linking_equality_constraint = model.solver.interface.Constraint(
            y_var - yl_var, # y - yl = 0
            lb=0,
            ub=0,
            name="linking_equality_constraint_" + yl_num + "_" + str(count)
        )
        model.solver.add(linking_equality_constraint)
        count += 1
        
    
    return yl_var



def add_strong_duality_constraint(model, dual_problem, constrained_dual_vars):
    """
    set equality constraint between primal and dual objective

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.
    dual_problem : cobra.Model
        Dual model.
    constrained_dual_vars : list
        variables of the dual constraints.

    Returns
    -------
    None.

    """
    # get interface of the solver
    interface = model.solver.interface
    
    primal_objective = model.solver.objective
    dual_objective = interface.Objective.clone(dual_problem.objective,
                                               model=model.solver)
    
    # construct dual objective
    reduced_expression = optlang.symbolics.Add(
        *((c * v) for v, c in
          dual_objective.expression.as_coefficients_dict().items()
          if v not in constrained_dual_vars)
    )
    dual_objective = interface.Objective(reduced_expression,
                                         direction=dual_objective.direction)
    
    # set equality constraint between primal and dual objective
    optimality_constraint = model.solver.interface.Constraint(
        primal_objective.expression - dual_objective.expression,
        lb=0, ub=0, name="inner_optimality")
    model.solver.add(optimality_constraint)
    
    
def add_intervention_number_constraints(model, intervention_numbers,
                                        y_vars, y_vars_kind):
    """
    Add constraints to limit the number of interventions

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.
    intervention_numbers : dict
        maximum number of interventions for various intervention types (deletions, insertions etc).
    y_vars : list
        decision variables.
    y_vars_kind : list
        type of intervention a decision variable controls.

    Returns
    -------
    None.

    """
    
    # add total intervention number constraints
    total_intervention_constraint = model.solver.interface.Constraint(
        optlang.symbolics.Add(*y_vars), lb=0,
        ub=intervention_numbers["n_total_interventions"],
        name="number_of_total_interventions_constraint"
    )
    model.solver.add(total_intervention_constraint)
    # print information
    print("\tMaximum number of interventions: {0} out of {1}".format(
                intervention_numbers["n_total_interventions"],
                len(y_vars)))
    
    
    
    # sort intervention types
    y_var_sort = {}
    for y_var, kind in y_vars_kind.items():
        if kind in y_var_sort.keys():
            y_var_sort[kind][y_var] = y_vars[y_var]
        else:
            y_var_sort[kind] = {y_var: y_vars[y_var]}
            
    # add intervention kind specific constraints
    for kind in y_var_sort:
          
        # check if intervention number was specified
        intervention_key = "n_" + kind
        if intervention_key not in intervention_numbers:
            intervention_numbers[intervention_key] \
                = intervention_numbers["n_total_interventions"]
        # enforce choice of a source target
        if kind == "source":
            lower_bound = 1
        else:
            lower_bound = 0
        # add constraint              
        intervention_constraint = model.solver.interface.Constraint(
            optlang.symbolics.Add(*y_var_sort[kind]), lb=lower_bound,
            ub=intervention_numbers[intervention_key],
            name="number_of_" + kind + "s_constraint"
            )
        model.solver.add(intervention_constraint)
        
        # print information
        print("\tRange of {0} targets: {3}-{1} out of {2}".format(kind,
                intervention_numbers[intervention_key],
                len(y_var_sort[kind]), lower_bound))
        
        
 
def add_objective(model, outer_objective, outer_objective_direction):
    """
    Add the outer objective function

    Parameters
    ----------
    model : cobra.Model
        MILP model.
    outer_objective : str
        ID of the target reaction of the outer objective.
    outer_objective_direction : str
        Direction of the outer objective function (max, min).

    Returns
    -------
    None.

    """
    # add objective
    model.objective = outer_objective
    # change direction
    model.objective_direction = outer_objective_direction
    
    # print information
    print("\tOuter objective function {0}imizes {1}".format(
        outer_objective_direction, outer_objective))
    