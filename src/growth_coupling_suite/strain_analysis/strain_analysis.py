# Analyzes metabolic models for key phenotypic parameters

import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from warnings import warn


# Disable gurobi logging output
try:
    import gurobipy
    gurobipy.setParam("OutputFlag", 0)
except ImportError:
    pass


class StrainAnalyzer():
    
    def __init__(self, model, copy_model=True, save_model=True,
                 target_reaction=None):
        
        
        if copy_model:
            self.model = model.copy()
        else:
            self.model = model
            
        # save original model 
        if save_model:
            self._model_original = model.copy()
        else:
            self._model_original = None
                
        # create hidden data storage for methods
        self._flux_space = {}
        self._plot_2D = {}
        
        # initialize parameters
        self.target_reaction = target_reaction
        
        
    def maximum_theoretical_yield(self, target_rxn, substrate_rxn):
        # calculate maximum theoretical flux yield of a reaction on a specified substrate

        # check if substrate exchange reaction is open
        if self.model.reactions.get_by_id(substrate_rxn).lower_bound >= 0:
            raise Exception("Substrate uptake reaction " + substrate_rxn + "is blocked")
            
        with self.model as model_m:
            model_m.objective = target_rxn
            # optimize model
            sol = model_m.optimize()
            # calculate theoretical yield [mol/mol]
            return sol.objective_value / -sol.fluxes[substrate_rxn]
    
        
    def flux_at_maximum_objective(self, target_rxns=None):
        # calculate flux rate at maximum objective
        
        if isinstance(target_rxns, str):
            target_rxns = [target_rxns]           
        
        try:
            sol = cobra.flux_analysis.pfba(self.model)
        except:
            # infeasible models
            return {}
        
        # 
        if isinstance(target_rxns, list):
            # extract fluxes of specified reactions
            fluxes = {}
            for t in target_rxns:
                fluxes[t] = sol.fluxes[t]
                
        else:
            # return all fluxes
            fluxes = sol.fluxes
                            
        return fluxes
    
    
    def flux_space_projection(self, target_rxn=None, reference_rxn=None, grid_size=15,
                              yield_space=False, yield_reference_reaction=None,
                              unit_transformation_factors={},
                              plot=True, fill_flux_space=True, save_name=None,
                              xlabel=None, ylabel=None, **kwargs):
        # plots (2D) flux space projection
        # default reference reaction is the objective function (biomass equation)
        
        # check input
        if not(isinstance(target_rxn, str)):
            if not(isinstance(self.target_reaction, str)):
                raise Exception("Provide target reaction as a string")
            else:
                target_rxn = self.target_reaction
            
        # protect model
        # objective_save = self.model.objective
        # objective_direction_save = self.model.objective_direction
        # model_m = self.model

        # with self.model as model_m:    
        
        # Determine reference reaction from objective function
        if not(reference_rxn):
            # identifiy objective function
            objective_rxn = []
            for rxn in self.model.reactions:
                if rxn.objective_coefficient != 0:
                    objective_rxn.append(rxn.id)
            if len(objective_rxn) != 1:
                raise Exception("Objective function does not include a single reaction")
            else:
                reference_rxn = objective_rxn[0]
                
        else:
            # reference reaction provided
            if reference_rxn not in [rxn.id for rxn in self.model.reactions]:
                raise Exception("Reference reaction " + reference_rxn + " is not in model")
            
        # check target reaction
        if target_rxn not in [rxn.id for rxn in self.model.reactions]:
                raise Exception("Target reaction " + target_rxn + " is not in model")
        else:
            target_rxn_obj = self.model.reactions.get_by_id(target_rxn)

        # protect model and compute range of reference reaction flux
        with self.model as model_m:             
            # calculate maximum/minimum reference reaction flux
            model_m.objective = reference_rxn
            model_m.objective_direction = "max"
            ref_max = model_m.slim_optimize()
            model_m.objective_direction = "min"
            ref_min = model_m.slim_optimize()
            
        # check if reference reaction can be active
        if np.isnan(ref_max) or np.isnan(ref_min):
            warn("Reference reaction " + reference_rxn + " cannot carry any flux!", UserWarning)
            return None
        elif abs(ref_max) < 1e-4:
            warn("Reference reaction " + reference_rxn + " cannot carry any flux!", UserWarning)
            return None
        
        # scan flux space for target reaction
        reference_grid = np.linspace(ref_min, ref_max, grid_size)

        # protect model and compute flux space
        with self.model as model_m:  
        
            ref_rxn = model_m.reactions.get_by_id(reference_rxn)
            # ref_rxn_bounds_save = ref_rxn.bounds
            model_m.objective = target_rxn #  set objective function
            
            # define factors for unit transformation
            if target_rxn in unit_transformation_factors:
                target_unit_t = unit_transformation_factors[target_rxn]
            else:
                target_unit_t = 1
                
            if reference_rxn in unit_transformation_factors:
                ref_unit_t = unit_transformation_factors[reference_rxn]
            else:
                ref_unit_t = 1

            
            fs_lb = np.empty([0,1])
            fs_ub = np.empty([0,1])
            
            # distinguish between yield and flux space
            if yield_space:
                # calculate yield space from reaction yields
                # check yield reference reaction, normally the substrate uptake reaction
                if not(yield_reference_reaction) or (yield_reference_reaction not in model_m.reactions):
                    raise Exception("Yield reference reaction " + yield_reference_reaction + " is not available")
                else:
                    ref_yield_rxn = model_m.reactions.get_by_id(yield_reference_reaction)
                    # define unit transformation
                    if yield_reference_reaction in unit_transformation_factors:
                        yield_ref_unit_t = unit_transformation_factors[yield_reference_reaction]
                    else:
                        yield_ref_unit_t = 1
                
                
                # save target reaction bounds
                target_rxn_bounds = target_rxn_obj.bounds
                
                ref_yield_lb = np.empty([0,1])  
                ref_yield_ub = np.empty([0,1]) 
                
                for grid_point in reference_grid:
                    # fix reference flux
                    ref_rxn.bounds = (grid_point, grid_point)
                    # ref_rxn.lower_bound = grid_point
                    # ref_rxn.upper_bound = grid_point

                    # calculate upper flux space bound
                    model_m.objective = target_rxn #  set objective function
                    model_m.objective_direction = "max"
                    sol_t = model_m.slim_optimize(error_value="infeasible")
                    if sol_t == "infeasible":
                        continue
                    # calculate minimal corresponding yield reference flux
                    target_rxn_obj.bounds = (sol_t, sol_t)
                    # target_rxn_obj.lower_bound = sol_t
                    # target_rxn_obj.upper_bound = sol_t
                    model_m.objective = yield_reference_reaction
                    model_m.objective_direction = "min"
                    sol_r = model_m.slim_optimize()
                    # reset bounds
                    target_rxn_obj.bounds = target_rxn_bounds
                    # target_rxn_obj.lower_bound = target_rxn_bounds[0]
                    # target_rxn_obj.upper_bound = target_rxn_bounds[1]
                    # check solution
                    if abs(sol_r) < 1e-4:
                        # dismiss solution to avoid dividing by zero
                        continue
                    # save
                    fs_ub = np.append(fs_ub, (sol_t*target_unit_t) / (-sol_r*yield_ref_unit_t))
                    ref_yield_ub = np.append(ref_yield_ub, (grid_point*ref_unit_t) / (-sol_r*yield_ref_unit_t))
                    
                    # calculate lower flux space bound
                    model_m.objective = target_rxn #  set objective function
                    model_m.objective_direction = "min"      
                    sol_t = model_m.slim_optimize()
                    # calculate minimal corresponding yield reference flux
                    target_rxn_obj.bounds = (sol_t, sol_t)
                    # target_rxn_obj.lower_bound = sol_t
                    # target_rxn_obj.upper_bound = sol_t
                    model_m.objective = yield_reference_reaction
                    model_m.objective_direction = "min"
                    sol_r = model_m.slim_optimize()
                    # reset bounds
                    target_rxn_obj.bounds = target_rxn_bounds
                    # target_rxn_obj.lower_bound = target_rxn_bounds[0]
                    # target_rxn_obj.upper_bound = target_rxn_bounds[1]
                    # save
                    fs_lb = np.append(fs_lb, (sol_t*target_unit_t) / (-sol_r*yield_ref_unit_t))
                    ref_yield_lb = np.append(ref_yield_lb, (grid_point*ref_unit_t) / (-sol_r*yield_ref_unit_t))
                                                
                # combine to a closed flux space
                x_data = np.concatenate((ref_yield_ub, np.flip(ref_yield_lb), ref_yield_ub[:1]), axis=0)
                y_data = np.concatenate((fs_ub, np.flip(fs_lb), fs_ub[:1]), axis=0)
                    
            else:
                # calculate flux space from flux rates
                for grid_point in reference_grid:
                    # fix reference flux

                    ref_rxn.bounds = (grid_point, grid_point)
                    # ref_rxn.lower_bound = grid_point
                    # ref_rxn.upper_bound = grid_point
                    
                    # calculate lower and upper flux space bound
                    model_m.objective_direction = "max"
                    fs_ub = np.append(fs_ub, model_m.slim_optimize())
                    model_m.objective_direction = "min"
                    fs_lb = np.append(fs_lb, model_m.slim_optimize())
                    
                # combine to a closed flux space
                x_data = np.concatenate((reference_grid, np.flip(reference_grid), reference_grid[:1]), axis=0)
                y_data = np.concatenate((fs_ub, np.flip(fs_lb), fs_ub[:1]), axis=0)
                    
            # save and return data
            flux_space = {"x_data": x_data,
                            "y_data": y_data,
                            "x_label": reference_rxn,
                            "y_label": target_rxn
                }  
            
            if yield_space:
                flux_space["x_label"] = flux_space["x_label"] + " yield"
                flux_space["y_label"] = flux_space["y_label"] + " yield"
            
            self._flux_space = flux_space
            
            # plot projection
            if plot:
                # parse keyword arguments
                kwargs_parse = {}
                for key, arg in kwargs.items():
                    if key in ["dpi", "figsize"]:
                        kwargs_parse[key] = arg
                # plot data
                self.plot_2D(x_data,
                                y_data,
                                x_label=reference_rxn,
                                y_label=target_rxn,
                                save_name=save_name,
                                **kwargs_parse)
                # # fill accessible flux space with color
                # if fill_flux_space:
                #     self.axes.fill_between(x_data, 0, y_data, alpha=0.7)
                # # change labels of x and y axes
                # if xlabel:
                #     self.axes.set_xlabel(xlabel)
                # if ylabel:
                #     self.axes.set_ylabel(ylabel)
                
                # return self.axes

                # revert model changes
                # model_m.objective = objective_save
                # model_m.objective_direction = objective_direction_save
                # ref_rxn.bounds = ref_rxn_bounds_save
                
        return flux_space
        

    
    def plot_2D(self, x_data, y_data, x_label="", y_label="", save_name=None,
                dpi=300, figsize=(3.5, 3.5)):
        
        # change font to Arial
        try:
            plt.rcParams['font.sans-serif'] = "Arial"
        except:
            print('Font "Arial" not found for plotting')
        

        # 2D plot from flux data
        self.figure, self.axes = plt.subplots(1, 1,
                                              figsize=figsize,
                                              dpi=dpi)
        # self.plot = self.figure.add_subplot(111)
        self.axes.set_ylabel(y_label, fontsize=12)
        self.axes.set_xlabel(x_label, fontsize=12)
        # ticks
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
    
        # plot
        self.axes.plot(x_data, y_data,
                       linewidth=4)  
        
        self.axes.set_xticks(np.arange(0, np.max(x_data), np.round(np.max(x_data)/5, 2)))
        
        plt.tight_layout(pad=.8)

        # save plot as PNG
        if save_name:
            self.figure.savefig(save_name, dpi=400, edgecolor="white")

        
    def coupling_strength(self, target_rxn, reference_rxn=None, grid_size=15,
                          flux_space = {}):
        """
        calculate coupling strength between two fluxes based on a comparison
        between the accessible and inaccessible flux space
        """
        # calculate flux space
        if not flux_space:
            flux_space = self.flux_space_projection(
                target_rxn=target_rxn,
                reference_rxn=reference_rxn,
                grid_size=grid_size,
                plot=False)
            if not flux_space:
                # flux space could not be calculated
                return None
        
        
        def approximate_integral(x1, x2, y1, y2):
            return ((abs(x2-x1) * (abs(y2))) + (abs(x2-x1) * (abs(y1)))) / 2
        
        # get data
        x_data = flux_space["x_data"]
        y_data = flux_space["y_data"]
        
        # get maximum reference flux
        ref_flux_max = np.max(x_data)
        ref_flux_max_pos = np.where(x_data == ref_flux_max)
        # calculate integrals
        upper_boundary_integral = 0
        lower_boundary_integral = 0
        # approximate upper boundary integral
        for i in range(ref_flux_max_pos[0][0]):

            upper_boundary_integral += approximate_integral(x_data[i],
                                                            x_data[i+1],
                                                            y_data[i],
                                                            y_data[i+1])
        # approximate lower boundary integral
        for i in range(ref_flux_max_pos[0][1], len(x_data)-2):

            lower_boundary_integral += approximate_integral(x_data[i],
                                                            x_data[i+1],
                                                            y_data[i],
                                                            y_data[i+1])   
        if upper_boundary_integral <= 0 or not(upper_boundary_integral) \
            or not(lower_boundary_integral):
                # invalid integrals
                coupling_strength = None
        else:
            coupling_strength = lower_boundary_integral / upper_boundary_integral
            
        return coupling_strength
        
          
    def precursor_availability(self, target_rxn=None, disable_reactions=[]):
        # check if precursor of target reaction can be synthesized and return maximum (exchange) fluxes
        
        if not(target_rxn):
            # no target reaction provided, use objective
            for rxn in self.model.reactions:
                if rxn.objective_coefficient != 0:
                    target_rxn = rxn.id
                    print("Target reaction: " + target_rxn)
                    break                
        else:
            # check if target reaction exist in model
            if not(target_rxn in [rxn.id for rxn in self.model.reactions]):
                raise Exception("Target reaction " + target_rxn + " not in model")
        
        with self.model as model:
            # disable reactions
            model_rxns = [rxn.id for rxn in model.reactions]
            for rxn in disable_reactions:
                if rxn in model_rxns:
                    # disable reactions
                    model.reactions.get_by_id(rxn).bounds = (0, 0)
                else:
                    print("Reaction " + rxn + " not in model")
        
               
            # calculate target reaction precursor sink fluxes
            precursor_fluxes = {}
            target_rxn_rescue_flux = {}
            rxn = self.model.reactions.get_by_id(target_rxn)
            for met, coeff in rxn.metabolites.items():
                if coeff < 0:
                    # is a precursor or educt, create sink for metabolite
                    sink_rxn_id = "sink_" + met.id
                    sink_rxn = cobra.Reaction(
                        sink_rxn_id,
                        lower_bound=0,
                        upper_bound=1000)
                    sink_rxn.add_metabolites({met: -1})

                    with model:     
                        model.add_reactions([sink_rxn])
                        # calculate maximum sink flux
                        model.objective = sink_rxn_id 
                        precursor_fluxes[met.id] = model.slim_optimize()

                    # if precursor synthesis is blocked try to rescue target reaction 
                    if precursor_fluxes[met.id] == 0:   
                        with model: 
                            model.add_reactions([sink_rxn])
                            sink_rxn.lower_bound = -1000
                            model.objective = target_rxn
                            target_rxn_rescue_flux[met.id] = model.slim_optimize()
                        
        return precursor_fluxes, target_rxn_rescue_flux
    
    def precursor_bottleneck(self, target_rxn=None, disable_reactions=[]):
       # determine bottleneck in precursor supply by removing each precursor from (biomass) reaction 
       
        if not(target_rxn):
            # no target reaction provided, use objective
            for rxn in self.model.reactions:
                if rxn.objective_coefficient != 0:
                    target_rxn = rxn.id
                    print("Target reaction: " + target_rxn)
                    break                
        else:
            # check if target reaction exist in model
            if not(target_rxn in [rxn.id for rxn in self.model.reactions]):
                raise Exception("Target reaction " + target_rxn + " not in model")   
        
        
        with self.model as model:
            # disable reactions
            model_rxns = [rxn.id for rxn in model.reactions]
            for rxn in disable_reactions:
                if rxn in model_rxns:
                    # disable reactions
                    model.reactions.get_by_id(rxn).bounds = (0, 0)
                else:
                    print("Reaction " + rxn + " not in model")
                
                
            # one by one delete precursor from target reaction
            rxn = self.model.reactions.get_by_id(target_rxn)
            target_flux_after_removal = {}
            for met, coeff in rxn.metabolites.items(): 
                if coeff < 0:
                    with model:
                        # erase metabolite
                        model.reactions.get_by_id(target_rxn).add_metabolites({met: 0}, combine=0)
                        target_flux_after_removal[met.id] = model.slim_optimize()
                
        return target_flux_after_removal
                
                
    def flux_bound_limitations(self, fluxes=pd.Series(dtype="float64")):
        # detect limitations in the fluxes or flux distribution enforced by bounds
        
        # parameter
        flux_significant = 1e-4
        
        # check input
        if len(fluxes) == 0:
            # no fluxes provided, calculate fluxes by FBA
            sol = self.model.optimize()
            fluxes = sol.fluxes
         
        flux_limit_non_zero = pd.Series(dtype="float64")
        flux_limit_zero = pd.Series(dtype="float64")
        rxn_id_model = [rxn.id for rxn in self.model.reactions]
        with self.model as model:
            for rxn_id, flux in zip(fluxes.index, fluxes.values):
                # check if reaction exist in model
                if rxn_id in rxn_id_model:
                    # get reaction
                    rxn = model.reactions.get_by_id(rxn_id)
                    # get bounds
                    if abs(flux-rxn.lower_bound) < flux_significant \
                        or abs(flux-rxn.upper_bound) < flux_significant:
                            # flux is at the boundary limit
                            if abs(flux) < flux_significant:
                                flux_limit_zero[rxn_id] = flux
                            else:
                                flux_limit_non_zero[rxn_id] = flux
                                
        return flux_limit_non_zero, flux_limit_zero
    
    def flux_bound_limitation_changes(self, target_rxn,
                                      objective_direction="max",
                                      reference_rxn=None,
                                      grid_size=15):
        # detect changes in flux limitations for a range of fixed reference reaction flux rates
        
        # check inputs
        if not(target_rxn in [rxn.id for rxn in self.model.reactions]):
            raise Exception("Target reaction not in model")
            
        
        if not(reference_rxn):
            # no reference reaction provided, use objective
            for rxn in self.model.reactions:
                if rxn.objective_coefficient != 0:
                    reference_rxn = rxn.id
                    print("Reference reaction: " + reference_rxn)
                    break 
        
        with self.model as model:
            # get reference reaction range
            model.objective = reference_rxn
            model.objective_direction = "max"
            ref_max = model.slim_optimize()
            model.objective_direction = "min"
            ref_min = model.slim_optimize()
            ref_range = np.linspace(ref_min, ref_max, grid_size)
            
            # prepare model for target reaction optimization
            model.objective = target_rxn
            model.objective_direction = objective_direction
            # scan reference flux and compute flux limitations
            flux_limit_non_zero_prev = pd.Series(dtype="float64")
            flux_limit_zero_prev = pd.Series(dtype="float64")
            flux_limit_zero_changes = []
            ref_rxn = model.reactions.get_by_id(reference_rxn)
            for ref_flux in ref_range:
                ref_rxn.bounds = (ref_flux, ref_flux)
                sol = model.optimize()
                # compute flux limitations
                flux_limit_non_zero, flux_limit_zero \
                    = self.flux_bound_limitations(sol.fluxes)
                # determine changes in limitations
                changes_zero = pd.Series(dtype="float64")
                changes_non_zero = pd.Series(dtype="float64")
                if len(flux_limit_non_zero_prev) != 0 \
                    and len(flux_limit_zero_prev) != 0:
                        for idx in flux_limit_zero.index.difference(flux_limit_zero_prev.index):
                            changes_zero[idx] = flux_limit_zero[idx]
                        for idx in flux_limit_non_zero.index.difference(flux_limit_non_zero_prev.index):
                            changes_non_zero[idx] = flux_limit_non_zero[idx]
                 
                return_arg = {"reference": reference_rxn,
                              "reference_flux": ref_flux,
                              "target": target_rxn,
                              "target_flux": sol.objective_value,
                              "flux_limit_zero_changes": changes_zero,
                              "flux_limit_non_zero_changes": changes_non_zero}            
                 
                flux_limit_zero_prev = flux_limit_zero
                flux_limit_non_zero_prev = flux_limit_non_zero
                flux_limit_zero_changes.append(return_arg)
                
        return flux_limit_zero_changes
    
    def determine_auxotrophies(self, disable_reactions=[]):
        """
        Determine auxotrophies for an infeasible model or a model with an objective value of zero
        Method: Calculate shadow prices and try to rescue optimization problem 
                by adding metabolites with negative shadow prices

        Parameters
        ----------
        disable_reactions : TYPE, optional
            DESCRIPTION. The default is [].

        Returns
        -------
        auxotrophies : TYPE
            DESCRIPTION.
        sp_neg : TYPE
            DESCRIPTION.

        """
        #  
        # 
        
        with self.model as model:
            
            # disable reactions
            model_rxns = [rxn.id for rxn in model.reactions]
            for rxn in disable_reactions:
                if rxn in model_rxns:
                    # disable reactions
                    model.reactions.get_by_id(rxn).bounds = (0, 0)
                else:
                    print("Reaction " + rxn + " not in model")
                    
            # check objective value and retrieve shadow prices
            sol = model.optimize()
            if sol.status != "infeasible" and sol.objective_value != 0:
                print("Model is still feasible and objective value is non zero -> No auxotrophies")
                max_objective = sol.objective_value
            else:
                max_objective = 0
                
            # get negative shadow prices, indicates that metabolite source may resolve model infeasibility
            sp_neg = sol.shadow_prices[sol.shadow_prices < 0]
            # check if corresponding metabolite source enables model feasibility
            auxotrophies = {}
            for met in sp_neg.index:
                # protect model
                with model:
                    # create source reaction
                    source_rxn = cobra.Reaction("source_rxn",
                                                lower_bound=sp_neg[met],
                                                upper_bound=0)
                    source_rxn.add_metabolites({
                        model.metabolites.get_by_id(met): -1
                        })
                    # add to model
                    model.add_reactions([source_rxn])
                    # solve model
                    sol_aux = model.optimize()
                    # check solution
                    if not(sol_aux.objective_value):
                        # (still) infeasible solution
                        continue
                    if sol_aux.objective_value > max_objective:
                        # auxotrophy detected
                        auxotrophies[met] = {"ID": met,
                                             "sensitivity": 
                                                 (sol_aux.objective_value-max_objective)
                                                 /abs(sol_aux.fluxes["source_rxn"])}
                        
        return auxotrophies, sp_neg
                        
            
    def enforced_fluxes(self, objective=None):
        """
        determine reaction with enforced fluxes using flux variability analysis

        Parameters
        ----------
        objective : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        enforced_fluxes : DataFrame
            DESCRIPTION.

        """
        
        
        with self.model as model:
            # set objective
            if objective:
                model.objective = objective
                
            # conduct FVA
            sol_fva = cobra.flux_analysis.flux_variability_analysis(
                model,
                processes=1)
            
            # determine enforced fluxes
            enforced_fluxes = pd.DataFrame(columns=["minimum", "maximum", "minimum_enforced", "minimum_net"])
            for i, s in sol_fva.iterrows():
                if s["minimum"]*s["maximum"] > 1e-05:
                    if s["minimum"] > 0:
                        row = [s["minimum"], s["maximum"],s["minimum"], s["minimum"]]
                    elif s["maximum"] < 0:
                        row = [s["minimum"], s["maximum"],s["maximum"], -s["maximum"]]
                        
                    # add to dataframe
                    enforced_fluxes.loc[i, :] = row
                    
            # order by minimum net flux       
            enforced_fluxes = enforced_fluxes.sort_values(
                                by=["minimum_net"],
                                axis=0,
                                ascending=False
                                )
            
        return enforced_fluxes, sol_fva

            
                
               
     
                           
                
                    
        
            