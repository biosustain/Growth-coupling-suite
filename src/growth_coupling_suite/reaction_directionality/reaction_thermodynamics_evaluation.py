"""
Evaluate thermodynamics of reactions given in a COBRA model format
"""
import numpy as np
import pandas as pd
# import cobra
from equilibrator_api import ComponentContribution, Q_

from growth_coupling_suite.reaction_directionality \
    import config_reaction_directionality_default as config_default

from growth_coupling_suite.external_APIs import keggAPI

from growth_coupling_suite.util.util import check_config

# %%
cc = ComponentContribution()  # initiate group contribution method
# set standard parameters
cc.ionic_strength = Q_(config_default.ionic_strength)
cc.temperature = Q_(config_default.temperature)
cc.p_h = Q_(config_default.p_h)

# %% initialize kegg
kegg = keggAPI.KeggAPI()


# %%
def get_reaction_direction(rxn, config=config_default, eq_rxn_in=[], r_measure_in=[]):
    """
    
    Parameters
    ----------
    rxn : core.reaction.Reaction
        COBRA model reaction
    config : module, optional
        configuration file/module including options. The default is config_default.
    eq_rxn_in : phased_reaction.PhasedReaction, optional
        phased reaction format from equilibrator API. The default is [].
    r_measure_in : float or list, optional
        default reversibility measure (dG or RI). The default is [].
    
    Returns
    -------
    direction: int
        predicted allowable reaction directionality
            1: irreversible in forward direction
            -1: irreversible in backward direction
            0: reversible
    eq_rxn: phased_reaction.PhasedReaction
        phased reaction format from equilibrator API
    r_measure: float or list
        reversibility measure on which direction was evaluated
    
    """
    
    

    # if check_config:
    check_config(config, config_default)

    # set up equilibrator reaction format
    if not(eq_rxn_in):
        eq_rxn = []
    else:
        eq_rxn = eq_rxn_in.clone()
        
    # set provided parameter: ionic strength, pH, temperature
    set_eq_parameters(config)
    
    # extract parameter
    
    
    if not(len(r_measure_in)>0):
        # no reversibility measure values provided
        
        
        
        if config.reaction_direction_eval_type == "RI":
            """ 
            calculate and evaluiate reaction reversibility measure
            """
            # use reversibility index as reversibility measure
            # calculate reversibility index 
            if config.use_experimental_concentrations:
                # load metabolite concentrations
                met_concentrations \
                    = load_metabolite_concentrations(config.metabolite_concentration_filename)
                # calculate reversibility index and standard deviation
                ln_ri, ln_ri_std, ln_ri_eq, eq_rxn \
                    = get_physiological_reversibility_index(rxn,
                                        met_concentrations=met_concentrations,
                                        config=config, eq_rxn_in=eq_rxn)
            else:
                # calculate reversibility index and standard deviation
                ln_ri, ln_ri_std, ln_ri_eq, eq_rxn \
                    = get_reversibility_index(rxn, config=config, eq_rxn_in=eq_rxn)
        
            # evaluate standard deviation and reliability of dG
            if ln_ri_std>=abs(ln_ri):
                # standard deviation greater than actual value
                raise TypeError("reversibility index is not reliable: " 
                                + str(ln_ri) + " +- " + str(ln_ri_std))
        
        
            # consider standard deviation for approximating worst case szenario
            if config.consider_standard_deviation:
                ln_ri += ln_ri_std-(2*ln_ri_std*(ln_ri>0))

            # check "worst" case reversibility indices
            if ln_ri<-config.ri_threshold:
                # reaction is irreversible in forward direction
                direction = 1
            elif ln_ri>config.ri_threshold:
                # reaction is irreversible in backward direction
                direction = -1
            else:
                # reaction is reversible
                direction = 0
                
            r_measure = ln_ri
        
        elif config.reaction_direction_eval_type == "dG":
            """ 
            calculate and evaluate Gibbs free energy of reaction
            """
            # use Gibbs free energy of reaction as reversibility measure
            # calculate Gibbs free energy of reaction 
            if config.use_experimental_concentrations:
                # load metabolite concentrations
                met_concentrations \
                    = load_metabolite_concentrations(config.metabolite_concentration_filename)
                    
            else:
                # no physiological metabolite concentrations provided
                met_concentrations = pd.DataFrame() # empty data frame
                        
                    
            if config.evaluate_metabolite_range:
                # evaluate min/max dG prime according to expected metabolite ranges
                dG_range, dG_range_std, dG_range_eq, eq_rxn \
                   = get_dG_prime_range(rxn, met_conc=met_concentrations,
                                   config=config, eq_rxn_in=eq_rxn) 
                 
                   
                # evaluate standard deviation and reliability of dG
                if dG_range_std[0]>=abs(dG_range[0]) and dG_range_std[1]>=abs(dG_range[1]):
                    # standard deviation greater than actual value
                    print("Gibbs free energy of " + rxn.id +  " is not reliable" )
                    raise TypeError("Gibbs free energy of reaction is not reliable" 
                                    + str(dG_range) + " +- " + str(dG_range_std))
                
                # construct a reversibility measure
                dG_range = np.array(dG_range)
                dG_range_std = np.array(dG_range_std)
                
                # consider standard deviation?
                if config.consider_standard_deviation:
                    dG_range_std = dG_range_std[np.argsort(dG_range)]
                    dG_range = np.sort(dG_range)
                    dG_range[0] -= dG_range_std[0]
                    dG_range[1] += dG_range_std[1]
                    
                
                if np.all(dG_range<0):
                    # reaction is irreversible in forward direction
                    direction = 1
                elif np.all(dG_range>0):
                    # reaction is irreversible in backward direction
                    direction = -1
                else:
                    # reaction is reversible
                    direction = 0
                    
                    # evaluate general dG threshold
                    dG_threshold_general = config.dG_threshold
                    if sum(dG_range)<-dG_threshold_general:
                        direction = 1
                    elif sum(dG_range)>dG_threshold_general:
                        direction = -1
             
                    # consider heuristics
                    # energy cofactors, position denotes hierarchy of energy equivalent      
                    ntp_ids = ["atp", 
                               "utp", 
                               "ctp", 
                               "xtp", 
                               "itp", 
                               "gtp",
                               "ditp",
                               "dgtp"]
                    coeff_ntp = np.zeros(len(ntp_ids))

                    coeff_quinone = 0 # quinone coefficient
                    coeff_co2 = 0 # co2 coefficient
                    coeff_hco3 = 0 #hco2 coefficient
                    for met in rxn.metabolites:
                        met_id = met.id[:-2]
                        met_name = met.name
                        # is nucleosid triphosphate in reaction formula?
                        for i in range(len(ntp_ids)):
                            if met_id == ntp_ids[i]:
                                coeff_ntp[i] = rxn.metabolites.get(met)

                        # quinone in reaction formula?
                        if "quinone" in met_name:
                            coeff_quinone = rxn.metabolites.get(met) 
                            
                        if "co2" in met_id:
                            coeff_co2 = rxn.metabolites.get(met)
                            
                        if "hco3" in met_id:
                            coeff_hco3 = rxn.metabolites.get(met) 
                    
            
                    dG_cutoff_energy = 0 
                    if np.any(coeff_ntp<0) and np.any(coeff_ntp>0):
                        # nucleotid triphosphates are produced and consumed
                        direction = 0
                    else:
                        for coeff in coeff_ntp:
                            if coeff!=0:
                                direction = evaluate_cofactor_heuristics(coeff,
                                                                         dG_range,
                                                                         dG_cutoff_energy)
                                break
                            
                    # co2 heuristics
                    dG_cutoff_co2 = 0
                    if coeff_co2 != 0:
                        if coeff_hco3 != 0:
                            direction = 0
                        else:
                            if coeff_co2 < 0:
                                direction = -1
                                if sum(dG_range) < -dG_cutoff_co2:
                                    direction = 0
                                    
                            else:
                                direction = 1
                                if sum(dG_range) > dG_cutoff_co2:
                                    direction = 0
                            
                        
                    
                    
                    # quinones: reactions are irreversible in direction of the quinol
                    # if thermodynamic strongly points in the other direction, reverse
                    dG_cutoff_quinone = 30
                    if coeff_quinone < 0:
                        if sum(dG_range)<dG_cutoff_quinone:
                            direction = 1
                        else:
                            direction = 0
                    elif coeff_quinone > 0:
                        if sum(dG_range)<-dG_cutoff_quinone:
                            direction = 0
                        else:
                            direction = -1
                           
                r_measure = list(dG_range)
                        

                
                
            else:
                dG_prime, dG_prime_std, dG_prime_eq, eq_rxn \
                    = get_physiological_dg_prime(rxn,
                                                met_concentrations=met_concentrations,
                                                config=config, eq_rxn_in=eq_rxn)
                    
                
                # consider standard deviation?
                if config.consider_standard_deviation:
                    dG_prime += dG_prime_std-(2*dG_prime_std*(dG_prime>0))
                    
                # evaluate dG
                if abs(dG_prime)<config.dG_threshold:
                    direction = 0
                else:
                    if dG_prime<0:
                        direction = 1
                    else:
                        direction = -1
                        
                r_measure = dG_prime
                
        else:
            raise TypeError("Reaction evaluation type is unknown. Choose RI or dG")
              
    
            
    return direction, eq_rxn, r_measure

    
def create_equilibrator_reaction(rxn):
    # define parameter
    bigg_parser = "bigg.metabolite:{0}"
    kegg_parser = "kegg:{0}"
    # create default equilibrator reaction
    eq = cc.parse_reaction_formula("=")
    # make metabolite ids conform to equilibrator
    for met in rxn.metabolites:
        # strip compartment sign
        met_id = met.id[:-2]
        # check if id is known to equilibrator
        eq_met = cc.get_compound(bigg_parser.format(met_id))
        if not(eq_met):
            print("\tMetabolite " + met.id + " not known in BIGG namespace")
            # if it is not known query kegg
            cmp_kegg = kegg.get("compound", met.name, prioritize=True, best_guess=True)
            if not(cmp_kegg):
                # query kegg unsuccessful
                print("\tUnkown compound encountered: " + met.id + " (" + met.name + ")")
                raise TypeError("Compound not found")
            else:
                if "IDs" in list(cmp_kegg.keys()):
                    if not(cmp_kegg["IDs"]):
                        # query kegg unsuccessful
                        print("\tUnkown compound encountered: " + met.id + " (" + met.name + ")")
                        raise TypeError("Compound not found")
                    else:
                        kegg_id = cmp_kegg["IDs"]
                        
                elif "ID" in list(cmp_kegg.keys()):
                    if not(cmp_kegg["ID"]):
                    # query kegg unsuccessful
                        print("\tUnkown compound encountered: " + met.id + " (" + met.name + ")")
                        raise TypeError("Compound not found")
                    else:
                        kegg_id = cmp_kegg["ID"]
                    
                        
                # get metabolite
                print("\t" + met.id + " found in KEGG namespace")
                eq_met = cc.get_compound(kegg_parser.format(kegg_id))        
            
        # add stoichiometry to equilibrator reaction
        eq.add_stoichiometry(eq_met, rxn.get_coefficient(met.id))



    return eq


def get_physiological_reversibility_index(rxn, met_concentrations = pd.DataFrame(),
                                             config=config_default, eq_rxn_in=[]):
    
    
    # check config module
    check_config(config, config_default)
    
    # set up equilibrator reaction format
    if not(eq_rxn_in):
        eq_rxn = create_equilibrator_reaction(rxn)
        # # check reaction
        # if eq_rxn == -1:
        #     raise TypeError("Creation of equilibraotr reaction format failed!")
    else:
        eq_rxn = eq_rxn_in.clone()
    
    
    # calculate physiological dG prime considering metabolite concentrations
    phys_dG_prime_val, phys_dG_prime_std, physiological_dg_prime, eq_rxn_dump \
        = get_physiological_dg_prime(rxn, met_concentrations=met_concentrations,
                                     config=config,
                                     eq_rxn_in=eq_rxn)
    # get absolute sum of coefficients
    abs_sum_coeff = eq_rxn._sum_absolute_coefficients()
    
    # print(physiological_dg_prime)
    # calculate ln reversibility index
    ln_RI_eq = (2.0 / abs_sum_coeff) * physiological_dg_prime / cc.RT
    
    # process reversibility index
    ln_ri_val, ln_ri_std, ln_ri_min_val = process_ln_ri(ln_RI_eq)

    
    return ln_ri_val, ln_ri_std, ln_RI_eq, eq_rxn


    
def get_physiological_dg_prime(rxn, met_concentrations = pd.DataFrame(),
                               config=config_default, eq_rxn_in=[]):
    
    # (manually) determine physiological logarithmic reversibility index
    
    # check config module
    check_config(config, config_default)
    
    # set up equilibrator reaction format
    if not(eq_rxn_in):
        eq_rxn = create_equilibrator_reaction(rxn)
        # # check reaction
        # if eq_rxn == -1:
        #     raise TypeError("Creation of equilibraotr reaction format failed!")
    else:
        eq_rxn = eq_rxn_in.clone()

    # load standard E. coli metabolite concentrations if none provided
    if len(met_concentrations) == 0:
        met_names = []
    else:
        # get names of metabolites in concetrations data frame
        met_names = list(met_concentrations.Met)

    
    # set provided physiological metabolite concentrations and calculate absolute sum of coefficients
    for cmp, coeff in eq_rxn.items(protons=False, water=False):
                          
        # compare compound ID with provided concentration list and extract concentration
        # use standard concentration if physiological concentration is not provided
        concentration = config.standard_physiological_concentration
        identifier_list = cmp.identifiers
        for namespace in identifier_list:
            if namespace.accession in met_names:
                selection = met_concentrations.loc[met_concentrations.Met == namespace.accession]
                concentration = selection.iloc[0, 1]
                break
         
  
        # set concentration
        # print(namespace.accession + " " + str(concentration))
        eq_rxn.set_abundance(cmp, Q_(str(concentration) + " mM"))
        
        
    
    # calculate physiologica dG
    physiological_dg_prime  = cc.dg_prime(eq_rxn)
        
    # process dG prime
    phys_dG_prime_val, phys_dG_prime_std = process_dG(physiological_dg_prime)
        
    return phys_dG_prime_val, phys_dG_prime_std, physiological_dg_prime, eq_rxn
     

def get_dG_prime_range(rxn, met_conc=pd.DataFrame(), config=config_default, eq_rxn_in=[]):   
    
    # manually calculate dG prime range assuming a range of metabolite concentrations
     
    # check config module
    check_config(config, config_default)
    
    # set up equilibrator reaction format
    if not(eq_rxn_in):
        eq_rxn = create_equilibrator_reaction(rxn)
        # # check reaction
        # if eq_rxn == -1:
        #     raise TypeError("Creation of equilibraotr reaction format failed!")
    else:
        eq_rxn = eq_rxn_in.clone()

    # load metabolite concentration if provided
    # get names of metabolites in concetrations data frame
    if len(met_conc) > 0:
        met_names = list(met_conc.Met)
    else:
        met_names = []


    # parameter
    met_conc_max = max(config.metabolite_concentration_range)
    met_conc_min = min(config.metabolite_concentration_range)
    # concentration ranges for dissolved gasses [mM] (Feist et al. 2008)
    gasses_conc_max = {"co2": 1.4,
                       "o2": 0.055,
                       "h2": 0.034} 
    gasses_conc_min = {"co2": 0.00001,
                       "o2": 0.00001,
                       "h2": 0.00001} 
    gasses_id = list(gasses_conc_max.keys())
    
    # set concentration for maximum driving force towards products
    cmp_exclude_id = []
    for cmp, coeff in eq_rxn.items(protons=False, water=False):
        # check and use if metabolite concentration is provided
        identifier_list = cmp.identifiers
        concentration = -1
        for namespace in identifier_list:
            if namespace.accession in met_names:
                selection = met_conc.loc[met_conc.Met == namespace.accession]
                concentration = selection.iloc[0, 1]
                break
        if concentration >= 0:
            eq_rxn.set_abundance(cmp, Q_(str(concentration) + " mM"))
            # save compound ID
            cmp_exclude_id.append(cmp.id)
            continue
           
        # set concentrations for dissolved gasses (Feist et al. 2008)
        is_gas = False
        for namespace in identifier_list:
            if namespace.accession in gasses_id:
                is_gas = True
                if coeff < 0:
                    eq_rxn.set_abundance(cmp,
                        Q_(str(gasses_conc_max[namespace.accession]) + " mM"))
                elif coeff > 0:
                    eq_rxn.set_abundance(cmp,
                        Q_(str(gasses_conc_min[namespace.accession]) + " mM"))
        if is_gas:
            continue
        
        # set concentration for maximum driving force towards products and educts
        if coeff < 0:
            # set educt concentration
            eq_rxn.set_abundance(cmp, Q_(str(met_conc_max) + " mM"))
        elif coeff > 0:
            # set product concentration
            eq_rxn.set_abundance(cmp, Q_(str(met_conc_min) + " mM"))
            
    # calculate dG prime
    try:
        dG_prime_forward = cc.dg_prime(eq_rxn)
    except:
        print("\tAn error occured during calculation of dG_prime")
        raise TypeError("dG_prime calculation failed")
        
     
    # set concentration for maximum driving force towards educts
    for cmp, coeff in eq_rxn.items(protons=False, water=False):
        # exclude compound if concentration was set before
        if cmp.id in cmp_exclude_id:
            continue
        
        # set concentrations for dissolved gasses (Feist et al. 2008)
        is_gas = False
        for namespace in identifier_list:
            if namespace.accession in gasses_id:
                is_gas = True
                if coeff < 0:
                    eq_rxn.set_abundance(cmp,
                        Q_(str(gasses_conc_min[namespace.accession]) + " mM"))
                elif coeff > 0:
                    eq_rxn.set_abundance(cmp,
                        Q_(str(gasses_conc_max[namespace.accession]) + " mM"))
        if is_gas:
            continue
        
        # set concentration for maximum driving force towards products and educts
        if coeff < 0:
            # set educt concentration
            eq_rxn.set_abundance(cmp, Q_(str(met_conc_min) + " mM"))
        elif coeff > 0:
            # set product concentration
            eq_rxn.set_abundance(cmp, Q_(str(met_conc_max) + " mM"))   
            
    # calculate dG prime
    dG_prime_backward = cc.dg_prime(eq_rxn)
     
    # transform results
    dG_prime_forward_val, dG_prime_forward_std = process_dG(dG_prime_forward)
    dG_prime_backward_val, dG_prime_backward_std = process_dG(dG_prime_backward)
     
    return [dG_prime_forward_val, dG_prime_backward_val],\
        [dG_prime_forward_std, dG_prime_backward_std],\
            [dG_prime_forward, dG_prime_backward], eq_rxn
 
    

def get_reversibility_index(rxn, config=config_default, eq_rxn_in=[]):
    # if check_config:
    check_config(config, config_default)


    # set up equilibrator reaction format
    if not(eq_rxn_in):
        eq_rxn = create_equilibrator_reaction(rxn)
        # # check reaction
        # if eq_rxn == -1:
        #     raise TypeError("Creation of equilibraotr reaction format failed!")
    else:
        eq_rxn = eq_rxn_in.clone()



    # determine reversibility index
    if config.check_pH_range:
        # test for a range of pH values
        num_pH = 10
        pH_range = np.linspace(config.pH_limits[0], config.pH_limits[1], num_pH)
        ln_ri = np.empty(num_pH)
        for i in range(num_pH):
            # change pH value
            cc.p_h = Q_(str(pH_range[i]))
            # calculate reversibility index
            ln_ri_out = cc.ln_reversibility_index(eq_rxn)
            ln_ri_val, ln_ri_std, ln_ri_min_val = process_ln_ri(ln_ri_out)
            # print([ln_ri_val, ln_ri_std, ln_ri_min_val])
            # consider standard deviation?
            if config.consider_standard_deviation:
                ln_ri[i] = ln_ri_min_val
            else:
                ln_ri[i] = ln_ri_val
            
    else:
        ln_ri_eq = cc.ln_reversibility_index(eq_rxn)
        ln_ri_val, ln_ri_std, ln_ri_min_val = process_ln_ri(ln_ri_eq)            
            
    return ln_ri_val, ln_ri_std, ln_ri_eq, eq_rxn


# %% Additional functions
def set_eq_parameters(config):
    # set important equilibrator parameters
    cc.ionic_strength = Q_(config.ionic_strength)
    cc.temperature = Q_(config.temperature)
    cc.p_h = Q_(config.p_h)

def load_metabolite_concentrations(filename):
    
    if ".xlsx" in filename:
        # load Excel file as data frame
        df = pd.read_excel(filename, engine='openpyxl')
    elif ".csv" in filename:
        # load .csv file as data frame
        df = pd.read_csv(filename)
    else:
        raise TypeError("Unknown file format. Use .xlsx or .csv format")

    # strip compartment sign from metabolite ID
    for i in range(len(df.index)):
        if df.iloc[i,0][-2] == "_":
            # "_c" format
            df.iloc[i,0] = df.iloc[i,0][:-2]
        elif df.iloc[i,0][-3] == "[" and df.iloc[i,0][-1] == "]":
            # "[c]" format
            df.iloc[i,0] = df.iloc[i,0][:-3]
    
    return df

def process_ln_ri(ln_ri):
    ln_ri = ln_ri.m_as("")
    ln_ri_val = ln_ri.nominal_value
    ln_ri_std = ln_ri.std_dev
    # calculate smallest necessary concentration fold change to change reaction direction
    # within an 95% probability range
    if ln_ri_val>0:
        ln_ri_min = ln_ri_val - ln_ri_std
    else:
        ln_ri_min = ln_ri_val + ln_ri_std
        
    return ln_ri_val, ln_ri_std, ln_ri_min

def process_dG(dG):
    dG = dG.m_as("kJ/mol")
    dG_val = dG.nominal_value
    dG_std = dG.std_dev

    return dG_val, dG_std
            
def evaluate_cofactor_heuristics(coeff, dG_range, dG_cutoff):
    # evaluate heursitics fpr cofactors such as ATP, GTP, quinone etc
    if coeff < 0:
        direction = 1
        if sum(dG_range)>-dG_cutoff:
            direction = 0
    else:
        direction = -1
        if sum(dG_range)<dG_cutoff:
            direction = 0
      
    return direction
