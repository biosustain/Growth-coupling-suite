"""
default config file for heterologous reaction extraction
"""

# %% heterologous reaction extraction

# assess directions of heterologous reaction by thermodynamic and flux variability analysis
directionality_assessment = True

# if thermodynamic assessment of reactions is done, only consider reactions for which
# Gibbs free energy of reaction is accessable
consider_thermodynamic_valid_reactions_only = True

num_carbons_threshold = 35

currency_mets = {'h2o', 'co2', 'o2', 'h2o2', 'nh4', 'no2', 'no3', 'no', 'h2s',
                 'so3','so4','h','h2','pi','ppi','coa','accoa','ppcoa','aacoa',
                 'butcoa','succoa','atp','gtp','adp','gdp','amp','gmp','nad',
                 'nadp','nadh','nadph','fad','fadh','na1','ahcys','amet','thf',
                 'mlthf', 'q8h2','q8','mql8','mqn8','2dmmql8','2dmmq8'}

# only consider heterologous reactions with exclusive compartments in host?
only_host_compartments = True

# %% 
