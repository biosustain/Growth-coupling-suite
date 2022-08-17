# Growth Coupling Suite
Framework for computing and analyzing strain designs that couple a target reaction to growth. ([gcOpt](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2946-7)) is used as the underlying optimization algorithm for deriving growth-coupled strain design solutions.  

## Installation
For a detailed installation manual, refer to the "docs" folder in this repository. The Growth Coupling Suite has been tested in Python >3.8. It is recommended to create a virtual python environment with, e.g., conda or virtualenv, for installing and using the Growth Coupling Suite.
1. Install COBRApy (run `pip install cobra`) ([documentation](https://cobrapy.readthedocs.io/en/latest/))
2. Install the Gurobi solver software under an (academic) license (https://www.gurobi.com/)
3. Install gurobipy (run `pip install gurobipy`) ([documentation](https://www.gurobi.com/documentation/9.1/quickstart_mac/cs_using_pip_to_install_gr.html))
4. If a heterologous reactiond database model needs to be built, install the Equilibrator API (run `pip install equilibrator-api`) ([documentation](https://equilibrator.readthedocs.io/en/latest/index.html))
5. Clone/fork/download this repository
6. Browse to the main directory of the repository and run `python setup.py install` or `python setup.py develop` (for code development purposes)

## Use
Refer to and run the example scripts in 'examples' (conduct_gcOpt_optimization_parallel, conduct_gcOpt_optimization_single) for an introduction setting up and using the Growth Coupling Suite.  
Once design solutions were found and analyzed with the StrainDesignAnalyzer (part of the Growth Coupling Suite), a summary including unique, valid strain designs and their metadata is saved in an Excel format to the specified location (cf. the parameter 'results_dir' of the 'growth_coupling_summary' function). 

## Note  
The first computation with a new model can be quite time consuming due to the curation of the heterologous reaction database, if heterologous insertions are allowed (num_addins>0 in the gcOpt_config_file). The database will automatically be saved for later applications of that model.  
A heterologous reaction database model can be manually built and saved beforehand by using the 'heterologous_reaction_processing' module. Refer to 'tests/create_heterologous_reaction_database_model.ipynb' for a respective example.
