# Growth-Coupling Suite
The Growth-Coupling Suite is a framework for computing and analyzing strain designs that couple a target reaction to growth. [gcOpt](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2946-7) is used as the underlying optimization algorithm for deriving growth-coupled strain design solutions.  

## Installation
For a detailed installation manual, refer to the "docs" folder in this repository. The Growth Coupling Suite has been tested in Python >3.8. It is recommended to create a virtual python environment with, e.g., conda or virtualenv, for installing and using the Growth Coupling Suite.
1. Clone/fork/download this repository
2. Browse to the main directory of the repository and run `pip install .` or `pip install -e .` (for code development purposes). This will install all packages (cf. requirements.txt)
3. (optional) If a heterologous reactiond database model needs to be built, run `python -c "from equilibrator_api import ComponentContribution; cc = ComponentContribution()"` to initialize the equilibrator API ([documentation](https://equilibrator.readthedocs.io/en/latest/index.html))

## Use
Refer to and run the example scripts in 'examples' (conduct_gcOpt_optimization_parallel, conduct_gcOpt_optimization_single) for an introduction setting up and using the Growth Coupling Suite.  
Once design solutions were found and analyzed with the StrainDesignAnalyzer (part of the Growth Coupling Suite), a summary including unique, valid strain designs and their metadata is saved in an Excel format to the specified location (cf. the parameter 'results_dir' of the 'growth_coupling_summary' function). 

## Note  
The first computation with a new model can be quite time consuming due to the curation of the heterologous reaction database, if heterologous insertions are allowed (num_addins>0 in the gcOpt_config_file). The database will automatically be saved for later applications of that model.
A heterologous reaction database model can be manually built and saved beforehand by using the 'heterologous_reaction_processing' module. Refer to 'tests/create_heterologous_reaction_database_model.ipynb' for a respective example.

## Citation

If you use the Growth-Coupling Suite in your research or find it helpful for your work, please cite the following paper:

> **Metabolic growth-coupling strategies for in vivo enzyme selection systems**  
> Tobias B. Alter, Pascal A. Pieters, Colton J. Lloyd, Adam M. Feist, Emre Özdemir, Bernhard O. Palsson, Daniel C. Zielinski  
> *Metabolic Engineering Communications, 2025*  
> DOI or URL: [https://doi.org/10.1016/j.mec.2025.e00257](https://doi.org/10.1016/j.mec.2025.e00257)
