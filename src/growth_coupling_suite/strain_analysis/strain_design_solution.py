# simple class for storing design solutions
from types import SimpleNamespace

class Solution():
    def __init__(self, DesignAnalyzer, design_solution, config=SimpleNamespace(),
                 filename="", solution_key="", key_in_file=""):
        """
        Object for strain design solutions

        Parameters
        ----------
        DesignAnalyzer : 
            Instance of parent StrainDesignAnalyzer class.
        design_solution : dict
            Specifications of design solutions.
        config : SimpleNamespace, optional
            config parameter of strain design optimization algorithm
        filename : str, optional
            filename where design is saved. The default is "".
        solution_key : str, optional
            key of design solution in  parent namespace. The default is "".

        Returns
        -------
        None.

        """
        
        # save key of solution from DesignAnalyzer class
        self._solution_key = solution_key
        # save of of solution from file
        self._solution_key_in_file = key_in_file
        # get instance of StrainDesignAnalyzer
        self._DesignAnalyzer = DesignAnalyzer
        # save config
        self._config = config
        # save and check design solution
        self._design_solution = design_solution
        self._check_solution_consistency()
        # save filename
        self.filename = filename
        
        
        
        
    def load(self):
        # load interventions of solution into StrainDesignAnalyzer
        # check if medium defiition exists
        if "medium" in list(self._design_solution.keys()):
            medium = self._design_solution["medium"]
        else:
            medium = {}
        
        # load design in DesignAnalyzer
        self._DesignAnalyzer.load_strain_design(
            design=self._design_solution["interventions"],
            medium=medium,
            name=self._solution_key)
        
    def _check_solution_consistency(self):
        # check consistency of design solution, adapt if necessary
        
        # set of interventions have to be provided
        is_interventions = False
        if "interventions" in self._design_solution:
            if len(self._design_solution["interventions"]) > 0:
                is_interventions = True
        if not(is_interventions):
            print("\tNo interventions in strain design " + self._solution_key)
            self._design_solution["interventions"] = {}
            
        # check for provided medium
        is_medium = False
        if "medium" in self._design_solution:
            if len(self._design_solution["medium"]) > 0:
                is_medium = True
        if not(is_medium):
            print("\tNo medium specifications in strain design " + self._solution_key)
            self._design_solution["medium"] = {}
        