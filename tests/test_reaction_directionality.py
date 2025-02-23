import cobra
import pytest
from pathlib import Path
from growth_coupling_suite.reaction_directionality.reaction_thermodynamics_evaluation import get_reaction_direction

# load model
model_path = Path(__file__).parent / "Files" / "ECC2.json"
model = cobra.io.load_json_model(model_path)

# pytest parameters
@pytest.mark.parametrize("model, rxn_id, expected", [
    (model, "G6PDH2r", 0), 
    (model, "PFK", 1)
])

def test_reaction_direction(model, rxn_id, expected):
    direction_thermodynamic, eq_rxn, reversibility_measure \
        = get_reaction_direction(model.reactions.get_by_id(rxn_id))       
    assert direction_thermodynamic == expected
        
    
