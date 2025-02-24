import pytest
import cobra
from pathlib import Path
from growth_coupling_suite.gcOpt_algorithm import gcOpt

# load model
model_path = Path(__file__).parent / "Files" / "ECC2.json"
model = cobra.io.load_json_model(model_path)


def test_gcopt_initialization():
    gcopt_instance = gcOpt.GCOpt(
        model,
        model.reactions.get_by_id("EX_ac_e"),
        build_gcopt_problem = False        
        )
    assert gcopt_instance.model == model
    
    
    
if __name__ == "__main__":
    test_gcopt_initialization()