import json
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def matlab_lotka_volterra_reference():
    """Reference trajectory captured from the real MATLAB gsua_dpmat/gsua_deval, not a live call.

    See tests/data/matlab_lotka_volterra_reference.json for provenance -- regenerate by running
    that same MATLAB snippet again if the reference ever needs updating.
    """
    with open(DATA_DIR / "matlab_lotka_volterra_reference.json") as f:
        raw = json.load(f)
    return {
        "params": np.array(raw["params"]),
        "t": np.array(raw["t"]),
        "x": np.array(raw["x"]),
        "y": np.array(raw["y"]),
    }
