import pytest

from pharmakin.parameters import PARAMETER_REGISTRY


@pytest.fixture
def all_parameters():
    pars = [par for par in PARAMETER_REGISTRY.values()]
    return pars