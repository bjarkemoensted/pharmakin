import pytest

from pharmakin import parameters


@pytest.fixture
def all_parameters():
    pars = [parameters.dose, parameters.auc, parameters.clearance]  # !!!
    return pars