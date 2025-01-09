import pytest

from pharmakin.utils import parameters


@pytest.fixture
def all_parameters():
    pars = [parameters.dose, parameters.auc, parameters.clearance]  # !!!
    return pars