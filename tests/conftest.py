import pytest

from pharmakin import formulas
from pharmakin import parameters
from pharmakin.utils.formula import Formulary


@pytest.fixture
def all_parameters():
    pars = parameters.get_all_parameters()
    return pars


@pytest.fixture
def all_formulas():
    res = formulas.get_all_formulas()
    return res


@pytest.fixture
def master_formulary(all_parameters, all_formulas):
    res = Formulary(
        parameters=all_parameters,
        formulas=all_formulas
    )
    
    return res


@pytest.fixture
def parameter_example_values(all_parameters):
    """Produces a list of dicts, mapping parameter names to examples of their values"""
    N = 100
    res = []
    for _ in range(N):
        d = {par.__name__:  par.example_values(with_units=False) for par in all_parameters}
        res.append(d)

    return res
