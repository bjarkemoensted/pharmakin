import pytest

from pharmakin import formulas
from pharmakin import parameters
from pharmakin.utils.registry import Formulary


@pytest.fixture
def all_parameters():
    pars = parameters.get_all_parameters()
    return pars


@pytest.fixture
def all_formulas():
    res = formulas.get_all_formulas()
    return res


# TODO MAKE SOME TESTS VERIFYING THAT WE HAVE DEFINITIONS FOR ALL PARAMS ETC!!!
def one_formulary_to_rule_them_all(all_parameters, all_formulas):
    Formulary(
        parameters=all_parameters,
        formulas=all_formulas
    )
