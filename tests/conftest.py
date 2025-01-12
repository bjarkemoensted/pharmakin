import pytest

from pharmakin import parameters


@pytest.fixture
def all_parameters():
    pars = parameters.get_all_parameters()
    return pars
