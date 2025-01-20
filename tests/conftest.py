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
def _parameter_examples(all_parameters):
    """Produces a list of dicts, mapping parameters to examples of their values (with units)"""
    N = 100
    res = []
    for _ in range(N):
        d = {par:  par.example_values(with_units=True) for par in all_parameters}
        res.append(d)

    return res


@pytest.fixture
def parameter_names_with_example_vals(_parameter_examples):
    """Produces a list of dicts, mapping parameter names to examples of their values"""
    
    return [{par.__name__: par.ensure_float(val) for par, val in d.items()} for d in _parameter_examples]


@pytest.fixture
def example_calculations_from_formulas(master_formulary):
    n = 5
    with_units = True
    res = []
    
    for formula in master_formulary.formulae:
        for _ in range(n):
            input_parameters = [master_formulary.parameters[p] for p in formula.inputs]
            d = {par: par.example_values(with_units=with_units) for par in input_parameters}
            
            kwds = {p.__name__: v for p, v in d.items()}
            computed_par = formula(**kwds)
            d[formula.result_class] = computed_par
            
            res.append(d)
    
    return res