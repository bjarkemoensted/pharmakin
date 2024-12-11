from functools import wraps
import inspect
import sympy
from typing import Callable

from pharmakin.parameters import (
    Parameter,
    ParameterMeta,
    dose,
    auc,
    clearance,
    PARAMETER_REGISTRY
)
from pharmakin.units import has_units
from pharmakin.registry import Formulary


# TODO have parameter classes redirect lookups here, to get like clearence.compute(auc=5.0, dose=70.0) or something
formulary = Formulary()


def _convert_to_sympy_expression(func: Callable):
    """Takes a callable and converts it into a sympe expression for its result"""
    
    sig = inspect.signature(func)
    symbol_kws = {k: sympy.Symbol(k) for k in sig.parameters.keys()}
    expr = func(**symbol_kws)
    return expr



def formula(parameter: ParameterMeta):
    """Decorator for formulas expressing a parameter in terms of other parameters.
    Takes the class of the parameter computed by the formula as an input, i.e.
    @formula(my_parameter)
    def f(...):
        ...
        
    Adds an optional keyword with_units. If True/False this determines whether the result has units.
    If not set, result will have units only if all inputs have units.
    If not all or no inputs have unit, an error is raised.
    Also validates all input parameters"""
    
    # Only allow parameter classes as arguments to the decorator
    assert issubclass(parameter, Parameter)
    
    def outer(fun):
        
        # Register the function to the formulary
        rhs = _convert_to_sympy_expression(func=fun)
        lhs = sympy.Symbol(parameter.__name__)
        eq = sympy.Equality(lhs, rhs)
        formulary.register(eq=eq)
        
        @wraps(fun)
        def inner(*args, with_units: bool=None, **kwargs):
            # Convert the args and kwargs into a single dictionary, including any keyword defaults
            bound = inspect.signature(fun).bind(*args, **kwargs)
            bound.apply_defaults()
            kw_only = bound.arguments
            
            # Validate input parameters
            for k, v in kw_only.items():
                parcls = PARAMETER_REGISTRY[k]
                parcls.validate(v)
            
            # Ensure that either all or none of the inputs have declared types
            n_units = sum(has_units(val) for val in kw_only.values())
            consistent = n_units in (0, len(kw_only))
            if not consistent:
                raise RuntimeError(f"Units must be specified for all or no inputs (got {n_units}/{len(kw_only)})")
            
            # Compute the result
            res = fun(**kw_only)
            parameter.validate(res)
            
            # Enforce units/no units if with_unit keyword is used.
            if with_units is True:
                res = parameter.ensure_units(value=res)
            elif with_units is False:
                res = parameter.ensure_float(value=res)
            
            return res
        
        return inner
    return outer


@formula(clearance)
def clearence_from_dose_auc(dose, auc):
    """Determines clearence from dose and AUC"""
    res = dose / auc
    return res


if __name__ == '__main__':
    pass
