from functools import wraps
import inspect

from pharmakin import parameters
from pharmakin.units import has_units


def formula(parameter: parameters.ParameterMeta):
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
    assert issubclass(parameter, parameters.Parameter)
    
    def outer(fun):
        @wraps(fun)
        def inner(*args, with_units: bool=None, **kwargs):
            # Convert the args and kwargs into a single dictionary, including any keyword defaults
            bound = inspect.signature(fun).bind(*args, **kwargs)
            bound.apply_defaults()
            kw_only = bound.arguments
            
            # Validate input parameters
            for k, v in kw_only.items():
                parcls = parameters.PARAMETER_REGISTRY[k]
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


@formula(parameters.clearance)
def clearence_from_dose_auc(dose, auc):
    """Hi I am a docstring"""
    res = dose / auc
    return res



if __name__ == '__main__':
    pass
