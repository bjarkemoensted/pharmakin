from __future__ import annotations
from functools import wraps
import inspect
import sympy
from typing import Callable, Iterable, TYPE_CHECKING

if TYPE_CHECKING:
    from pharmakin.utils.parameter_base import ParameterMeta, Parameter
from pharmakin.utils.units import has_units


def _convert_to_sympy_expression(func: Callable):
    """Takes a callable and converts it into a sympe expression for its result"""
    
    sig = inspect.signature(func)
    symbol_kws = {k: sympy.Symbol(k) for k in sig.parameters.keys()}
    expr = func(**symbol_kws)
    return expr


def format_sympy_equation(eq: sympy.core.relational.Equality) -> str:
    """Formats a sympy equation into a string."""

    if not isinstance(eq, sympy.core.relational.Equality):
        raise TypeError
    
    res = f"{eq.lhs} = {eq.rhs}"
    return res


def make_symbol_key(symbols: Iterable):
    """Takes an iterable of sympy symbols. Returns the symbols as a tuple, ordered
    by their string values. This can be used as a hash for a given collection of free symbols
    in an equation."""
    res = tuple(sorted(symbols, key=str))
    return res


class Formula:
    def __init__(self, func, result_class: ParameterMeta):
        self.func = func
        self.result_class = result_class
        self.unit = self.result_class.unit
        
        rhs = _convert_to_sympy_expression(func=self.func)
        lhs = sympy.Symbol(self.result_class.__name__)
        self.eq = sympy.Equality(lhs, rhs)
        
        wraps(func)(self)
        
        self.__doc__ = self._make_new_docstring()
    
    @property
    def inputs(self) -> tuple:
        """Returns a tuple representing the arguments for the underlying function."""
        res = tuple(inspect.signature(self.func).parameters.keys())
        return res
    
    def _make_new_docstring(self):
        """Constructs a new docstring which describes this object's __call__ method,
        the additional arguments introduced in it, and the docstring of the underlying functino."""
        
        argstring = ", ".join(self.inputs)
        parts = [
            f"Formula for {self.result_class.__name__}({argstring}).",
            self.__call__.__doc__,
            "Docstring for underlying function:",
            self.func.__doc__,
        ]
        
        res = "\n".join(parts)
        return res
    
    def validate_input(self, kwds: dict):
        """Ensures that either all or none of the inputs have declared types"""
        n_units = sum(has_units(val) for val in kwds.values())
        consistent = n_units in (0, len(kwds))
        if not consistent:
            raise RuntimeError(f"Units must be specified for all or no inputs (got {n_units}/{len(kwds)})")
    
    def validate_output(self, res):
        """Uses the result data model to validate the result of a computation"""
        self.result_class.validate(res)
    
    def _to_keywords(self, *args, **kwargs) -> dict:
        """Takes arbitrary args and kwargs and returns a dictionary representing both.
        This is done by inspecting the signature of the underlying function."""
        
        bound = inspect.signature(self.func).bind(*args, **kwargs)
        bound.apply_defaults()
        kw_only = bound.arguments
        return kw_only
    
    def compute(self, kwds: dict, with_units: bool|None):
        """If with_units is True or False, the result will have units/be unitless.
        If with_units is None (default), the result will have units only if all inputs do."""
        
        # Compute and run sanity checks
        self.validate_input(kwds=kwds)
        res = self.func(**kwds)
        self.validate_output(res=res)
        
        # Enforce units/no units if with_unit keyword is used.
        if with_units is True:
            res = self.result_class.ensure_units(value=res)
        elif with_units is False:
            res = self.result_class.ensure_float(value=res)
            
        return res
    
    def __call__(self, *args, with_units: bool=None, **kwargs):
        """Computes a quantity given the input values.
        If with_units is True or False, the result will have units/be unitless.
        If with_units is None (default), the result will have units only if all inputs do."""

        kwds = self._to_keywords(*args, **kwargs)        
        res = self.compute(kwds=kwds, with_units=with_units)
        
        return res
    
    def __repr__(self):
        s = format_sympy_equation(self.eq)
        return s
    #
    
    def __str__(self):
        res = f"Formula: {repr(self)}"
        return res
    #


if __name__ == '__main__':
    pass
