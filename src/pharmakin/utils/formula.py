from __future__ import annotations
from collections import defaultdict
from functools import wraps
import inspect
import sympy
from typing import Callable, Iterable, TYPE_CHECKING

if TYPE_CHECKING:
    from pharmakin.utils.parameter_base import ParameterMeta
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
        
        self.rhs = _convert_to_sympy_expression(func=self.func)
        self.lhs = sympy.Symbol(self.result_class.__name__)
        self.eq = sympy.Equality(self.lhs, self.rhs)
        
        self._call = self._make_callable()
        
        wraps(func)(self)
        
        self.__doc__ = self._make_new_docstring()
    
    def _make_callable(self) -> Callable:
        """Returns a the input function, but with any sympy operations evaluated numerically"""
        expr = self.rhs
        if expr.atoms(sympy.Function):
            expr = expr.evalf()
        
        input_symols = [sympy.Symbol(x) for x in self.inputs]
        f = sympy.lambdify(input_symols, expr)
        return f
    
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
        
        res = "\n".join([part for part in parts if part is not None])
        return res
    
    def is_symbolic(self) -> bool:
        """Returns a bool indicating whether the formula contains symbolic operations (like sympy.log).
        We need to check for this because pint.Units don't play nice with such operations."""
        res = self.rhs.atoms(sympy.Function)
        return res
    
    @staticmethod
    def validate_input(kwds: dict):
        """Ensures that either all or none of the inputs have declared units"""
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
    
    def compute(self, kwds: dict, with_units: bool|None, evaluate_numerically=False):
        """If with_units is True or False, the result will have units/be unitless.
        If with_units is None (default), the result will have units only if all inputs do.
        evaluate_numerically indicates whether sympy expressions should be evaluated numerically."""
        
        # Compute and run sanity checks
        self.validate_input(kwds=kwds)
        res = self._call(**kwds)
        self.validate_output(res=res)
        
        # Enforce units/no units if with_unit keyword is used.
        if with_units is True:
            res = self.result_class.ensure_units(value=res)
        elif with_units is False:
            res = self.result_class.ensure_float(value=res)
    
        if evaluate_numerically and isinstance(res, sympy.core.expr.Expr):
            res = res.evalf()

        return res
    
    def __call__(self, *args, with_units: bool=None, **kwargs):
        """Computes a quantity given the input values.
        If with_units is True or False, the result will have units/be unitless.
        If with_units is None (default), the result will have units only if all inputs do."""

        kwds = self._to_keywords(*args, **kwargs)
        res = self.compute(kwds=kwds, with_units=with_units, evaluate_numerically=True)
        
        return res
    
    def __repr__(self):
        s = format_sympy_equation(self.eq)
        return s
    #
    
    def __str__(self):
        res = f"Formula: {repr(self)}"
        return res
    #


class Formulary:
    """Helper class for keeping track of a collection of pharmacokinetic parameters,
    and relations between them."""

    def __init__(self, parameters: Iterable[ParameterMeta], formulas: Iterable[Formula]):
        # For storing callable formulas and sympy equations
        self.formulae = []
        self.equations = []

        # Maps from parameter names to parameter classes and sympy symbols
        self.parameters = dict()
        self.symbols = dict()
        # Map parameter names
        self.formulas_for_parameter = defaultdict(list)
        self.formulas_containing_parameter = defaultdict(list)

        # Register the parameters + formulas
        self.register_parameters(parameters)
        self.register_formulas(formulas)

    @staticmethod
    def _parkey(param: ParameterMeta) -> str:
        """Creates a name (string) from a parameter class"""
        return param.__name__
    
    def get_formulas_containing_parameter(self, parameter: ParameterMeta) -> list:
        k = self._parkey(parameter)
        res = self.formulas_containing_parameter[k]
        return res

    def _is_registered(self, param: ParameterMeta|str) -> bool:
        """Denotes whether a parameter has been registered. Works on classes and strings (parameter names)"""
        k = param if isinstance(param, str) else self._parkey(param=param)
        return k in self.parameters

    def ensure_registered(self, pars: Iterable[ParameterMeta|str]):
        failed = set([])
        for par in pars:
            k = par if isinstance(par, str) else self._parkey(par)
            if k not in self.parameters:
                failed.add(k)
            #

        if failed:
            msg = f"The following parameters were not registered: {', '.join(sorted(failed))}."
            raise RuntimeError(msg)

    def register_parameters(self, params: Iterable[ParameterMeta]):
        """Registers parameter classes to the formulary"""
        for param in params:
            # Make sure parameters are only registered once
            if self._is_registered(param):
                raise RuntimeError(f"Parameter {param} has already been registered")

            # Create a name and sympy symbol for the parameter
            k = self._parkey(param=param)
            self.parameters[k] = param
            self.symbols[k] = sympy.Symbol(k)
        #

    def register_formulas(self, formulas: Iterable[Formula]):
        """Registers formula instances to the formulary"""
        for formula in formulas:
            self.formulae.append(formula)
            self.equations.append(formula.eq)

            formula_for = self._parkey(formula.result_class)
            vars = list(formula.inputs) + [formula_for]
            self.ensure_registered(vars)

            parkey = self._parkey(formula.result_class)
            self.formulas_for_parameter[parkey].append(formula)

            for var in vars:
                self.formulas_containing_parameter[var].append(formula)
            #
        #

    def determine_parameter(self, parameter: ParameterMeta, **kwargs):
        """Attempts to determine a parameter value from other parameters, specified as keyword arguments.
        The first argument is the target parameter class, which should be computed.
        The remaining keyword arguments are key-values for other parameters.
        The formulary then attempts to solve its equations for the target parameter in terms of
        the provided parameter values."""

        Formula.validate_input(kwargs)

        # Check if inputs have units
        res_units = all(has_units(val) for val in kwargs.values())

        # Determine target sympy symbol
        target_name = self._parkey(parameter)
        target = self.symbols[target_name]

        # Use keyword args to make a dict of other symbols and their values
        vals = dict()
        for otherpar, val in kwargs.items():
            k = self.symbols[otherpar]
            class_ = self.parameters[otherpar]
            # Drop units to avoid sympy issues
            val = class_.ensure_float(val)
            vals[k] = val

        # Insert into the formulary's equations and solve
        substituted = [eq.subs(vals) for eq in self.equations]
        solved = sympy.solve(substituted, target)
        res = solved[target]

        # If input had units, add to result
        if res_units:
            res = parameter.ensure_units(res)

        parameter.validate(res)

        return res

    def __repr__(self):
        s = f"Formulary ({len(self.parameters)} parameters, {len(self.formulae)} formulas)"
        return s

    def __str__(self):
        s = f"Formulary - parameters: {', '.join(sorted(self.parameters.keys()))}"
        for formula in self.formulae:
            s += "\n"+"  " + repr(formula)
        return s
    #


if __name__ == '__main__':
    pass
