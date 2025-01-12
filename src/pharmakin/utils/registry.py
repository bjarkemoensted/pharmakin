import sympy
from typing import Iterable

from pharmakin.utils.parameter_base import ParameterMeta
from pharmakin.utils.formulas import Formula
from pharmakin.utils.units import has_units


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
        
        # Register the parameters + formulas
        self.register_parameters(parameters)
        self.register_formulas(formulas)        

    def _parkey(self, param: ParameterMeta) -> str:
        """Creates a name (string) from a parameter class"""
        return param.__name__

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
            
            vars = list(formula.inputs) + [formula.result_class]
            self.ensure_registered(vars)
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


if __name__ == '__main__':
    pass
