from __future__ import annotations
from collections import defaultdict
import sympy

from typing import Callable, Iterable, TYPE_CHECKING

if TYPE_CHECKING:
    from pharmakin.utils.parameter_base import ParameterMeta
from pharmakin.utils.formulas import Formula



class Formulary:
    """Helper class for keeping track of a collection of pharmacokinetic parameters,
    and relations between them."""

    def __init__(self):
        self.formulae = []
        self.parameters = dict()
        self.equations = []
        self.symbols = dict()
        
        
        self._var2parclass = dict()
        self.symbol_to_eqs = defaultdict(lambda: [])
        self.free_symbols_to_eqs = defaultdict(lambda: [])

    def _parkey(self, param: ParameterMeta):
        return param.__name__

    def _is_registered(self, param: ParameterMeta|str):
        if isinstance(param, ParameterMeta):
            k = self._parkey(param=param)

    def _register_parameter(self, param: ParameterMeta):
        k = self._parkey(param=param)
        
        # TODO FIX!!!
        #_ = self._is_registered(param)
        
        if k in self.parameters:
            raise RuntimeError(f"Parameter {param} has already been registered")
        
        self.parameters[k] = param
        self.symbols[k] = sympy.Symbol(k)

    def register_parameters(self, params: Iterable[ParameterMeta]):
        for param in params:
            self._register_parameter(param=param)
        #

    def _register_formula(self, formula: Formula):
        self.formulae.append(formula)
        self.equations.append(formula.eq)
        print(formula.inputs)

    def register_formulas(self, formulas: Iterable[Formula]):
        for formula in formulas:
            self._register_formula(formula=formula)
        #
    
    def determine_parameter(self, parameter: ParameterMeta, **kwargs):
        pass
    
    def __str__(self):
        s = f"Formulary - parameters: {', '.join(sorted(self.parameters.keys()))}"
        for formula in self.formulae:
            s += "\n"+"  " + repr(formula)
        return s


if __name__ == '__main__':
    pass
