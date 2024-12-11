from collections import defaultdict
import sympy
from typing import Iterable


class Registry(dict):
    """Registry class for uniquely registering quantities under a given name.
    Similar to a dict, but throws an error if a name has already been registered."""

    def __setitem__(self, key, value):
        if key in self:
            raise RuntimeError(f"Key {key} is already set (to {value})")
        return super().__setitem__(key, value)
    #


class Formulary:
    def __init__(self):
        self.formulae = []
        self.symbol_to_eqs = defaultdict(lambda: [])
        self.free_symbols_to_eqs = defaultdict(lambda: [])
    
    @staticmethod
    def _format_equation(eq: sympy.core.relational.Equality) -> str:
        if not isinstance(eq, sympy.core.relational.Equality):
            raise TypeError
        
        res = f"{eq.lhs} = {eq.rhs}"
        return res
    
    @staticmethod
    def _key(symbols: Iterable):
        res = tuple(sorted(symbols, key=str))
        return res
    
    def lookup_free_symbols(self, symbols: Iterable):
        """Looks up equations with the specified free symbols"""
        k = self._key(symbols=symbols)
        for eq in self.free_symbols_to_eqs[k]:
            yield eq
    
    def register(self, eq: sympy.core.relational.Equality):
        """Register an equation to the formulary"""
        
        if not isinstance(eq, sympy.core.relational.Equality):
            raise TypeError
        
        self.formulae.append(eq)
        
        # Map each free symbol in the equation to it
        symbols = eq.free_symbols
        for symbol in symbols:
            self.symbol_to_eqs[symbol].append(eq)
        
        # Map the collection of free symbols to the equation as well
        symbols_key = self._key(symbols=symbols)
        self.free_symbols_to_eqs[symbols_key].append(eq)
    
    def __repr__(self):
        n = len(self.formulae)
        res = f"{self.__class__.__name__} (holding {n} equation{'s'*int(n != 1)})"
        return res
    
    def __str__(self):    
        header = f"*{repr(self)}*"
        eqs = [self._format_equation(eq) for eq in self.formulae]
        lines = [header] + ["  "+s for s in eqs]
        res = "\n".join(lines)
        return res
    #


if __name__ == '__main__':
    pass
