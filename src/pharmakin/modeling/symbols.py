import sympy
from sympy.core.function import AppliedUndef, Derivative, UndefinedFunction
from typing import Callable, Generic, Iterator, TypeVar

# Time variable
_time_symbol = "t"
t = sympy.Symbol(_time_symbol)

# Amount
_amount_symbol = "A"
A = sympy.Function(_amount_symbol)
A_t = A(t)
Ap = sympy.Derivative(A(t), t)


def _symbol_for_specific(label: str) -> str:
    res = f"{_amount_symbol}_{label}"
    return res


S = TypeVar('S')


class SymbolMaker(Generic[S]):
    """Handles the creation of often-needed sympy symbols from strings.
    Sympy uses memory address to identify symbols, so similar symbols defined in different places will not
    be interpreted as the same entity by sympy. This is desired in some situations, but not in others.
    This class serves as a housekeeping tool for such symbols. The symbol maker exposes methods for
    getting a symbol by a 'label', creating it if it hasn't been yet.
    Hence, one can use some_symbol_maker["x"] in different locations, and as long as the same symbol maker
    instance is used, the resulting variables will refer to the same symbolic entity."""

    def __init__(self, symbol_factory: Callable[[str], S], create_hook: Callable[[str], None]|None=None):
        """symbol_factory: callable which takes a string and returns a resulting sympy object
        create_hook: optional callable which is called whenever a new symbol with a given label is created.
            This can be used to ensure that for every such label, a number of different symbols have been defined"""
        
        self._cache: dict[str, S] = dict()
        self.f = symbol_factory
        self.create_hook = create_hook
    
    def _ensure_exists(self, label: str) -> None:
        """Checks if the label already maps to a sympy object, and creates the object if not"""
        if label not in self._cache:
            sym = self.f(label)
            self._cache[label] = sym
        #
    
    def __contains__(self, key: str) -> bool:
        return key in self._cache

    def ensure_exists(self, label: str) -> None:
        """Checks that the label is mapped to a sympy object. Calls the creation hook if one is set, alerting
        other objects ot the creation of a new sympy object."""

        if label in self:
            return
        
        self._ensure_exists(label)
        if self.create_hook is not None:
            self.create_hook(label)

    def __getitem__(self, key: str) -> S:
        self.ensure_exists(key)
        res = self._cache[key]
        return res
    #


class SymbolRegistry:
    """Registry for sympy objects which are useful for PK models.
    This holds a number of symbol makers which map a label (corresponding to some compound),
    to various symbols for concentration, rate of change, etc.
    The symbol makers are linked such that the same labels should be available in all symbol makers, i.e.
    as soon as one calls e.g.
        my_amount = sym.A["some_drug"]
    to get a symbol for the concentration A, symbols will immediately be available for related quantities,
    such as the concentration as a function of time A(t), in sym_A_t["some_drug"]."""

    t: sympy.Symbol = t  # Time symbol
    A: SymbolMaker[UndefinedFunction]  # Amount (concentration)
    A_t: SymbolMaker[AppliedUndef]  # Amount as function of time A(t)
    Ap: SymbolMaker[Derivative]  # Rate of change  A'(t)

    def __init__(self):
        d = dict(create_hook=self._check_symbol_exists)
        self.A = SymbolMaker(symbol_factory=self._make_A, **d)
        self.A_t = SymbolMaker(symbol_factory=self._make_A_t, **d)
        self.Ap = SymbolMaker(symbol_factory=self._make_Ap, **d)

    def _check_symbol_exists(self, label: str):
        
        symbol_makers = (a for a in vars(self).values() if isinstance(a, SymbolMaker))
        for symbol_maker in symbol_makers:
            symbol_maker._ensure_exists(label)
            #
        #

    @staticmethod
    def _symbol_with_suffix(label: str):
        return f"{_amount_symbol}_{label}"
    
    def _make_A(self, label: str) -> UndefinedFunction:
        varname = self._symbol_with_suffix(label)
        res = sympy.Function(varname)
        return res
    
    def _make_A_t(self, label: str) -> AppliedUndef:
        A = self.A[label]
        res = A(self.t)
        return res
    
    def _make_Ap(self, label:str) -> Derivative:
        res = sympy.Derivative(self.A_t[label], self.t)
        return res
    #
