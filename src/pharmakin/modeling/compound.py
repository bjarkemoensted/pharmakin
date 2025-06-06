from __future__ import annotations
from dataclasses import dataclass
import sympy
from sympy.core.function import AppliedUndef, Derivative, UndefinedFunction
from sympy.core.function import UndefinedFunction

from pharmakin.modeling.symbols import SymbolRegistry


_default_initial_amount = 0.0


# TODO would be cleaner to just use compounds as a dataclass and leave it to models/solvers to extract symbols etc!!!
@dataclass
class Compound:
    """A new quantity for modeling.
    label (str) - A label/name to describe the compound.
        This is just for readability, so you can use the full medicine name (e.g. lisdexamphetamine),
        Abreviation (LDX), brand name (Elvanse/Vyvanse) or whatever.
    A_0 (float, default 0.0): The initial quantity of the drug.
    A: A sympy function representing the amount
    A_t Sympy function of time, representing amount as function of time
    Ap: derivative of amount wrt. time"""

    label: str
    A: UndefinedFunction
    A_t: AppliedUndef
    Ap: Derivative
    A_0: float=_default_initial_amount

    @classmethod
    def create_with_symbol_registry(cls, symbol_registry: SymbolRegistry, label: str, A_0: float=_default_initial_amount) -> Compound:
        A = symbol_registry.A[label]
        A_t = symbol_registry.A_t[label]
        Ap = symbol_registry.Ap[label]

        res = cls(label=label, A_0=A_0, A=A, A_t=A_t, Ap=Ap)
        return res
    #


def get_initial_conditions(*compounds: Compound) -> dict[sympy.core.function.AppliedUndef,int|float]:
    res = {c.A(0): c.A_0 for c in compounds}
    return res
