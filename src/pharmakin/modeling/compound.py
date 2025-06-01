import sympy
from sympy.core.function import UndefinedFunction
from typing import cast

from pharmakin.modeling import symbols


class Compound:
    def __init__(self, label: str, A_0: float=0.0):
        """Make a new quantity for modeling.
        label (str) - A label/name to describe the compound.
            This is just for readibality, so you can use the full medicine name (e.g. lisdexamphetamine),
            Abreviation (LDX), brand name (Elvanse/Vyvanse) or whatever.
        A_0 (float, default 0.0): The initial quantity of the drug."""

        self.label = label
        self.A_0 = A_0
        # List of decays, indicating the rates and resulting compounds to which the drug metabolizes

        A_suffix_str = f"{str(symbols.A)}_{self.label}"
        self.A = sympy.Function(A_suffix_str)
        self.A = cast(UndefinedFunction, self.A)
        self.A_t = self.A(symbols.t)
        self.Ap = sympy.Derivative(self.A(symbols.t), symbols.t)
    #


def get_initial_conditions(*compounds: Compound) -> dict[sympy.core.function.AppliedUndef,int|float]:
    res = {c.A(0): c.A_0 for c in compounds}
    return res
