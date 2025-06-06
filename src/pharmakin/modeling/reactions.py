from __future__ import annotations
import logging
logger = logging.getLogger(__name__)
import sympy
from typing import Callable, get_args, Literal, TypeAlias

from pharmakin.modeling.compound import Compound
from pharmakin.modeling import symbols


class Rate:
    # TODO FIX!!!!!
    # Think the clean solution here is to define different types of rates, and include class-level sttributes
    # indicating how they may be aggregated by solvers.
    # For instance, an instantaneous release (Dirac Delta) can be handled by analytical solvers, but not standard
    # numerical ones (at least not mine).
    # Migh make sense to have e.g. a 'can_be_lambdified' flag, so numerical solvers can treat non-lambdifiable
    # rates separately. Those rates can then maybe expose methods for how A(t) changes when integrating across
    # a time region (Dirac Deltas bumping the amount iff we integrate across t0, for instance).
    
    def __init__(self, gradient: sympy.Derivative, expr: sympy.Expr, t_start:float=0.0, t_stop: float=float('inf')):
        self.gradient = gradient
        self.expr = expr
        self.t_start = t_start
        self.t_stop = t_stop


class Reaction:
    def __init__(self, *rates: tuple[sympy.Derivative, sympy.Expr], t_start:float=0.0, t_stop: float=float('inf')) -> None:
        logger.debug(f"Created reaction: {rates}.")
        self.rates = rates
        self.t_start = t_start
        self.t_stop = t_stop
    
    def is_active_at_time(self, t: float) -> bool:
        return self.t_start <= t < self.t_stop


class Administration(Reaction):
    # TODO get this working with solvers!!!
    
    #TODO also!!! https://chatgpt.com/c/684155a9-7aac-8001-8215-3f2e57e58e8e
    @classmethod
    def instantaneous(cls, gradient: sympy.Derivative, amount: float, t0:float=0.0) -> Administration:
        expr = amount*sympy.DiracDelta(symbols.t - t0)
        rate = (gradient, expr)
        return cls(rate, t_start=t0, t_stop=t0)
    #
        

def summarize_reactions(*reactions: Reaction) -> dict[sympy.Derivative, sympy.Expr]:
    res: dict[sympy.Derivative, sympy.Expr] = dict()
    for reaction in reactions:
        for gradient, expr in reaction.rates:
            try:
                res[gradient] += expr
            except KeyError:
                res[gradient] = expr
            #
        #

    return res


class ReactionSummary:
    def __init__(self, *reactions: Reaction):
        self.reactions: list[Reaction] = []
        self.critical_times = set([])
        for r in reactions:
            self.critical_times.add(r.t_start)
            self.critical_times.add(r.t_stop)
            self.reactions.append(r)
        #
    
    def get_active(self, t: float) -> list[Reaction]:
        res = [r for r in self.reactions if r.is_active_at_time(t)]
        return res
        

administration_type: TypeAlias = Literal["instantaneous"]
valid_administration_methods: tuple[administration_type, ...] = get_args(administration_type)


if __name__ == '__main__':
    from pharmakin.modeling import example_models
    from pharmakin.modeling.base import FirstOrderModel
    import typing
    ex = example_models.dexamphetamine_example()
    m = typing.cast(FirstOrderModel, ex.model)
    
    r = m.reactions[0]
    c = m.compounds["AMP"]
    print(c.Ap)
    
    m.add_administration("AMP", amount=60.0, kind="instantaneous")
    m.reactions.append(r)
    
    rs = ReactionSummary(*m.reactions)
    print(rs)
    
    print(rs.get_active(0.7))
    import numpy as np
    t = np.linspace(0.0, 100.0, num=1000)
    sol = m.solve(how="simple", t_vals=t)
    print(sol)