import logging
logger = logging.getLogger(__name__)
import sympy
from typing import TypeAlias

from pharmakin.modeling.compound import Compound
from pharmakin.modeling import symbols



class Reaction:
    def __init__(self, *rates: tuple[sympy.Derivative, sympy.Expr], t_start:float=0.0, t_stop: float=float('inf')) -> None:
        logger.debug(f"Created reaction: {rates}.")
        self.rates = rates
        self.t_start = t_start
        self.t_stop = t_stop
    #    


class Instant(Reaction):
    # TODO get this working with solvers!!!
    
    #TODO also!!! https://chatgpt.com/c/684155a9-7aac-8001-8215-3f2e57e58e8e
    def __init__(self, gradient: sympy.Derivative, amount: float, t0:float=0.0, **kwargs):
        r = amount*sympy.DiracDelta(symbols.t - t0)
        super().__init__((gradient, r), t_start=t0, t_stop=t0, **kwargs)


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


class ReactionSummay:
    def __init__(self, *reactions: Reaction):
        pass


if __name__ == '__main__':
    from pharmakin.modeling import example_models
    ex = example_models.dexamphetamine_example()
    m = ex.model
    
    r = m.reactions[0]
    c = m.compounds["AMP"]
    print(c.Ap)
    r = Instant(gradient=c.Ap, t0=1.0, amount=60.0)
    print(r, r.rates[0])
