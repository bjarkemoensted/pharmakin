import logging
logger = logging.getLogger(__name__)
import sympy


class Reaction:
    def __init__(self, *rates: tuple[sympy.Derivative, sympy.Expr]):
        logger.debug(f"Created reaction: {rates}.")
        self.rates = rates


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