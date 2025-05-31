""" Contains some example models, for e.g. comparing if solutions match expected solutions.
Analyitcal first=order kinetics solutions for a single drug, and for a prodrug model are taken from:
Cho S, Yoon YR. Understanding the pharmacokinetics of prodrug and metabolite.
Translational and Clinical Pharmacology. 2018"""

from __future__ import annotations
from dataclasses import dataclass
import sympy
from typing import Callable

from pharmakin.kinetics.first_order import k_el
from pharmakin.modeling.base import Model, FirstOrderModel
from pharmakin.modeling import symbols


# Define half lifes and labels for a few compounds
AMP_T_HALF = 10.5
AMP_LABEL = "AMP"
LDX_T_HALF = 1.0
LDX_LABEL = "LDX"


def _lambdify_solutions(solutions: dict[str, sympy.Expr]) -> dict[str, Callable[[float], float]]:
    """Helper function for turning symbolic (sympy) solutions into callables, which take a numpy array representing
    time. Both input and result are dictionaries mapping strings (labels representing a shorthand name
    for the compound) to the solution for that compound as a function of time."""
    res = {k: sympy.lambdify(symbols.t, v) for k, v in solutions.items()}
    return res


@dataclass
class Example:
    """Helper class for exposing some example/toy models, along with established analytical
    solutions from the literature. The idea is to provide a simple way of instantiating common models,
    both for rapid prototyping of e.g. visualization, and for simplifying tests."""
    
    model: Model
    solution: dict[str, sympy.Expr]

    @classmethod
    def linear_single_drug(cls, t_half: float|int, initial_amount: float=100.0, label: str="drug") -> Example:
        """Sets up a model for a single drug being metabolized (with first=order kinetics).
        Returns the model and the analyitical solution. For consistency with the prodrug model,
        the solution is provided as a list of a single sympy expression, representing the
        drug concentration as a function of time."""
        
        k = k_el(half_life=t_half)
        
        # Set up the model
        model = FirstOrderModel().add_compound(
            label, initial_amount=initial_amount
        ).add_reaction(label, k=k)
        
        # Make the analytical solution
        solution_symbolic = {label: initial_amount*sympy.exp(-k*symbols.t)}
        solution = _lambdify_solutions(solution_symbolic)
        res = cls(model=model, solution=solution)
        return res
    
    @classmethod
    def linear_prodrug(
        cls,
        t_half_prodrug: float|int,
        t_half_active: float|int,
        initial_amount: float=100.0,
        label_prodrug: str="prodrug",
        label_active: str="active"
    ) -> Example:
        """Makes model and solution for a prodrug system, with a prodrug metabolized into an active drug."""
        
        k1 = k_el(half_life=t_half_prodrug)
        k2 = k_el(half_life=t_half_active)
        
        # Set up the model
        model = FirstOrderModel().add_compound(
            label_prodrug, initial_amount=initial_amount
        ).add_compound(
            label_active, initial_amount=0.0
        ).add_reaction(
            label_prodrug, label_active, k=k1
        ).add_reaction(
            label_active, k=k2
        )
        
        # Construct analytic solutions for prodrug and active drug, respectively
        solution_symbolic = dict()
        solution_symbolic[label_prodrug] = initial_amount*sympy.exp(-k1*symbols.t)
        if k1 != k2:
            norm = (k1*initial_amount)/(k2 - k1)
            solution_symbolic[label_active] = norm*(sympy.exp(-k1*symbols.t) - sympy.exp(-k2*symbols.t))
        else:
            # Special case when the rates are identical (think this can be obtained via L'Hopital from the general case)
            solution_symbolic[label_active] = k1*initial_amount*symbols.t*sympy.exp(-k1*symbols.t)
            
        solution = _lambdify_solutions(solution_symbolic)
        res = cls(model=model, solution=solution)
        return res


def dexamphetamine_example(
        t_half: float=AMP_T_HALF,
        initial_amount=20.0,
        label=AMP_LABEL
    ) -> Example:
    """Return example model for dexamphetamine"""
    res = Example.linear_single_drug(t_half=t_half, initial_amount=initial_amount, label=label)
    return res


def lisdexamphetamine_example(
        t_half_prodrug: float=LDX_T_HALF,
        t_half_active: float=AMP_T_HALF,
        initial_amount: float=60.0,
        label_prodrug: str=LDX_LABEL,
        label_active: str=AMP_LABEL
    ) -> Example:
        """Return example model of lisdexamphetamine (i.e. brand names Vyvanse/Elvanse)."""
        
        res = Example.linear_prodrug(
            t_half_active=t_half_active,
            t_half_prodrug=t_half_prodrug,
            initial_amount=initial_amount,
            label_active=label_active,
            label_prodrug=label_prodrug
        )
        
        return res
    #
