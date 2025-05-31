from __future__ import annotations
import logging
logger = logging.getLogger(__name__)
import networkx as nx
import numpy as np
import sympy
from sympy import Derivative
from sympy.core.function import AppliedUndef, UndefinedFunction
from typing import Any, Callable, cast, Iterable, Literal, TypeAlias

from pharmakin.kinetics.first_order import k_el, t_half_from_k
from pharmakin.modeling import symbols
from pharmakin.modeling.solvers import Solver


solve_methods: TypeAlias = Literal["analytic", "numeric", "simple", "auto"]


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


class Model:
    def __init__(self) -> None:
        self.G: nx.DiGraph = nx.DiGraph()
        
        self.compounds: dict[str, Compound] = dict()
        # Map concentration(time) to compound objects
        self._A_t_to_compound: dict[AppliedUndef, Compound] = dict()
        self.reactions: list[Reaction] = []

    def add_compound(self, label: str, initial_amount = 0.0):
        if label in self.compounds:
            raise RuntimeError(f"Compound {label} has already been added to model.")
        
        compound = Compound(label=label, A_0=initial_amount)
        self.G.add_node(label)
        self.compounds[label] = compound
        self._A_t_to_compound[compound.A_t] = compound
        return self
    
    def _ensure_added(self, *labels: str|None):
        for label in labels:
            if label is not None and label not in self.compounds:
                self.add_compound(label)


class Reaction:
    def __init__(self, *rates: tuple[sympy.Derivative, sympy.Expr]):
        logger.debug(f"Created reaction: {rates}.")
        self.rates = rates
    #

def _summarize_reactions(reactions: Iterable[Reaction]) -> dict[sympy.Derivative, sympy.Expr]:
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


param: TypeAlias = float|int|sympy.Basic
solution_type: TypeAlias = sympy.Expr|Callable[[float], float]


class FirstOrderReaction(Reaction):
    def __init__(self, k: param, reactant: Compound, metabolite: Compound|None=None, ratio: float=1.0):
        """Reaction, first order kinetics
        k: elimination constant (function of half-life)
        reactant: Compound - the compound which is converted in the reaction.
        metabolite: Compound (optional, defaults to None). The resulting metabolite.
            Can be left as None to not include any metabolite in the reaction.
            This can be done if the metabolite is not of interest in the model, for example
            if the reaction represents a drug broken down into components without any
            pharmacological effect.
        ratio: float (default 1.0). The metabolite to reactant mass ratio.
            This can be used for reactions where the reactant is not converted to the metabolite
            in a 1:1 mass ratio.
        """
        
        self.k = k
        self.t_half = t_half_from_k(k)
        self.ratio = ratio
        self.rate = self.k*reactant.A_t
        
        rates = [(Derivative(reactant.A_t, symbols.t), -self.rate)]
        if metabolite:
            rates.append((Derivative(metabolite.A_t, symbols.t), +self.rate*ratio))
        super().__init__(*rates)
    

class FirstOrderModel(Model):
    def add_reaction(
            self,
            reactant_label: str,
            metabolite_label: str|None=None,
            k: param|None=None,
            t_half: param=None,
            ratio: float=1.0):
        """Adds a reaction to the model.
        reactant_label: Label for the reactant (drug which is converted into another)
        metabolite_label: Label for metabolite (the resulting compound). Can be None to ignore/represent excretion.
        k, t_half - elimination constant and half_life. One must be specified. Defines the speeed of the reaction,
            as a function of the concentrations of compounds involved in the reaction.
        ratio: The ratio by which to multiply the resulting compound concentration. Defaults to 1.0.
            The purpose of this parameter is to take stoichiometry into account, so compounds aren't
            necessarily converted in a 1:1 mass ratio."""

        if k is None:
            k = k_el(t_half)
            
        self._ensure_added(reactant_label, metabolite_label)
        
        reactant = self.compounds[reactant_label]
        metabolite = self.compounds[metabolite_label] if metabolite_label is not None else None
        
        reaction = FirstOrderReaction(k=k, reactant=reactant, metabolite=metabolite, ratio=ratio)
        
        self.reactions.append(reaction)
        if metabolite:
            self.G.add_edge(reactant_label, metabolite_label, reaction=reaction)
        
        return self
    
    def get_initial_conditions(self) -> dict[sympy.core.function.AppliedUndef,int|float]:
        """Get a dict containing the initial conditions for the system."""
        res = dict()
        for compound in self.compounds.values():
            res[compound.A(0)] = compound.A_0
        
        return res
    
    def get_equations(self) -> list[sympy.Eq]:
        """Get a list of equations describing the temporal dynamics of the system."""
        d = _summarize_reactions(self.reactions)
        res = []
        for gradient, expr in d.items():
            eq = sympy.Eq(gradient, expr)
            res.append(eq)
        
        return res

    def _make_solver(self, t_vals: np.ndarray) -> Solver:
        """Sets up a solver for solving the system"""
        funcs = [c.A_t for c in self.compounds.values()]
        gradients = _summarize_reactions(self.reactions)
        ics = self.get_initial_conditions()
        solver = Solver(
            funcs = funcs,
            gradients = gradients,
            ics = ics,
            t_vals=t_vals
        )
        return solver
    
    def solve(self, t_vals: np.ndarray, how: solve_methods="auto", **kwargs):
        """Solves the system for the specified time values.
        t_vals: numpy array representing time values.
        how: The method to be used when solving. Can be:
            "analytic": Uses sympy to solve the ODEs symbolically.
            "numeric": Uses scipy's ivp method for numerical intergration.
            "simple": Uses my own (super inefficient probably) solver which repeatedly iterates to the next time step
                by using a second-order approximation given the current values.
            "auto": Attempts the aforementioned methods in the listed order.
        **kwargs will be forwarded to the corresponding solve method in the solver."""
        
        # Attempt multiple solve methods if auto is selected
        if how == "auto":
            priorities = ("analytic", "numeric", "simple")
            for method in priorities:
                try:
                    return self.solve(t_vals=t_vals, how=cast(solve_methods, method), **kwargs)
                except (ValueError, IndexError):
                    continue
                #
            raise RuntimeError(f"All methods {priorities} failed to reach a solution.")
        
        # Make a solver and select solve method corresponding to the 'how' parameter
        solver = self._make_solver(t_vals=t_vals)
        methods: dict[solve_methods, Callable] = dict(
            analytic=solver.solve_analytic,
            simple=solver.solve_simple,
            numeric=solver.solve_numerical
        )
        
        try:
            solve_func = methods[how]
        except KeyError:
            raise ValueError(f"Invalid solve method: '{how}'")
    
        # Solve the system
        solutions = solve_func(**kwargs)
        
        # Use compound labels as key, instead of the symbol for concentration(time) (A(t)).
        res = {self._A_t_to_compound[k].label: v for k, v in solutions.items()}
        return res
    #


if __name__ == '__main__':
    pass
