from __future__ import annotations
import logging

from pharmakin.modeling.compound import Compound
from pharmakin.modeling.reactions import Reaction, summarize_reactions, Administration
logger = logging.getLogger(__name__)
import networkx as nx
import numpy as np
import sympy
from sympy import Derivative
from sympy.core.function import AppliedUndef
from typing import cast, get_args, Literal, Type, TypeAlias

from pharmakin.modeling.reactions import administration_type
from pharmakin.kinetics.first_order import k_el, t_half_from_k
from pharmakin.modeling import symbols
from pharmakin.modeling.symbols import SymbolRegistry
from pharmakin.modeling import solvers


# Valid arguments for the 'how' keyword, indicating how to solve ODEs
solve_methods: TypeAlias = Literal["analytic", "numeric", "simple"]
# Solve methods plus an "auto" options
valid_solve_args: TypeAlias = solve_methods|Literal["auto"]

# Get valid solvers and map their names to the solver classes
valid_solve_methods: tuple[solve_methods, ...] = get_args(solve_methods)
_solvers: dict[str, Type[solvers.Solver]] = {
    "analytic": solvers.AnalyticalSolver,
    "numeric": solvers.NumericSolver,
    "simple": solvers.SimpleSolver
}


class Model:
    def __init__(self) -> None:
        self.G: nx.DiGraph = nx.DiGraph()

        # TODO use symbol registry instead of compounds to lookup A_t, Ap, etc  # !!!
        self.symbols = SymbolRegistry()
        self.compounds: dict[str, Compound] = dict()
        # Map concentration(time) to compound objects
        self._A_t_to_compound: dict[AppliedUndef, Compound] = dict()
        self.reactions: list[Reaction] = []

    def _lookup_A_t(self, label: str) -> AppliedUndef:
        c = self.compounds[label]
        res = c.A_t
        return res
    
    def _lookup_Ap(self, label: str) -> Derivative:
        c = self.compounds[label]
        res = c.Ap
        return res

    def add_compound(self, label: str, initial_amount = 0.0):
        if label in self.compounds:
            raise RuntimeError(f"Compound {label} has already been added to model.")
        
        compound = Compound.create_with_symbol_registry(
            label=label,
            A_0=initial_amount,
            symbol_registry=self.symbols
            )
        
        self.G.add_node(label)
        self.compounds[label] = compound
        self._A_t_to_compound[compound.A_t] = compound
        return self
    
    def _ensure_added(self, *labels: str|None):
        for label in labels:
            if label is not None and label not in self.compounds:
                self.add_compound(label)
            #
        #
    
    def solve(self, t_vals: np.ndarray, how: valid_solve_args="auto", **kwargs):
        raise NotImplementedError


param: TypeAlias = float|int|sympy.Basic


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
    
    def add_administration(self, label: str, kind: administration_type, amount: float, **kwargs):
        gradient = self.symbols.Ap[label]
        d = dict(gradient=gradient, amount=amount, **kwargs)
        match kind:
            case "instantaneous":
                adm = Administration.instantaneous(**d)
            case _:
                raise ValueError(f"Unrecognized drug administration: {kind}")
            #
        
        self.reactions.append(adm)
    
    def get_initial_conditions(self) -> dict[sympy.core.function.AppliedUndef,int|float]:
        """Get a dict containing the initial conditions for the system."""
        res = dict()
        for compound in self.compounds.values():
            res[compound.A(0)] = compound.A_0
        
        return res
    
    def get_equations(self) -> list[sympy.Eq]:
        """Get a list of equations describing the temporal dynamics of the system."""
        d = summarize_reactions(*self.reactions)
        res = []
        for gradient, expr in d.items():
            eq = sympy.Eq(gradient, expr)
            res.append(eq)
        
        return res

    def _make_solver(self, how: solve_methods) -> solvers.Solver:
        """Sets up a solver for solving the system"""
        
        if how not in _solvers:
            raise ValueError(f"Invalid solve method: '{how}'")
        
        class_ = _solvers[how]
        res = class_(compounds=self.compounds.values(), reactions=self.reactions)
        
        return res
    
    def solve(self, t_vals: np.ndarray, how: valid_solve_args="auto", **kwargs):
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
            priorities = valid_solve_methods
            for method in priorities:
                try:
                    return self.solve(t_vals=t_vals, how=cast(solve_methods, method), **kwargs)
                except (ValueError, IndexError):
                    continue
                #
            raise RuntimeError(f"All methods {priorities} failed to reach a solution.")
        
        # Solve the system
        solver = self._make_solver(how=how)
        solutions = solver.solve(t=t_vals)
        
        # Use compound labels as key, instead of the symbol for concentration(time) (A(t)).
        res = {self._A_t_to_compound[k].label: v for k, v in solutions.items()}
        return res
    #


if __name__ == '__main__':
    pass
