from __future__ import annotations
import logging
logger = logging.getLogger(__name__)
import networkx as nx
import numpy as np
import sympy
from sympy import Function, dsolve, Derivative, Eq
from sympy.core.function import AppliedUndef
from typing import Any, Iterable, TypeAlias

from pharmakin.kinetics.first_order import k_el, t_half_from_k
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
        self.A_t = self.A(symbols.t)
        self.Ap = sympy.Derivative(self.A(symbols.t), symbols.t)


class Model:
    def __init__(self):
        self.G = nx.DiGraph()
        
        self.compounds: dict[str, Compound] = dict()
        # Map concentration(time) to compound objects
        self._A_t_to_compound: dict[sympy.core.function.AppliedUndef, Compound] = dict()
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
    def __init__(self, *rates: tuple[Any, sympy.Expr]):
        for gradient, expr in rates:
            pass  # TODO maybe do some type checking here?
        
        logger.debug(f"Created reaction: {rates}.")
        self.rates = rates
        #TODO implement stoichiometry stuff (taking differences in molecular mass into account)


def _summarize_reactions(reactions: Iterable[Reaction]) -> dict:
    res = dict()
    for reaction in reactions:
        for gradient, expr in reaction.rates:
            try:
                res[gradient] += expr
            except KeyError:
                res[gradient] = expr
            #
        #
    
    return res


def solve_numerical(compounds: Iterable[Compound], reactions: Iterable[Reaction], delta_t: float, T:float|int):
    
    # TODO make this cleaner and faster. Make test for prodrug system, then optimize this!!!
    compounds = sorted(compounds, key=lambda c: c.label)
    x = [c.A_0 for c in compounds]

    temp = [[val for val in x]]

    conc_vars = tuple(c.A_t for c in compounds)
    f_elems = []
    
    gradients = _summarize_reactions(reactions)
    for c in compounds:
        slope = gradients[c.Ap]
        f_dA = sympy.lambdify(conc_vars, slope)
        f_elems.append(f_dA)
    
    t = 0.0
    while t < T:
        dA_dt = [f(*x) for f in f_elems]
        x = [a+b*delta_t for a, b in zip(x, dA_dt)]
        temp.append(x)
        t += delta_t


    seqs = list(zip(*temp))
    res = {c.label: s for c, s in zip(compounds, seqs)}
    return res    


param: TypeAlias = float|int|sympy.Basic
solution_type: TypeAlias = dict[str, sympy.Expr]


class FirstOrderReaction(Reaction):
    def __init__(self, k: param, reactant: Compound, metabolite: Compound|None=None):

        self.k = k
        self.t_half = t_half_from_k(k)
        self.rate = self.k*reactant.A_t
        
        rates = [(Derivative(reactant.A_t, symbols.t), -self.rate)]
        if metabolite:
            rates.append((Derivative(metabolite.A_t, symbols.t), +self.rate))
        super().__init__(*rates)
    

class FirstOrderModel(Model):
    def add_reaction(self, reactant_label: str, metabolite_label: str|None=None, k: param|None=None, t_half: param=None):
        if k is None:
            k = k_el(t_half)
            
        self._ensure_added(reactant_label, metabolite_label)
        
        reactant = self.compounds[reactant_label]
        metabolite = self.compounds.get(metabolite_label)
        
        reaction = FirstOrderReaction(k=k, reactant=reactant, metabolite=metabolite)
        
        self.reactions.append(reaction)
        if metabolite:
            self.G.add_edge(reactant_label, metabolite_label, reaction=reaction)
        
        return self
    
    def get_initial_conditions(self):
        res = dict()
        for compound in self.compounds.values():
            res[compound.A(0)] = compound.A_0
        
        return res
    
    def get_equations(self):
        d = _summarize_reactions(self.reactions)
        res = []
        for gradient, expr in d.items():
            eq = sympy.Eq(gradient, expr)
            res.append(eq)
        
        return res

    def solve_analytic(self) -> dict[str, sympy.Expr]:
        eqs = self.get_equations()
        ics = self.get_initial_conditions()
        sols = dsolve(
            eqs,
            ics=ics
        )
        
        res = {self._A_t_to_compound[solution.lhs].label: solution.rhs for solution in sols}
        
        return res

    def solve_numerical(self, delta_t: float, T: float|int):
        res = solve_numerical(compounds=self.compounds.values(), reactions=self.reactions, delta_t=delta_t, T=T)
        return res


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    hmm = t_half_from_k(2.0)
    
    
    model = FirstOrderModel()
    model.add_compound("LDX", initial_amount=100.0)
    model.add_compound("AMP")
    model.add_reaction(reactant_label="LDX", metabolite_label="AMP", t_half=1.0)
    model.add_reaction(reactant_label="AMP", t_half=10.5)

    c = model.compounds["AMP"]
    
    a = model.solve_analytic()
    print(a)
    print()

    num = model.solve_numerical(delta_t=0.01, T=24)
    
    for k, v in num.items():
        print(k, v[:5])
    