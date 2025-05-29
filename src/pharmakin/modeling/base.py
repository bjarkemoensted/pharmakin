from __future__ import annotations
import logging
logger = logging.getLogger(__name__)
import networkx as nx
import numpy as np
import sympy
from sympy import Function, dsolve, Derivative, Eq
from sympy.core.function import AppliedUndef
from typing import Any, Callable, Iterable, TypeAlias

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


# TODO make this cleaner and faster. Make test for prodrug system, then optimize this!!!
def solve_numerical(compounds: Iterable[Compound], reactions: Iterable[Reaction], t_vals: np.ndarray):
    
    # Check we start at t=0 (because we initialize with the initial values A_0 for each compound)
    assert t_vals[0] == 0.0
    
    # Make sure compounds are ordered
    compounds = sorted(compounds, key=lambda c: c.label)
    
    # Determine the rate of change in concentration of each compound (dA/dt), as a function of current concentrations
    conc_vars = tuple(c.A_t for c in compounds)
    gradients = _summarize_reactions(reactions)
    slopes = [gradients[c.Ap] for c in compounds]
    fp = sympy.lambdify(conc_vars, slopes, modules="numpy")
    
    # Wrap in a function which takes and returns a numpy array (each element corresponding to a concentration)
    def f_prime(s: np.array) -> np.array:
        res_list = fp(*s)
        return np.array(res_list)
    
    # Make a matrix for holding the results, and set the first row to initial values
    x = np.array([c.A_0 for c in compounds])
    m = np.empty(shape=(len(t_vals), len(compounds)))
    m.fill(np.nan)
    ind = 0
    m[ind, :] = x
    
    dA_dt_prev = np.full(shape=x.shape, fill_value=np.nan)
    
    # Fill up the matrix using the gradients to interpolate from previous points
    for dt in np.diff(t_vals):
        ind += 1
        dA_dt = f_prime(x)
        step = dt*dA_dt
        x += step
        m[ind, :] = x
    
    # Map the label for each compound to an array of its concentrations at the input times
    d = {compound.label: col for compound, col in zip(compounds, m.T, strict=True)}

    return d


param: TypeAlias = float|int|sympy.Basic
solution_type: TypeAlias = sympy.Expr|Callable[[float], float]


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

    def solve_analytic(self, lambdify: bool=True) -> dict[str, solution_type]:
        """Tries to solve the system analytically.
        Returns a dict mapping each compound label to its solution.
        If lambdify is True, each solution is a callable, representing the concentration as a function of time.
        Otherwise, each solution is a sympy expression representing the solution."""
        
        eqs = self.get_equations()
        ics = self.get_initial_conditions()
        solutions = dsolve(
            eqs,
            ics=ics
        )
        
        res = dict()
        for solution in solutions:
            expr = solution.rhs
            compound = self._A_t_to_compound[solution.lhs]
            this_res = sympy.lambdify(symbols.t, expr) if lambdify else expr
            res[compound.label] = this_res
            
        return res

    def solve_numerical(self, t_vals: np.ndarray):
        res = solve_numerical(compounds=self.compounds.values(), reactions=self.reactions, t_vals=t_vals)
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
    
    t_vals = np.linspace(0.0, 24.0, num=10000)
    num = model.solve_numerical(t_vals)
    
    for k, v in num.items():
        print(k, v[:5])
    