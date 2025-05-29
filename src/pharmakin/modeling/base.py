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


def solve_numerical(
        compounds: Iterable[Compound],
        reactions: Iterable[Reaction],
        t_vals: np.ndarray,
        max_relative_step_size: float=0.01):
    """Solves numerically a system consisting of the specified compounds and reactions.
    This takes a simple approach and computes at each time step:
        * The current reaction rates from the current drug concentrations.
        * How quickly the reaction rates are changing, given the current rates.
    The above steps are equivalent to computing a velocity given a current position, and an acceleration given
    the current velocities. For this reason, standard symbols representing position (x), velocity (v), and
    acceleration (a) are used in the code.
    compounds: Iterable of compunds in the model.
    reactions: Iterable of reactions in the model.
    t_vals: numpy array of time values. This must begin at 0.0, as the solver assumes the initial values
        stored under .A_0 in each compund are initially valid. The reason for explicitly specifying an
        array of time values (rather than e.g. a stop-time, and step size), is to make it easier to
        compare with analytical solutions, where concentrations are already available as a function of time,
        and must be passed a time array to obtain a sequence of concentrations.
    max_relative_step_size (float) - The maximum allowed change relative to current concentration before
        a warning is raised. For example, if this value is 0.01, and if at any time step sum of absolute changes
        in concentrations exceed 1% of the current concentration, the logger will issue a warning.
        The logger will only warn the first time this threshold is exceeded."""
    
    # Check we start at t=0 (because we initialize with the initial values A_0 for each compound)
    assert t_vals[0] == 0.0
    warned = False
    
    # Make sure compounds are ordered
    compounds = sorted(compounds, key=lambda c: c.label)
    
    # Determine the rate of change in concentration of each compound (dA/dt), as a function of current concentrations
    x_vars = tuple(c.A_t for c in compounds)
    gradients = _summarize_reactions(reactions)
    slopes = [gradients[c.Ap] for c in compounds]
    fp = sympy.lambdify(x_vars, slopes, modules="numpy")
    
    v_vars = tuple(c.Ap for c in compounds)
    second_derivatives = [sympy.Derivative(s, symbols.t).doit() for s in slopes]
    fpp = sympy.lambdify(v_vars, second_derivatives, modules="numpy")
    
    # Make a matrix for holding the results, and set the first row to initial values
    x = np.array([c.A_0 for c in compounds])
    m = np.empty(shape=(len(t_vals), len(compounds)))
    m.fill(np.nan)
    ind = 0
    m[ind, :] = x
    
    # Fill up the matrix using the gradients to interpolate from previous points
    for dt in np.diff(t_vals):
        ind += 1
        v = np.array(fp(*x))
        a = np.array(fpp(*v))
        step = dt*v +0.5*a*dt**2
        
        total_change = np.sum(np.abs(step))
        total_conc = np.sum(np.abs(x))
        toofast = total_change >= total_conc*max_relative_step_size
        if toofast and not warned:
            logger.warning(f"Step size exceeded threshold {max_relative_step_size}. Try greater temporal resolution.")
            warned = True
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

    t_vals = np.linspace(0.0, 1.0, num=1000)
    num = model.solve_numerical(t_vals)
