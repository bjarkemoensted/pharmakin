from __future__ import annotations
from abc import ABC, abstractmethod
import logging
logger = logging.getLogger(__name__)
import numpy as np
from scipy.integrate import solve_ivp
import numpy as np
import sympy
from typing import Iterable, TypeAlias

from pharmakin.modeling import symbols
from pharmakin.modeling.compound import Compound, get_initial_conditions
from pharmakin.modeling.reactions import Reaction, summarize_reactions


result_type: TypeAlias = dict[sympy.core.function.AppliedUndef, np.typing.NDArray[np.float64]]


class Solver(ABC):
    def __init__(self, compounds: Iterable[Compound], reactions: Iterable[Reaction]):
        # Convert functions to a list to ensure consistent ordering
        self.compounds = list(compounds)
        self.reactions = reactions
        
        self.funcs = [c.A_t for c in self.compounds]
        self.ics = get_initial_conditions(*self.compounds)
        
        # Symbols and expressions for time derivatives
        self.gradients = summarize_reactions(*self.reactions)
        self.dadt = [sympy.Derivative(a, symbols.t) for a in self.funcs]
        self.slopes = [self.gradients.get(Ap, 0.0) for Ap in self.dadt]
    #
    
    @abstractmethod
    def solve(self, t: np.ndarray) -> result_type:
        raise NotImplementedError
    #


class AnalyticalSolver(Solver):
        
    @property
    def eqs(self) -> list[sympy.Eq]:
        res = [sympy.Eq(Ap, slope) for Ap, slope in zip(self.dadt, self.slopes, strict=True)]
        return res

    def _get_analytic_symbolic_solutions(self) -> dict[sympy.core.function.AppliedUndef, sympy.Expr]:
        """Map each variable (f(t)) to a symbolic expression for its solution."""
        
        solutions = sympy.dsolve(self.eqs, ics=self.ics)
        res = {s.lhs: s.rhs for s in solutions}
        return res
    
    def solve(self, t: np.ndarray) -> result_type:
        """Solves the problem analytically.
        The exact solutions are then converted into numeric functions and evaluated, to
        obtain a result consistent with numerical methods."""
        
        res = dict()
        solutions = self._get_analytic_symbolic_solutions()
        for variable, expr in solutions.items():
            f = sympy.lambdify(symbols.t, expr)
            vals = f(t)
            
            # sympy lambdified functions return scalars for contant functions, so convert those to arrays
            if not isinstance(vals, np.ndarray):
                vals = np.array([float(vals) for _ in range(len(t))])
            
            res[variable] = vals
        
        return res
    #


class SimpleSolver(Solver):
    def solve(self, t, max_relative_step_size: float=0.01):
        """Solves numerically a system of differential equations.
        This takes a simple approach and computes at each time step:
            * The current rate of change given the current values.
            * How quickly said rates are changing, given the current rates.
        The above steps are equivalent to computing a velocity given a current position, and an acceleration given
        the current velocities. For this reason, standard symbols representing position (x), velocity (v), and
        acceleration (a) are used in the code.
        
        max_relative_step_size (float) - The maximum allowed change relative to current values before
            a warning is raised. For example, if this value is 0.01, and if at any time step sum of absolute changes
            in values exceed 1% of the current concentration, the logger will issue a warning.
            The logger will only warn the first time this threshold is exceeded."""

        warned = False

        x_vars = tuple(self.funcs)

        # Determine the rate of change in concentration of each compound (dA/dt), as a function of current concentrations
        fp = sympy.lambdify(x_vars, self.slopes, modules="numpy")
        
        # Express acceleration of change (d^2A/dt^2) as a function of current rate of change
        v_vars = tuple(self.dadt)
        second_derivatives = [sympy.Derivative(s, symbols.t).doit() for s in self.slopes]
        fpp = sympy.lambdify(v_vars, second_derivatives, modules="numpy")

        # Make a matrix for holding the results, and set the first row to initial values
        x = np.array([self.ics[A_t.subs(symbols.t, 0)] for A_t in self.funcs])
        m = np.empty(shape=(len(t), len(self.funcs)))
        m.fill(np.nan)
        ind = 0
        m[ind, :] = x

        # Fill up the matrix using the gradients to interpolate from previous points
        for dt in np.diff(t):
            ind += 1
            v = np.array(fp(*x))
            a = np.array(fpp(*v))
            step = dt*v +0.5*a*dt**2
            
            # Throw a warning if the rate of change in concentration is relatively large (compared with current value)
            total_change = np.sum(np.abs(step))
            total_conc = np.sum(np.abs(x))
            # Add a tiny threshold to avoid warnings when concentrations are ~zero
            threshold = 10e-8
            toofast = total_change >= total_conc*max_relative_step_size + threshold
            if toofast and not warned:
                logger.warning(f"Step size exceeded threshold {max_relative_step_size}. Try greater temporal resolution.")
                warned = True
            
            x += step
            m[ind, :] = x

        # Map the label for each compound to an array of its concentrations at the input times
        d = {A_t: col for A_t, col in zip(self.funcs, m.T, strict=True)}

        return d
    #


class NumericSolver(Solver):
    def solve(self, t, method="DOP853", rtol=0.00001, **kwargs):
        """Uses numerical integration for solving the system."""
        
        time_derivatives = sympy.lambdify((symbols.t, self.funcs), self.slopes, modules="numpy")

        def system(t, z):
            return time_derivatives(t, z)
        
        z0 = [self.ics[A_t.subs(symbols.t, 0)] for A_t in self.funcs]
        t_span = (min(t), max(t))

        kw = dict(method=method, rtol=rtol, **kwargs)
        sol = solve_ivp(system, t_span, z0, t_eval=t, **kw)
        
        res = dict(zip(self.funcs, sol.y, strict=True))
        return res
    #
