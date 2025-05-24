from functools import partial
import logging
logger = logging.getLogger(__name__)
import sympy
from sympy import Function, dsolve, Derivative, Eq
from sympy.core.function import AppliedUndef

from pharmakin.modeling import symbols


class Meh:
    target_symbol = "A"
    time_symbol = "t"
    slope_rel_threshold = 0.01
    
    def __init__(self, k: float):
        self.k = k
        
        self.t = sympy.Symbol(self.time_symbol)
        
        self.y = Function(self.target_symbol)
        self.yp = -self.k*self.y(self.t)
    
    @property
    def eq(self):
        return sympy.Eq(sympy.Derivative(self.y(self.t), self.t), self.yp)
    
    def solve(self, a_initial):
        solution = dsolve(self.eq, self.y(self.t), ics={self.y(0): a_initial})
        expr = solution.rhs
        f = sympy.lambdify((self.t), expr)
        return f
    
    def numerical(self, a_initial: float, t: float, n_steps: int):
        res = []
        dt = t/n_steps
        running = a_initial
        
        issued_warning = False
        
        grad = sympy.lambdify((self.y(self.t)), self.yp)
        for _ in range(n_steps):
            res.append(running)
            slope = grad(running)
            delta = slope*dt
            if abs(delta) > running*self.slope_rel_threshold and not issued_warning:
                logger.warning(f"{self} quantity changed more than allowed fraction ({self.slope_rel_threshold})!")
                issued_warning = True
            running += delta
        
        return res
    #


hmm = Meh(k=.01)
#print(hmm.yp.subs())


class Rate:
    """Represents metabolism or elimination - the temporal evolution of some quantity"""
    
    def __init__(self, expr: sympy.Expr):
        """expr is a sympy expression for the rate, i.e. dA/dt, where A is the quantity of interest"""
        
        self.expr = expr
    
    def get_vars(self) -> set:
        """Returns a set of function variables (like A_SUF(t)).
        Using a custom method for this because expr.free_symbols() would return t as the free symbol,
        rather than A(t)."""
        
        vars_ = self.expr.atoms(AppliedUndef)
        return vars_
    
    def substitute_quantities(self, quantities: dict):
        """Takes a dict mapping quantities (e.g. A_SUF(t)) to their current values.
        Returns the decay expression after substituting the values."""
        
        in_expression = self.get_vars()
        replace = {k: v for k, v in quantities.items() if k in in_expression}
        print(in_expression, replace)
        res = self.expr.subs(replace)
        return res
    
    def first_order_step(self, current_values: dict, dt: float) -> float:
        slope = self.substitute_quantities(current_values)
        delta = slope*dt
        return delta


class Quantity:
    def __init__(self, label: str, A_0: float=0.0):
        self.label = label
        self.A_0 = A_0
        self.decays: list[tuple[str|None, sympy.Expr]] = []
        
        A_suffix_str = f"{str(symbols.A)}_{self.label}"
        self.A = sympy.Function(A_suffix_str)

    def replace_with_suffixed_var(self, expr: sympy.Expr) -> sympy.Expr:
        """Takes a sympy expression with a generic symbol for drug amount (A).
        Replaces A with the variable representing this quantity.
        Note that A must be expressed as a function of t (symbols.A(symbols.t) can be used,
        for instance)."""
        
        old = symbols.A(symbols.t)
        new = self.A(symbols.t)
        res = expr.subs(old, new)
        return res

    def add_decay(self, decay: sympy.Expr, metabolite: str|None=None):
        rate = self.replace_with_suffixed_var(decay)
        conversion = (metabolite, rate)
        self.decays.append(conversion)


if __name__ == '__main__':
    q = Quantity("LDX")
    
    t = sympy.Symbol("t")
    A = sympy.Function("A")

    decay_expr_base = -0.01*A(t)
    decay_expr = q.replace_with_suffixed_var(decay_expr_base)
    
    #total_out = sum(tuple(zip(*q.decays))[1])
    
    rate = Rate(decay_expr)
    print(rate)
    
    d = {q.A(t): 10}
    print(d)
    
    grad = rate.substitute_quantities(d)
    print(grad, type(grad), float(grad))
    print(rate.first_order_step(d, dt=0.12))