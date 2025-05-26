""" Contains some tests for models, checking if analytical and numerical solutions obtained match expected solutions.
Analyitcal first=order kinetics solutions for a single drug, and for a prodrug model are taken from:
Cho S, Yoon YR. Understanding the pharmacokinetics of prodrug and metabolite. Translational and Clinical Pharmacology. 2018"""

import numpy as np
import sympy
from unittest import TestCase

from pharmakin.kinetics.first_order import k_el
from pharmakin.modeling.base import FirstOrderModel
from pharmakin.modeling import symbols


def single_drug_example(
        t_half: float|int,
        initial_amount: float=100.0,
        label: str="drug") -> tuple[FirstOrderModel, list[sympy.Expr]]:
    """Sets up a model for a single drug being metabolized (with first=order kinetics).
    Returns the model and the analyitical solution. For consistency with the prodrug model,
    the solution is provided as a list of a single sympy expression, representing the
    drug concentration as a function of time."""
    
    k = k_el(half_life=t_half)
    # Make the analyitical solution
    solutions = [initial_amount*sympy.exp(-k*symbols.t)]
    
    # Set up the model
    model = FirstOrderModel().add_compound(
        label, initial_amount=initial_amount
    ).add_reaction(label, k=k)
    
    return model, solutions


def prodrug_example(
        t_half_prodrug: float|int,
        t_half_active: float|int,
        A_0: float=100.0,
        label_prodrug: str="prodrug",
        label_active: str="active"
    ) -> tuple[FirstOrderModel, list[sympy.Expr]]:
    """Makes model and solution for a prodrug system, with a prodrug metabolized into an active drug."""
    
    k1 = k_el(half_life=t_half_prodrug)
    k2 = k_el(half_life=t_half_active)
    
    # Construct analytic solutions for prodrug and active drug, respectively
    s1 = A_0*sympy.exp(-k1*symbols.t)
    s2 = None
    if k1 != k2:
        norm = (k1*A_0)/(k2 - k1)
        s2 = norm*(sympy.exp(-k1*symbols.t) - sympy.exp(-k2*symbols.t))
    else:
        # Special case when the rates are identical (think this can be obtained via L'Hopital from the general case)
        s2 = k1*A_0*symbols.t*sympy.exp(-k1*symbols.t)
        
    solutions = [s1, s2]
    
    # Set up the model
    model = FirstOrderModel().add_compound(
        label_prodrug, initial_amount=A_0
    ).add_compound(
        label_active, initial_amount=0.0
    ).add_reaction(
        label_prodrug, label_active, k=k1
    ).add_reaction(
        label_active, k=k2
    )
    
    return model, solutions


class TestFirstOrderSolve(TestCase):
    
    def setUp(self):
        t_half_amp = 1.0
        t_half_ldx = 10.5
        # Example - metabolization of dextroamphetamine
        self.amp_model, self.amp_solutions = single_drug_example(label="AMP", t_half=t_half_amp, initial_amount=20.0)
        
        # Example - lisdexamphetamine prodrug
        self.ldx_model, self.ldx_solutions = prodrug_example(
            label_prodrug="LDX", label_active="AMP",
            t_half_active=t_half_amp, t_half_prodrug=t_half_ldx,
            A_0=60.0)
    
    def _check_solutions_equiv(self, s1: sympy.Expr, s2: sympy.Expr, t_vals: np.ndarray|None=None):
        """Checks if the 2 provided sympy expressions are (approximately) the same.
        Converts both into functions and checks that a number of values for t result in very close results."""
        if not t_vals:
            t_vals = np.linspace(0.0, 100.0, num=1000)
        
        f1, f2 = tuple(sympy.lambdify(symbols.t, e) for e in (s1, s2))
        np.testing.assert_almost_equal(f1(t_vals), f2(t_vals))
    
    def test_analytic_single(self):
        """Solves the single-drug model analytically, and compares with the expected result"""
        correct = self.amp_solutions[0]
        model_solutions = self.amp_model.solve_analytic()
        for model_sol, correct in zip(model_solutions, self.amp_solutions):
            self._check_solutions_equiv(model_sol.rhs, correct)
    
    def test_analytic_prodrug(self):
        """Solves the pro-drug model analytically and compares with expected result"""
        model_solutions = self.ldx_model.solve_analytic()
        for model_sol, correct in zip(model_solutions, self.ldx_solutions):
            self._check_solutions_equiv(model_sol.rhs, correct)
        #
    #
