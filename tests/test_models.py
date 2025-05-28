""" Contains some tests for models, checking if analytical and numerical solutions obtained match expected solutions.
Analyitcal first=order kinetics solutions for a single drug, and for a prodrug model are taken from:
Cho S, Yoon YR. Understanding the pharmacokinetics of prodrug and metabolite. Translational and Clinical Pharmacology. 2018"""

import numpy as np
import sympy
from typing import TypeAlias
from unittest import TestCase

from pharmakin.kinetics.first_order import k_el
from pharmakin.modeling.base import FirstOrderModel, solution_type
from pharmakin.modeling import symbols


def single_drug_example(
        t_half: float|int,
        initial_amount: float=100.0,
        label: str="drug") -> tuple[FirstOrderModel, solution_type]:
    """Sets up a model for a single drug being metabolized (with first=order kinetics).
    Returns the model and the analyitical solution. For consistency with the prodrug model,
    the solution is provided as a list of a single sympy expression, representing the
    drug concentration as a function of time."""
    
    k = k_el(half_life=t_half)
    
    # Set up the model
    model = FirstOrderModel().add_compound(
        label, initial_amount=initial_amount
    ).add_reaction(label, k=k)
    
    # Make the analyitical solution
    solutions = {label: initial_amount*sympy.exp(-k*symbols.t)}
    
    return model, solutions


def prodrug_example(
        t_half_prodrug: float|int,
        t_half_active: float|int,
        A_0: float=100.0,
        label_prodrug: str="prodrug",
        label_active: str="active"
    ) -> tuple[FirstOrderModel, solution_type]:
    """Makes model and solution for a prodrug system, with a prodrug metabolized into an active drug."""
    
    k1 = k_el(half_life=t_half_prodrug)
    k2 = k_el(half_life=t_half_active)
    
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
    
    # Construct analytic solutions for prodrug and active drug, respectively
    solutions = dict()
    solutions[label_prodrug] = A_0*sympy.exp(-k1*symbols.t)
    if k1 != k2:
        norm = (k1*A_0)/(k2 - k1)
        solutions[label_active] = norm*(sympy.exp(-k1*symbols.t) - sympy.exp(-k2*symbols.t))
    else:
        # Special case when the rates are identical (think this can be obtained via L'Hopital from the general case)
        solutions[label_active] = k1*A_0*symbols.t*sympy.exp(-k1*symbols.t)
        
    
    return model, solutions


class TestFirstOrderSolve(TestCase):
    
    def setUp(self):
        t_half_amp = 10.5
        t_half_ldx = 1.0
        # Example - metabolization of dextroamphetamine
        self.amp_model, self.amp_solutions = single_drug_example(label="AMP", t_half=t_half_amp, initial_amount=20.0)
        
        # Example - lisdexamphetamine prodrug
        self.ldx_model, self.ldx_solutions = prodrug_example(
            label_prodrug="LDX", label_active="AMP",
            t_half_active=t_half_amp, t_half_prodrug=t_half_ldx,
            A_0=60.0)
        
        self.models = (self.amp_model, self.ldx_model)

    def test_model_data_types(self):
        """Checks that the model uses expected data types. Adding this test because sympy types are a bit tricky
        to keep track of, e.g. instantiating a 'Function' doesn't give a 'Function' instance but rather an
        UndefinedFunction etc."""

        for model in self.models:
            for compound in model.compounds.values():
                self.assertIsInstance(compound.A, sympy.core.function.UndefinedFunction)
                self.assertIsInstance(compound.A_t, sympy.core.function.AppliedUndef)
    
    def _check_solutions_equiv(self, s1: solution_type, s2: solution_type, t_vals: np.ndarray|None=None):
        """Checks if the 2 provided sympy expressions are (approximately) the same.
        Converts both into functions and checks that a number of values for t result in very close results."""
        
        self.assertEqual(set(s1.keys()), set(s2.keys()))
        
        if not t_vals:
            t_vals = np.linspace(0.0, 100.0, num=1000)
        
        for compound_label in sorted(s1.keys()):
            expressions = (d[compound_label] for d in (s1, s2))
            funcs = (sympy.lambdify(symbols.t, e) for e in expressions)
            f1, f2 = funcs
            np.testing.assert_almost_equal(f1(t_vals), f2(t_vals))
    
    def test_analytic_single(self):
        """Solves the single-drug model analytically, and compares with the expected result"""
        correct = self.amp_solutions
        model_solutions = self.amp_model.solve_analytic()
        self._check_solutions_equiv(model_solutions, correct)
    
    def test_analytic_prodrug(self):
        """Solves the pro-drug model analytically and compares with expected result"""
        model_solutions = self.ldx_model.solve_analytic()
        self._check_solutions_equiv(model_solutions, self.ldx_solutions)
    
    def test_numeric_single(self):
        # TODO make less messy
        correct = self.amp_solutions["AMP"]
        n = 1_000_000
        T = 100.0
        t_vals = np.linspace(0.0, T, num=n)
        f = sympy.lambdify(symbols.t, correct)
        delta_t = t_vals[1] - t_vals[0]
        
        num = self.amp_model.solve_numerical(delta_t=delta_t, T=T)
        assert len(num) == 1
        vals_num = np.array([float(v) for v in list(num.values())[0]])
        print(vals_num[:3])
        vals_ana = np.array(f(t_vals))
        
        np.testing.assert_almost_equal(vals_num[:len(vals_ana)], vals_ana, decimal=4)

        pass
    #


if __name__ == '__main__':
    t = TestFirstOrderSolve()
    t.setUp()
    t.test_numeric_single()
