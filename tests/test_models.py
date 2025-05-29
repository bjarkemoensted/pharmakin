""" Contains some tests for models, checking if analytical and numerical solutions obtained match expected solutions.
Analyitcal first=order kinetics solutions for a single drug, and for a prodrug model are taken from:
Cho S, Yoon YR. Understanding the pharmacokinetics of prodrug and metabolite. Translational and Clinical Pharmacology. 2018"""

import numpy as np
import sympy
from typing import Callable
from unittest import TestCase

from pharmakin.kinetics.first_order import k_el
from pharmakin.modeling.base import FirstOrderModel, solution_type
from pharmakin.modeling import symbols


def _lambdify_solutions(solutions: dict[str, sympy.Expr]) -> dict[str, Callable[[float], float]]:
    res = {k: sympy.lambdify(symbols.t, v) for k, v in solutions.items()}
    return res


def single_drug_example(
        t_half: float|int,
        initial_amount: float=100.0,
        label: str="drug") -> tuple[FirstOrderModel, dict[str, sympy.Expr]]:
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
    ) -> tuple[FirstOrderModel, dict[str, sympy.Expr]]:
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
        self.amp_model, amp_solutions = single_drug_example(label="AMP", t_half=t_half_amp, initial_amount=20.0)
        self.amp_solutions = _lambdify_solutions(amp_solutions)
        
        # Example - lisdexamphetamine prodrug
        self.ldx_model, ldx_solutions = prodrug_example(
            label_prodrug="LDX", label_active="AMP",
            t_half_active=t_half_amp, t_half_prodrug=t_half_ldx,
            A_0=60.0)
        self.ldx_solutions = _lambdify_solutions(ldx_solutions)
        
        self.models = (self.amp_model, self.ldx_model)

    def get_tvals(self) -> np.ndarray:
        res = np.linspace(0.0, 100.0, num=100_000)
        return res

    def test_model_data_types(self):
        """Checks that the model uses expected data types. Adding this test because sympy types are a bit tricky
        to keep track of, e.g. instantiating a 'Function' doesn't give a 'Function' instance but rather an
        UndefinedFunction etc."""

        for model in self.models:
            for compound in model.compounds.values():
                self.assertIsInstance(compound.A, sympy.core.function.UndefinedFunction)
                self.assertIsInstance(compound.A_t, sympy.core.function.AppliedUndef)
    
    def _compare_numeric(self, *solutions: dict[str, np.ndarray], decimal: int|None=None, **kwargs):
        """Checks if the 2 provided sympy expressions are (approximately) the same.
        Converts both into functions and checks that a number of values for t result in very close results."""
        
        assert len(solutions) > 1
        assert all(isinstance(arr, np.ndarray) for d in solutions for arr in d.values())
        
        if decimal is not None:
            kwargs.update(decimal=decimal)
        for s1, s2 in zip(solutions[:-1], solutions[1:]):
            self.assertEqual(set(s1.keys()), set(s2.keys()))
            for k, arr1 in s1.items():
                arr2 = s2[k]
                np.testing.assert_almost_equal(arr1, arr2, **kwargs)
            #
        #
    
    def _compare_funcs(self, *solutions: dict[str, Callable[[float], float]], **kwargs):
        assert all(callable(f) for d in solutions for f in d.values())
        
        t = self.get_tvals()
        vals = ({label: func(t) for label, func in s.items()} for s in solutions)
        return self._compare_numeric(*vals, **kwargs)
    
    def test_analytic_single(self):
        """Solves the single-drug model analytically, and compares with the expected result"""
        correct = self.amp_solutions
        model_solutions = self.amp_model.solve_analytic()
        self._compare_funcs(model_solutions, correct)
    
    def test_analytic_prodrug(self):
        """Solves the pro-drug model analytically and compares with expected result"""
        model_solutions = self.ldx_model.solve_analytic()
        self._compare_funcs(model_solutions, self.ldx_solutions)
    
    def test_numeric_single(self):
        t = self.get_tvals()
        correct = {label: f(t) for label, f in self.amp_solutions.items()}
        
        num = self.amp_model.solve_numerical(t_vals=t)
        self._compare_numeric(num, correct, decimal=3)
    
    def test_numeric_prodrug(self):
        t = self.get_tvals()
        correct = {label: f(t) for label, f in self.ldx_solutions.items()}
        
        num = self.ldx_model.solve_numerical(t_vals=t)
        self._compare_numeric(num, correct, decimal=3)
    #


if __name__ == '__main__':
    t = TestFirstOrderSolve()
    t.setUp()
    t.test_numeric_single()
    t.test_numeric_prodrug()
