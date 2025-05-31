""" Contains some tests for models, checking if analytical and numerical solutions obtained match expected solutions.
Analyitcal first=order kinetics solutions for a single drug, and for a prodrug model are taken from:
Cho S, Yoon YR. Understanding the pharmacokinetics of prodrug and metabolite. Translational and Clinical Pharmacology. 2018"""

import numpy as np
import sympy
from sympy.core.function import AppliedUndef, Application, Derivative, UndefinedFunction
from typing import Callable
from unittest import TestCase

from pharmakin.kinetics.first_order import k_el
from pharmakin.modeling.base import FirstOrderModel, solution_type
from pharmakin.modeling import example_models
from pharmakin.modeling import symbols


AMP_T_HALF = 10.5
AMP_LABEL = "AMP"
LDX_T_HALF = 1.0
LDX_LABEL = "LDX"


class Base(TestCase):
    @staticmethod
    def get_tvals() -> np.ndarray:
        res = np.linspace(0.0, 100.0, num=100_000)
        return res

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
    #

class TestSingleDrugModel(Base):
    """Case for a single compound, first order kinetics."""
    
    @staticmethod
    def make_example() -> example_models.Example:
        return example_models.dexamphetamine_example()
    
    def setUp(self):
        self.example = self.make_example()
        self.model = self.example.model
        self.solution = self.example.solution

    def test_model_data_types(self):
        """Checks that the model uses expected data types. Adding this test because sympy types are a bit tricky
        to keep track of, e.g. instantiating a 'Function' doesn't give a 'Function' instance but rather an
        UndefinedFunction etc."""

        for compound in self.model.compounds.values():
            self.assertIsInstance(compound.A, UndefinedFunction)
            self.assertIsInstance(compound.A_t, AppliedUndef)
        
        for eq in self.model.get_equations():
            self.assertIsInstance(eq.lhs, Derivative)
            self.assertIsInstance(eq.rhs, sympy.Expr)

        for cond, val in self.model.get_initial_conditions().items():
            self.assertIsInstance(cond, Application)
            self.assertIsInstance(val, (int, float))
        #
    
    def test_analytic_solve(self):
        """Solves the single-drug model analytically, and compares with the expected result"""
        
        t = self.get_tvals()
        correct = {k: f(t) for k, f in self.solution.items()}
        model_solutions = self.model.solve(t_vals=t, how="analytic")
        self._compare_numeric(model_solutions, correct)
    
    def test_simple_numeric_solve(self):
        t = self.get_tvals()
        correct = {label: f(t) for label, f in self.solution.items()}
        
        model_solutions = self.model.solve(t_vals=t, how="simple")
        self._compare_numeric(model_solutions, correct, decimal=3)
    
    def test_numeric_solve(self):
        t = self.get_tvals()
        correct = {label: f(t) for label, f in self.solution.items()}
        
        model_solutions = self.model.solve(t_vals=t, how="numeric")
        self._compare_numeric(model_solutions, correct, decimal=3)
    #


class TestProDrugModel(TestSingleDrugModel):
    """Case for a prodrug model (first order kinetics)."""
    @staticmethod
    def make_example():
        return example_models.lisdexamphetamine_example()
    #


class ProdrugWithSingleRate(TestProDrugModel):
    """Case for a prodrug model (first order kinetics) in which the prodrug and active drug
    have the same half life. The analyitical solution looks different in this special case,
    so adding a separate test case to be sure it doesn't act weird."""
    
    @staticmethod
    def make_example():
        t_half = 2.0
        res = example_models.Example.linear_prodrug(
            t_half_prodrug=t_half,
            t_half_active=t_half
        )
        return res
    #
