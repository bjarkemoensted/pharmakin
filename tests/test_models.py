""" Contains some tests for models, checking if analytical and numerical solutions obtained match expected solutions.
Analyitcal first=order kinetics solutions for a single drug, and for a prodrug model are taken from:
Cho S, Yoon YR. Understanding the pharmacokinetics of prodrug and metabolite. Translational and Clinical Pharmacology. 2018"""

import numpy as np
import sympy
from sympy.core.function import AppliedUndef, Application, Derivative, UndefinedFunction
from typing import Iterable, Iterator
from unittest import TestCase

from pharmakin.modeling import example_models
from pharmakin.modeling.base import FirstOrderModel, valid_solve_methods, _solvers, solve_methods
from pharmakin.modeling.solvers import result_type, Solver


def test_all_solvers_have_a_keyword():
    """Check that all solvers are associated with a keyword (the 'how' parameter) so models can access them
    for solving ODEs."""
    
    assert set(valid_solve_methods) == set(_solvers.keys())
    assert all(issubclass(cls_, Solver) for cls_ in _solvers.values())


class Base(TestCase):
    # Which solvers to check results for. Default to all.
    check_solvers: Iterable[solve_methods] = valid_solve_methods
    
    t_low = 0.0
    t_high = 100.0
    t_n_values = 100_000
    
    @classmethod
    def get_tvals(cls) -> np.ndarray:
        res = np.linspace(cls.t_low, cls.t_high, num=cls.t_n_values)
        return res

    def _compare_numeric(self, *solutions: result_type, decimal: int|None=None, **kwargs):
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
    
    def iterate_solutions(self) -> Iterator[tuple[solve_methods, result_type]]:
        for method in self.check_solvers:
            t = self.get_tvals()
            solution = self.model.solve(t_vals=t, how=method)
            yield method, solution
        #
    #


class TestConstantModel(Base):
    """Case for a a single drug with zero amount and no reactions, just to make sure solvers behave
    reasonably when nothing is supposed to happen."""
    
    t_n_values = 100
    
    def setUp(self):
        self.model = FirstOrderModel()
        label = "some_drug"
        initial_amount = 0.0
        self.model.add_compound(label=label, initial_amount=initial_amount)
        self.correct = {label: np.array([initial_amount for _ in range(self.t_n_values)])}
        super().setUp()
    
    def test_vals(self):
        for _, solution in self.iterate_solutions():
            self._compare_numeric(solution, self.correct)
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
        t = self.get_tvals()
        self.correct_solution = {k: f(t) for k, f in self.solution.items()}

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
    
    def check_solver(self, how=solve_methods, decimal: int|None=None, **kwargs):
        t = self.get_tvals()
        solution = self.model.solve(t_vals=t, how=how)
        self._compare_numeric(solution, self.correct_solution, decimal=decimal, **kwargs)
    
    def test_analytic_solve(self):
        """Solves the single-drug model analytically, and compares with the expected result"""
        
        self.check_solver("analytic")
    
    def test_simple_numeric_solve(self):
        self.check_solver("simple", decimal=3)
        
    def test_numeric_solve(self):
        self.check_solver("numeric", decimal=3)


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


if __name__ == '__main__':
    t = TestConstantModel()
    t.setUp()
    t.test_vals()
