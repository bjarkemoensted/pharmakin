import numpy as np
from unittest import TestCase

from pharmakin.modeling import example_models
from pharmakin.modeling.base import Model, FirstOrderModel


class TestInstantaneousRelease(TestCase):
    def setUp(self):
        self.label = "drug"
        self.t0 = 10.0
        self.amount = 60.0
        self.t = np.linspace(0.0, 100.0, num=100)
        # Time point when the administration kicks in
        self.cut = next(i for i, t in enumerate(self.t) if t >= self.t0)
        self.model = FirstOrderModel().add_compound(label=self.label)
        self.model.add_administration(
            label=self.label,
            kind="instantaneous",
            amount=self.amount,
            t0=self.t0
        )
    
    def _check_solver(self, how: str):
        solution = self.model.solve(how=how, t_vals=self.t)
        series = list(solution.values())[0]

        np.testing.assert_almost_equal(0.0, series[:self.cut])
        np.testing.assert_almost_equal(series[self.cut:], self.amount)

    def test_release_analytical(self):
        self._check_solver(how="analytic")
    #
    
    # def test_release_num(self):
    #     self._check_solver(how="numeric")
    # #
    
    

if __name__ == '__main__':
    t = TestInstantaneousRelease()
    t.setUp()
    t.test_release_num()