import numpy as np
from unittest import TestCase

from pharmakin.modeling import example_models
from pharmakin.modeling.base import Model, FirstOrderModel


class TestInstantaneousRelease(TestCase):
    def setUp(self):
        self.t0 = 10.0
        self.amount = 60.0
        self.t = np.linspace(0.0, 100.0, num=100)
        self.model = FirstOrderModel().add_compound(label='drug')
    
    def test_release(self):
        solution = self.model.solve(how="simple", t_vals=self.t)
        series = list(solution.values())[0]
        cut = next(i for i, t in enumerate(self.t) if t >= self.t0)
        np.testing.assert_almost_equal(0.0, series[:cut])
        self.assertAlmostEqual(series[cut], self.amount, places=1)
        
        
        
        

if __name__ == '__main__':
    t = TestInstantaneousRelease()
    t.setUp()
    t.test_release()