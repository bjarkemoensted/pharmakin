import numpy as np
from unittest import TestCase

from pharmakin.modeling import example_models
from pharmakin.modeling.base import Model, FirstOrderModel


class TestInstantaneousRelease(TestCase):
    def setUp(self):
        self.t0 = 10.0
        self.t = np.linspace(0.0, 100.0, num=100)
        self.model = FirstOrderModel().add_compound(label='drug')
    
    def test_stationary(self):
        solution = self.model.solve(how="simple", t_vals=self.t)
        print(solution)
        

if __name__ == '__main__':
    t = TestInstantaneousRelease()
    t.setUp()
    t.test_stationary()