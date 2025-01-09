from collections import defaultdict
import inspect
import sympy
from typing import Callable

from pharmakin.quant import formulas


class Registry(dict):
    """Registry class for uniquely registering quantities under a given name.
    Similar to a dict, but throws an error if a name has already been registered."""

    def __init__(self):
        self.formulary = formulas.Formulary()

    def __setitem__(self, key, value):
        if key in self:
            raise RuntimeError(f"Key {key} is already set (to {value})")
        return super().__setitem__(key, value)
    #
    
    def formula(self, class_method: Callable):
        print("WEEEE registering!!!", class_method)
        bound = inspect.signature(class_method)
        print(bound)
        return class_method


if __name__ == '__main__':
    pass
