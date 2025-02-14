import pint
from types import UnionType
from typing import final, get_args

from pharmakin.utils.formula import Formula
from pharmakin.utils.units import ureg, Dim, Q_, coerce_float, coerce_unit
from pharmakin.utils.stats import Simulator


# Define required attributes and their types which must be set in all parameter classes
REQUIRED_ATTRIBUTES = (
    ("unit", None|pint.Unit),
)


def _type_is_correct(value, type_condition):
    """Checks whether a value matches a certain type.
    The type_condition can be either a class, None, or a UnionType, as defined
    in a standard type hint e.g. float|int."""

    if type_condition is None:
        return value is None
    
    allowed = type_condition
    # If using a type hint-like input, extract the allowed types from it
    if isinstance(type_condition, UnionType):
        allowed = get_args(type_condition)
    
    return isinstance(value, allowed)


class ParameterMeta(type):
    def __new__(cls, name, bases, dct):
        """Called when a subclass is defined. Checks that a unit is correctly declared and registers the class."""
        # Enforce that unit is specified explicitly in the class body (i.e. not missing or inherited)
        
        for attr, type_hint in REQUIRED_ATTRIBUTES:
            # Check that required attributes are explicitly declared in the class body
            if attr not in dct:
                raise RuntimeError(f"No unit set for parameter: {cls}.")
            
            # Check that attributes have required types
            attr_val = dct[attr]
            if not _type_is_correct(value=attr_val, type_condition=type_hint):
                raise TypeError(f"Class attribute {attr} ({attr_val}) must be of type {type_hint}.")
        
        res = super().__new__(cls, name, bases, dct)
        return res
    
    def __repr__(self):
        return f"{self.__name__} <parameter>"
    
    def __str__(self):
        return repr(self)


class Parameter(metaclass=ParameterMeta):
    """Base class for representing pharmacokinetic parameters.
    Contains various class attributes and methods for validating parameters, converting units, etc.
    Specifically, it declares a unit for the parameter, and has functionality for validating values for the parameter.
    Units must be specified explicitly when subclassing."""

    unit = None
    lower = 0.0
    upper = float("inf")
    
    @classmethod
    def _unit_is_valid(cls, value: float|pint.Quantity|pint.Unit) -> bool:
        """Checks if value has units compatible with the parameter"""
        
        if isinstance(value, float):
            return True
    
        unit = value if isinstance(value, pint.Unit) else value.units
        return unit.is_compatible_with(cls.unit)
    
    @classmethod
    def ensure_units(cls, value: float|pint.Quantity, force_standard_units=False) -> pint.Quantity:
        """Ensures that the input value has units compatible with the parameter.
        value: The value to be converted. Can be float or a pint.Quantity.
            If value is a float, it is converted to the parameter's default units.
            Otherwise, a dimension check is made and the value is optionally converted into default units.
        force_standard_units: Whether to convert an input with compatible units to the default units (e.g. L -> mL)"""

        # If input has units and we're not force-converting, just check for compatibility
        if isinstance(value, pint.Quantity) and not force_standard_units:
            if not cls._unit_is_valid(value):
                raise RuntimeError(f"Input ({value}) has units incompatible with {cls.unit}.")
            return value
        
        # Otherwise coerce to default units
        return coerce_unit(value=value, unit=cls.unit)
    
    @classmethod
    def ensure_float(cls, value: float|pint.Quantity) -> float:
        """Convert a parameter value into a float representing its value in the default units."""
        res = coerce_float(value=value, target_unit=cls.unit)
        return res

    @classmethod
    def is_valid(cls, value: float|pint.Quantity):
        if not cls._unit_is_valid(value):
            return False
    
        value = cls.ensure_float(value)
        res = cls.lower <= value <= cls.upper
        return res

    @final
    @classmethod
    def validate(cls, value: float|pint.Quantity):
        if not cls.is_valid(value):
            raise ValueError(f"Invalid value: {value}.")
        #
    
    @classmethod
    def formula(cls, func):
        """Implements a decorator for the class which adds unit handling and validation functionality
        to a wrapped function."""
        
        wrapped = Formula(func=func, result_class=cls)
        return wrapped
    #
    
    @classmethod
    def example_values(cls, size=None, with_units=False, seed: int=None):
        """Generates one or more examples of values the parameter could take.
        size is the number of values desired (None produces a single value).
        If with_units, the values will include values.
        seed is an optional random seed."""

        # Some parameters have no upper bound - just use some fixed value for those
        upper = cls.upper
        if upper == float("inf"):
            upper = 10
        
        # Simulate value(s)
        f = Simulator(distribution="uniform", seed=seed, low=cls.lower, high=upper)
        res = f(size=size)
        
        # Add units if required
        if with_units:
            if isinstance(res, float):
                res = cls.ensure_units(res)
            else:
                res = [cls.ensure_units(val) for val in res]
            #

        return res            


if __name__ == '__main__':
    pass
    