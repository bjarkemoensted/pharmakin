import pint
from types import UnionType
from typing import final, get_args

from pharmakin.utils.formulas import Formula
from pharmakin.utils.units import ureg, Dim, Q_, coerce_float, coerce_unit



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
    def _validate_units(cls, value: float|pint.Quantity|pint.Unit):
        """Raises an error if value has units incompatible with the parameter"""
        if isinstance(value, pint.Quantity) and not value.units.is_compatible_with(cls.unit):
            raise RuntimeError(f"Input unit ({value.units}) not compatible with {cls.unit}.")
    
    @classmethod
    def ensure_units(cls, value: float|pint.Quantity, force_standard_units=False) -> pint.Quantity:
        """Ensures that the input value has units compatible with the parameter.
        value: The value to be converted. Can be float or a pint.Quantity.
            If value is a float, it is converted to the parameter's default units.
            Otherwise, a dimension check is made and the value is optionally converted into default units.
        force_standard_units: Whether to convert an input with compatible units to the default units (e.g. L -> mL)"""

        # If input has units and we're not force-converting, just check for compatibility
        if isinstance(value, pint.Quantity) and not force_standard_units:
            cls._validate_units(value=value)
            return value
        
        # Otherwise coerce to default units
        return coerce_unit(value=value, unit=cls.unit)
    
    @classmethod
    def ensure_float(cls, value: float|pint.Quantity) -> float:
        """Convert a parameter value into a float representing its value in the default units."""
        res = coerce_float(value=value, target_unit=cls.unit)
        return res

    @classmethod
    def is_valid(cls, value: float):
        res = cls.lower <= value <= cls.upper
        return res

    @final
    @classmethod
    def validate(cls, value: float|pint.Quantity):
        val = cls.ensure_float(value=value)
        if not cls.is_valid(val):
            raise ValueError(f"Invalid value: {val} (validated {value}).")
        #
    
    @classmethod
    def formula(cls, func):
        """Implements a decorator for the class which adds unit handling and validation functionality
        to a wrapped function."""
        
        wrapped = Formula(func=func, result_class=cls)
        return wrapped
    #


if __name__ == '__main__':
    pass
    