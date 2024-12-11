import pint
from typing import final

from pharmakin.registry import Registry
from pharmakin.units import ureg, Dim, Q_, coerce_float, coerce_unit


PARAMETER_REGISTRY = Registry()


class Parameter():
    """Base class for representing pharmacokinetic parameters.
    Contains various class attributes and methods for validating parameters, converting units, etc.
    Specifically, it declares a unit for the parameter, and has functionality for validating values for the parameter.
    Units must be specified explicitly when subclassing.
    Subclasses are automatically registered to allow lookup from their class names.
    Note that this classes cannot be instantiated - it solely exists to offer the afforementioned functionality, and
    to employ some of python's inheritance machinery to simplify autoregistration of parameters.
    
    To override the autoregistration (if making a Parameter subclass which should not be registered), do
    class MyClass(Parameter, register=False)"""

    unit = None
    lower = 0.0
    upper = float("inf")
    
    @final
    def __init__(self):
        raise RuntimeError("Parameter classes are solely for namespacing and registering. Don't instantiate.")
    
    def __init_subclass__(cls, register=True):
        """Called when a subclass is defined. Checks that a unit is correctly declared and registers the class."""
        # Enforce that unit is specified explicitly in the class body (i.e. not missing or inherited)
        attr = "unit"
        if attr not in cls.__dict__:
            raise RuntimeError(f"No unit set for parameter: {cls}.")
        
        # Enforce that unit must by a pint.Unit
        unit = getattr(cls, attr)
        if not isinstance(unit, pint.Unit):
            raise TypeError
        
        # Register the parameter under its name
        if register:
            PARAMETER_REGISTRY[cls.__name__] = cls
    #
    
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
        if not cls.is_valid(value):
            raise ValueError(f"Invalid value: {val} (validated {value}).")
        #
    #
   

class dose(Parameter):
    """Dose parameter - the total dose (mass) administered"""
    unit = Dim.MASS


class clearance(Parameter):
    """Clearence parameter. Indicates the volume of plasma cleared per time unit."""

    lower = 0.0
    upper = float("inf")
    unit = Dim.VOLUME / Dim.TIME
    
    @classmethod
    def from_auc_dose(cls, auc: float, dose: float):
        res = dose/auc
        return res
    #


class auc(Parameter):
    """Represents the area under the curve (AUC) of a concentration/time curve.
    The AUC is the integral over the plasma concentration over time after a dose is administered.
    AUC is a measure of the total exposure (relative to plasma volume) of a dose."""
    
    unit = (Dim.MASS / Dim.VOLUME) * Dim.TIME


if __name__ == '__main__':
    pass
    