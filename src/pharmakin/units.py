import abc
from functools import wraps
import pint
from typing import final


ureg = pint.UnitRegistry()
Q_ = ureg.Quantity


class Dim:
    """Holds the default units for various dimensions (time, space, volume, etc)"""
    
    MASS = ureg.g
    VOLUME = ureg.ml
    TIME = ureg.hours
    
    @classmethod
    def to_unit_list(cls):
        """Returns a list of the units to be used as default units."""

        res = []
        for k in dir(cls):
            v = getattr(cls, k)
            if not k.startswith("_") and not callable(v):
                res.append(v)
            #
        
        return res
    #


ureg.default_preferred_units = Dim.to_unit_list()


def has_unit(value):
    return isinstance(value, pint.Quantity) and not value.dimensionless


def coerce_float(value, target_unit: pint.Unit=None):
    """Removes units from a value.
    If value doesn't have units, it is returned as-is.
    If value has units and target_unit is set, the value is converted to the target unit
    before dropping units. This is so the function can be used to take e.g. an arbitrary
    measure of volume (say, 4 dL) and return a dimensionless quantity representing it in another unit.
    For instance
        drop_units(value=4dL, target_unit=ml) will return 400."""

    if isinstance(value, pint.Quantity):
        converted = value.to(target_unit) if target_unit is not None else value
        res = converted.magnitude
    else:
        res = value
    return res


def coerce_unit(value, unit: pint.Unit):
    """Converts input value to the specified units.
    If the value already has units, it is converted."""

    if isinstance(value, pint.Quantity):
        res = value.to(unit)
    else:
        res = Q_(value, unit)

    return res


if __name__ == '__main__':
    pass
