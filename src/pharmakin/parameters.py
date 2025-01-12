import sys
_this_module = sys.modules[__name__]

from pharmakin.utils.parameter_base import Parameter, ParameterMeta
from pharmakin.utils.units import Dim
from pharmakin.utils.utils import BulkImporter


get_all_parameters = BulkImporter(from_=_this_module, child_of=Parameter)


class dose(Parameter):
    """Dose parameter - the total dose (mass) administered"""
    unit = Dim.MASS


class clearance(Parameter):
    """Clearence parameter. Indicates the volume of plasma cleared per time unit."""

    lower = 0.0
    upper = float("inf")
    unit = Dim.VOLUME / Dim.TIME


class auc(Parameter):
    """Represents the area under the curve (AUC) of a concentration/time curve.
    The AUC is the integral over the plasma concentration over time after a dose is administered.
    AUC is a measure of the total exposure (relative to plasma volume) of a dose."""
    
    unit = (Dim.MASS / Dim.VOLUME) * Dim.TIME
    

if __name__ == '__main__':
    pass
