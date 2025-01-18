import sys
_this_module = sys.modules[__name__]

from pharmakin.utils.parameter_base import Parameter, ParameterMeta
from pharmakin.utils.units import Dim
from pharmakin.utils.utils import BulkImporter


get_all_parameters = BulkImporter(from_=_this_module, child_of=Parameter)


class dose(Parameter):
    """Dose parameter - the total dose (mass) administered"""
    unit = Dim.MASS


class concentration(Parameter):
    """Concentration of drug in plasma"""
    
    unit = Dim.MASS / Dim.VOLUME


class clearance(Parameter):
    """Clearence parameter. Indicates the volume of plasma cleared per time unit."""

    unit = Dim.VOLUME / Dim.TIME


class auc(Parameter):
    """Represents the area under the curve (AUC) of a concentration/time curve.
    The AUC is the integral over the plasma concentration over time after a dose is administered.
    AUC is a measure of the total exposure (relative to plasma volume) of a dose."""
    
    unit = (Dim.MASS / Dim.VOLUME) * Dim.TIME
    

class volume_of_distribution(Parameter):
    """The 'apparent volume of distribution', given by the total amount of drug in the body divided
    by the plasma concentration.
    Hence, this is a theoretical 'effective' volume of plasma. If a person has 5L of blood,
    but only 50% of a drug is dissolved in the blood, the v_d is 10L."""
    
    unit = Dim.VOLUME


class half_life(Parameter):
    """The time taken in first-order processes for an amount of drug to reduce by 50%."""

    unit = Dim.TIME


if __name__ == '__main__':
    vd = volume_of_distribution.example_values(size=None, with_units=True)
    print(vd)
