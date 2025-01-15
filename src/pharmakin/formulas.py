import sys
_this_module = sys.modules[__name__]

from pharmakin import parameters
from pharmakin.utils.formulas import Formula
from pharmakin.utils.utils import BulkImporter


get_all_formulas = BulkImporter(from_=_this_module, instance_of=Formula)


@parameters.clearance.formula
def clearence_from_dose_auc(dose, auc):
    """Determines clearence from dose and AUC"""
    res = dose / auc
    return res


@parameters.volume_of_distribution.formula
def v_d(dose, concentration):
    res = dose / concentration
    return res


if __name__ == '__main__':
    pass
