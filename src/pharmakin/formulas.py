from __future__ import annotations

from pharmakin.utils.formulas import Formula
from pharmakin.utils import parameters


@parameters.clearance.formula
def clearence_from_dose_auc(dose, auc=0.3):
    """Determines clearence from dose and AUC"""
    res = dose / auc
    return res



if __name__ == '__main__':
    f = clearence_from_dose_auc
    dose = 0.5
    auc = 1.7
    res = f(dose, auc=auc)
    print(res)
    
    print(f.unit)
    print(f.eq)
    
    print(f)
    print(f.__doc__)
