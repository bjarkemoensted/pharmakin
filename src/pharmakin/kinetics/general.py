from pharmakin import parameters


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
