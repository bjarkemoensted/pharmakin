from pharmakin import parameters
from pharmakin.utils.registry import Formulary
from pharmakin import formulas

ALL_PARAMETERS = parameters.get_all_parameters()
ALL_FORMULAS = formulas.get_all_formulas()


formulary = Formulary(
    parameters=ALL_PARAMETERS,
    formulas=ALL_FORMULAS
)

if __name__ == '__main__':
    dose = 20
    dose = parameters.dose.ensure_units(dose)
    
    clearance = 1
    clearance = parameters.clearance.ensure_units(clearance)
    hmm = formulary.determine_parameter(
        parameter = parameters.auc,
        dose=dose,
        clearance=clearance
    )
    print(hmm)
    