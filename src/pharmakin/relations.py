from pharmakin import parameters
from pharmakin.utils.registry import Formulary
from pharmakin import formulas

ALL_PARAMETERS = parameters.get_all_parameters()
ALL_FORMULAS = formulas.get_all_formulas()


formulary = Formulary()
formulary.register_parameters(ALL_PARAMETERS)
formulary.register_formulas(ALL_FORMULAS)


if __name__ == '__main__':
    print(formulary)
    print(formulary.equations)