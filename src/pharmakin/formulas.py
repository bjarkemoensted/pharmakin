import sys
_this_module = sys.modules[__name__]

from pharmakin import parameters
from pharmakin.utils.formula import Formula
from pharmakin.utils.formula import Formulary
from pharmakin.utils.utils import BulkImporter


from pharmakin import kinetics
get_all_formulas = BulkImporter(from_=kinetics, instance_of=Formula, recurse_submodules=True)


ALL_PARAMETERS = parameters.get_all_parameters()
ALL_FORMULAS = get_all_formulas()


formulary = Formulary(
    parameters=ALL_PARAMETERS,
    formulas=ALL_FORMULAS
)


if __name__ == '__main__':
    pass
