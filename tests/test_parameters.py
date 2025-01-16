from pharmakin import parameters
from pharmakin.utils.units import ureg, Q_
from pharmakin.utils.utils import BulkImporter


def test_get_all_parameters(all_parameters):
    """Checks that every parameter subclass defined anywhere in the package is provided by get_all_parameters"""

    import pharmakin
    importer = BulkImporter(from_=pharmakin, child_of=parameters.Parameter, recurse_submodules=True)
    all_in_package = importer()
    
    assert set(all_in_package) == set(all_parameters)


def test_parameter_cannot_be_defined_without_unit():
    fail = False
    try:
        class ForgotUnit(parameters.Parameter):
            pass
    except RuntimeError:
        fail = True
    assert fail


def test_parameter_unit_conversion(all_parameters):
    for par in all_parameters:
        val = Q_(42.0, par.unit)
        try:
            par._validate_units(par.unit)
        except Exception:
            raise AssertionError(f"Invalid value ({val}) for parameter: {par}.")
    assert True


