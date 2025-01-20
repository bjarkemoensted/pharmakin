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
        assert par._unit_is_valid(par.unit)
    #


def test_example_values_valid(all_parameters):
    """Checks that parameter example values are valid with and without units"""
    for par in all_parameters:
        for with_units in (False, True):
            val = par.example_values(with_units=with_units)
            assert par.is_valid(val)
        #
    #


def test_validation_on_compatible_units(all_parameters):
    """Checks that valid values for each parameters are also valid after converting to compatible units"""
    for par in all_parameters:
        # Get the compatible units (this gives an empty list for dimentionless parameters)
        alternative_units = sorted(ureg.get_compatible_units(par.unit), key=str)
        val = par.example_values(with_units=True, seed=42)
        for alt in alternative_units:
            converted = val.to(alt)
            assert par.is_valid(converted)
        