from pharmakin import parameters
from pharmakin.units import ureg, Q_


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
    