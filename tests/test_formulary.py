from pharmakin.utils.units import has_units


def test_all_parameters_defined(master_formulary, all_parameters):
    """Checks that for each parameter, at least one formula containing the parameter has been defined"""
    
    missing = []
    for par in all_parameters:
        if len(master_formulary.get_formulas_containing_parameter(par)) == 0:
            missing.append(par)
        #
    assert not missing


def test_computations(master_formulary, parameter_example_values):
    for formula in master_formulary.formulae:
        for d in parameter_example_values:
            kws = {parname: d[parname] for parname in formula.inputs}
            res = formula(**kws)
            assert formula.result_class.is_valid(res)
            assert not has_units(res)
            
            kws_with_units = {k: master_formulary.parameters[k].ensure_units(v) for k, v in kws.items()}
            res_with_units = formula(**kws_with_units)
            assert formula.result_class.is_valid(res_with_units)
            