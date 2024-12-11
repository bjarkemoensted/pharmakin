from pharmakin.units import Q_, ureg, has_unit


def test_has_unit():
    a = Q_(42.0, ureg.m)
    b = Q_(1337.0)
    
    assert has_unit(a)
    assert not has_unit(b)


def test_preferred_units_dont_overlap():
    """Check that the specified preferred units are orthogonal, i.e. we don't prefer e.g. both meters and cm"""
    pref = ureg.default_preferred_units
    for i, unit in enumerate(pref):
        remaining = pref[i+1:]
        assert not any(unit.is_compatible_with(other) for other in remaining)
    #
