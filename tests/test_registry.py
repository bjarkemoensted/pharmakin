from pharmakin.registry import Registry


def test_double_registration_fails():
    """Check that it's not possible to overwrite an existing registration"""
    r = Registry()
    key = "foo"
    fail = False
    r[key] = "bar"
    try:
        r[key] = "baz"
    except RuntimeError:
        fail = True
    assert fail
