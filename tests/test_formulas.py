from pharmakin.utils.formula import Formula
from pharmakin.utils.utils import BulkImporter


def test_get_all_formulas(all_formulas):
    """Checks that all formulas defined anywhere is returned by get_all_formulas"""

    import pharmakin
    importer = BulkImporter(from_=pharmakin, instance_of=Formula, recurse_submodules=True)
    
    all_in_package = importer()
    assert set(all_formulas) == set(all_in_package)
