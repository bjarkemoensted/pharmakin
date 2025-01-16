import contextlib
from importlib import import_module
import os
from pathlib import Path
from types import ModuleType
from typing import Iterable


def _iterate_modules(module):
    """Takes a module and recursively iterates through submodules"""
    
    # Take the directory of the provided module as the root
    root_dir = Path(module.__file__).parent
    # Get the fully qualified name for the input module
    base = module.__name__.split(".")
    
    for path in root_dir.rglob("*.py"):
        relpath = path.relative_to(root_dir)
        filename = relpath.stem
        if filename == "__init__":
            continue
        
        # Construct the module path (e.g. "module.submodule.subsubmodule")
        parts = base + list(relpath.parent.parts) + [filename]
        modpath = ".".join(parts)
        
        # Do the import (silence any print statements)
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                submodule = import_module(modpath)
            #
        yield submodule
    #


class BulkImporter:
    """Helper callable for dynamically bulk importing specific classes or instances"""
    def __init__(
            self,
            from_: Iterable|ModuleType,
            instance_of: Iterable|type=None,
            child_of: Iterable|type=None,
            recurse_submodules=False
            ):
        """from_ is a module, or an iterable of modules.
        instance_of denotes one or more classes of which instances will be imported.
        child_of denotes one or more classes of which child classes will be imported.
        The class itself is not included. For either argument, None can be used to allow all.
        If recurse_submodules is True, all submodules are scanned as well."""
        
        # Cast all arguments to tuples for consistency
        if isinstance(from_, ModuleType):
            from_ = (from_,)
        if isinstance(instance_of, type):
            instance_of = (instance_of,)
        if isinstance(child_of, type):
            child_of = (child_of,)
        
        self._modules = tuple(elem for elem in from_)
        self.instance_of = instance_of
        self.child_of = child_of
        self.recurse_submodules = recurse_submodules
    
    @property
    def from_(self):
        mods = set()
        
        for mod in self._modules:
            mods.add(mod)
            if self.recurse_submodules:
                for submodule in _iterate_modules(mod):
                    mods.add(submodule)
                #
            #
        
        res = sorted(mods, key=str)
        return res
    
    def _instance_match(self, attribute):
        """Checks if the attribute is an instance of any specified classes"""
        if self.instance_of is None:
            return True
        
        return isinstance(attribute, self.instance_of)
    
    def _class_match(self, attribute):
        """Checks if the attribute is a child class of and specified class"""
        if self.child_of is None:
            return True
        
        # If attribute is not a class, it can't be a child class (and would cause issubclass to crash)
        if not isinstance(attribute, type):
            return False
        
        is_subclass = issubclass(attribute, self.child_of)
        exact_match = any(attribute is class_ for class_ in self.child_of)
        res = is_subclass and not exact_match
        
        return res
    
    def match(self, attribute):
        """Checks if an attribute in a module should be imported"""
        
        res = self._instance_match(attribute) and self._class_match(attribute)
        return res
    
    def __call__(self):
        hits = []
    
        # Look for instances of the desired class(es) in all modules
        for module in self.from_:
            for k, v in module.__dict__.items():
                elem = getattr(module, k)
                if self.match(elem):
                    hits.append(v)
                #
            #
        
        res = tuple(hits)
        return res
    #


if __name__ == '__main__':
    pass
