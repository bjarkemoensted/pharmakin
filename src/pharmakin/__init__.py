from importlib.metadata import version
import logging

__package_name__ = "pharmakin"
__version__ = version(__package_name__)

logging.getLogger(__name__).addHandler(logging.NullHandler())
