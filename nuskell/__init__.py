#
#  nuskell/__init__.py
#  NuskellCompilerProject
#
from .dsdcompiler import __version__ # Let's use one combined version, for now.

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

