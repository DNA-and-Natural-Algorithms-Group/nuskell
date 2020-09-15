#
#  nuskell/__init__.py
#  NuskellCompilerProject
#
from .dsdcompiler import __version__

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

