#
#  nuskell/verifier/__init__.py
#  NuskellCompilerProject
#
import sys
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from verifier import verify, modular_bisimulation
