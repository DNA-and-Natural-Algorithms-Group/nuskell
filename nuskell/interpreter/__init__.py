#
#  nuskell/interpreter/__init__.py
#  NuskellCompilerProject
#
import sys
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from interpreter import interpret
