#
#  nuskell/interpreter/__init__.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

import sys
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from nuskell.interpreter.interpreter import interpret
