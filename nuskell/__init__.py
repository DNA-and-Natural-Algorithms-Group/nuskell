#
#  __init__.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

__version__ = "v0.6"

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from nuskell.compiler import translate
from nuskell.verifier import verify
from nuskell.interpreter.environment import NuskellExit

