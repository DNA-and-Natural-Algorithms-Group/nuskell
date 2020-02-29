#
#  __init__.py
#  NuskellCompilerProject
#
__version__ = "v0.6"

import sys
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from nuskell.compiler import translate
from nuskell.verifier import verify
from nuskell.interpreter.environment import NuskellExit

