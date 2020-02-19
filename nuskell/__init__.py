##############################################
# Initiation of the nuskell compiler package #
##############################################

__version__ = "v0.6"

# Import the main interface to the compiler #
from nuskell.compiler import translate, genCRN, genCON
from nuskell.verifier import verify
from nuskell.interpreter.environment import NuskellExit

