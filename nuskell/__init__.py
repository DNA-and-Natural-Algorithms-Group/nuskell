##############################################
# Initiation of the nuskell compiler package #
##############################################

__version__ = "0.5"

# Import the main interface to the compiler #
from nuskell.compiler import translate, genCRN
from nuskell.verifier import verify
from nuskell.interpreter.environment import NuskellExit

