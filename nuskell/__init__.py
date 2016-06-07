##############################################
# Initiation of the nuskell compiler package #
##############################################

__version__ = "0.0.1"

# Import the compiler I/O base #
from nuskell.compiler import main, compile

import nuskell.parser

import nuskell.interpreter.interpreter as interpreter

# Import verification utilities #
# import nuskell.verifier.verifier as verifier
# import nuskell.verifier.basis_finder as basis_finder
# import nuskell.verifier.crn_pathway_equivalence as crn_pathway_equivalence 
# import nuskell.verifier.crn_bisimulation_equivalence as crn_bisimulation_equivalence

