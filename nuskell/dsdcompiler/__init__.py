#
#  nuskell/dsdcompiler/__init__.py
#  NuskellCompilerProject
#
"""
A standalone sub-package that compiles CRNs into DSD systems using translation
schemes.  The sub-package contains parsing and interpretation of the nuskell
programming language, which provides a standardized formulation CRN to DSD
translation schemes.
"""
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from .ts_parser import parse_ts_file, parse_ts_string
#from .crn_parser import parse_crn_file, parse_crn_string
from nuskell.crnutils import parse_crn_file, parse_crn_string
from .compiler import translate

