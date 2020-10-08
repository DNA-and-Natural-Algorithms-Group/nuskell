#
#  nuskell/dsdcompiler/__init__.py
#  NuskellCompilerProject
#
__version__ = "0.8"
"""
A sub-package that compiles CRNs into DSD systems using translation schemes.

This sub-package contains the parser and the interpreter for the nuskell
programming language (a language to formulate CRN to DSD translation schemes).
"""
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from .compiler import translate
from .interpreter import NuskellExit

