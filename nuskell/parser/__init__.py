#
#  nuskell/parser/__init__.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from nuskell.parser.ts_parser import parse_ts_file, parse_ts_string

# DEPRECATED
from nuskell.parser.crn_parser import parse_crn_file, parse_crn_string
