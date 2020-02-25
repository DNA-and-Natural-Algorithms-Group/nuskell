#!/usr/bin/env python
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@dna.caltech.edu)
#
# Compiler module.
#

""" Wrapper functions used by the nuskell compiler script. """

import os
import logging
import pkg_resources

from nuskell.parser import parse_ts_file
from nuskell.crnutils import parse_crn_string
from nuskell.interpreter import interpret

class InvalidSchemeError(Exception):
    """Raise Error: Cannot find translation scheme."""

    def __init__(self, ts_file, builtin=None):
        self.message = "Cannot find translation scheme: {}\n".format(ts_file)

        if builtin:
            self.message += "You may want to use one of the built-in schemes instead:\n"
            self.message += "Schemes in {}:\n".format(builtin)
            for s in sorted(os.listdir(builtin)):
                self.message += " * {}\n".format(s)

        super(InvalidSchemeError, self).__init__(self.message)

def translate(input_crn, ts_file, modular = False, verbose = False):
    """CRN-to-DSD translation wrapper function.

    A formal chemical reaction network (CRN) is translated into a domain-level
    strand displacement (DSD) system. The translation-scheme and the CRN are
    parsed into low-level instructions using the **nuskell.parser** module,
    passed on to the **nuskell.interpreter** and returned in form of a
    **nuskell.objects.TestTube()** object.

    Args:
      input_crn (str): An input string representation of the formal CRN.
      ts_file (str): The input file name of a translation scheme.
      modular (bool, optional): Split CRN into modules.
      verbose (bool, optional): Print logging information during translation.
        Defaults to False.

    Returns:
      [:obj:`TestTube()`,...]: A list of TestTube objects.
      The first object contains signal and fuel species of the full DSD
      system, followed by the modular system specifications.
    """
    if not os.path.isfile(ts_file):
        builtin = 'schemes/' + ts_file
        try:
            ts_file = pkg_resources.resource_filename('nuskell', builtin)
            logging.info("Using scheme: {}".format(ts_file))
        except KeyError:
            schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
            raise InvalidSchemeError(ts_file, schemedir)

    ts = parse_ts_file(ts_file)
    crn, fs = parse_crn_string(input_crn)

    solution, modules = interpret(ts, crn, fs, modular = modular, verbose = verbose)

    return solution, modules

