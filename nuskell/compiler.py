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
import pkg_resources

from nuskell.parser import parse_crn_string, parse_ts_file
from nuskell.parser import split_reversible_reactions
from nuskell.parser import combine_reversible_reactions

from nuskell.interpreter import interpret
from nuskell.objects import TestTube, TestTubeIO


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


def printCRN(crn, reversible=True, rates=True):
    """Pretty printing of CRNs.

    Args:
      crn (list of lists): A CRN in list of list format.
      reversible (bool, optional): True to combine reversible reactions into one
        line.  False to split reversible reactions into irreversible reactions.

    Note:
      Order of reactants/products may differ from the order in the crn argument.
    """
    if not rates:
        pcrn = map(lambda rxn: rxn[:2] + [[None]], crn)
    else:
        pcrn = [r for r in crn]

    if reversible:
        pcrn = combine_reversible_reactions(pcrn)
    else:
        pcrn = split_reversible_reactions(pcrn)

    for rxn in pcrn:
        assert len(rxn) == 3
        if len(rxn[2]) == 2:
            print ' + '.join(rxn[0]), '<=>', ' + '.join(rxn[1])
        else:
            print ' + '.join(rxn[0]), '->', ' + '.join(rxn[1])


def translate(input_crn, ts_file, modular=False,
              pilfile=None, dnafile=None, verbose=False):
    """CRN-to-DSD translation wrapper function.

    A formal chemical reaction network (CRN) is translated into a domain-level
    strand displacement (DSD) system. The translation-scheme and the CRN are
    parsed into low-level instructions using the **nuskell.parser** module,
    passed on to the **nuskell.interpreter** and returned in form of a
    **nuskell.objects.TestTube()** object.

    Args:
      input_crn (str): An input string representation of the formal CRN.
      ts_file (str): The input file name of a translation scheme.
      pilfile (str, optional): Prints the DSD system in form of a PIL file.
        Defaults to None.
      dnafile (str, optional): Prints the DSD system in form of a VisualDSD DNA
        file. Defaults to None.
      verbose (bool, optional): Print logging information during translation.
        Defaults to False.

    Returns:
      [:obj:`TestTube()`,...]: A list of TestTube objects.
      The first object contains signal and fuel species of the full DSD
      system, followed by *experimental* modular system specifications.

    """

    TestTube.warnings = verbose

    if not os.path.isfile(ts_file):
        builtin = 'schemes/' + ts_file

        try:
            ts_file = pkg_resources.resource_filename('nuskell', builtin)
            print "Using scheme:", ts_file
        except KeyError:
            schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
            raise InvalidSchemeError(ts_file, schemedir)

    ts = parse_ts_file(ts_file)
    crn, fs, signals, fuels = parse_crn_string(input_crn)

    solution, modules = interpret(ts, crn, fs, modular=modular)

    if pilfile:
        with open(pilfile, 'w') as pil:
            TestTubeIO(solution).write_pil_kernel(pil)
    if dnafile:
        with open(dnafile, 'w') as dna:
            TestTubeIO(solution).write_dnafile(
                dna, signals=fs, crn=crn, ts=os.path.basename(ts_file))

    return solution, modules
