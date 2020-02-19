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

from dsdobjects.utils import natural_sort
from crnsimulator import parse_crn_string

from nuskell.parser import parse_ts_file
from nuskell.parser import split_reversible_reactions
from nuskell.parser import combine_reversible_reactions
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

def genCON(species):
    for sp, data in natural_sort(species.items()):
        if None in data: continue
        yield '{} @{} {}'.format(sp, 'initial' if data[0][0] == 'i' else 'constant', data[1])

def genCRN(crn, reversible = True, rates = True):
    """ Pretty printing of CRNs.

    Args:
        crn (list of lists): A CRN in list of list format.
        reversible (bool, optional): True to combine reversible reactions into one
            line.  False to split reversible reactions into irreversible reactions.
        rates (bool, optional): 

    Note:
      Order of reactants/products may differ from the order in the crn argument.
    """
    if rates:
        pcrn = [r for r in crn]
    else:
        pcrn = map(lambda rxn: rxn[:2] + [[None]], crn)

    if reversible:
        pcrn = combine_reversible_reactions(pcrn)
    else:
        pcrn = split_reversible_reactions(pcrn)

    for rxn in pcrn:
        assert len(rxn) == 3
        if len(rxn[2]) == 2:
            yield '{} <=> {} [kf = {:g}, kr = {:g}]'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]), float(rxn[2][0]), float(rxn[2][1]))
        else:
            yield '{} -> {} [k = {:g}]'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]), float(rxn[2][0]))

    #for rxn in pcrn:
    #    assert len(rxn) == 3
    #    if len(rxn[2]) == 2:
    #        yield '{} <=> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]))
    #    else:
    #        yield '{} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]))


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
    crn, fs = post_process(parse_crn_string(input_crn, process = False))

    solution, modules = interpret(ts, crn, fs, modular = modular, verbose = verbose)

    return solution, modules


def post_process(crn):
    """Process a parsed CRN.
    Taken and modified from crnsimulator, drop for Python 3.x
    """
    def remove_multipliers(species):
        flat = []
        for s in species:
            if len(s) == 1:
                flat.append(s[0])
            elif len(s) == 2:
                ss = [s[1]] * int(s[0])
                flat.extend(ss)
        return flat

    new = []
    species = dict()
    for line in crn:
        if line[0] == 'concentration':
            spe = line[1][0]
            ini = 'initial' if line[2][0][0] == 'i' else 'constant'
            num = line[3][0]
            species[spe] = (ini, float(num))
            continue
        elif len(line) == 3:
            # No rate specified
            t, r, p = line
            r = remove_multipliers(r)
            p = remove_multipliers(p)
            if t == 'reversible':
                new.append([r, p, [1., 1.]])
            elif t == 'irreversible':
                new.append([r, p, [1.]])
            else:
                raise CRNParseError('Wrong CRN format!')
        elif len(line) == 4:
            t, r, p, k = line
            r = remove_multipliers(r)
            p = remove_multipliers(p)
            if t == 'reversible':
                assert len(k) == 2
                new.append([r, p, k])
            elif t == 'irreversible':
                assert len(k) == 1
                new.append([r, p, k])
            else:
                raise CRNParseError('Wrong CRN format!')
        else:
            raise CRNParseError('Wrong CRN format!')
        for s in r + p:
            if s not in species:
                species[s] = (None, None)
    return new, species


