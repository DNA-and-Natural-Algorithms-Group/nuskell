#!/usr/bin/env python
#
#  nuskell/verifier/verifier.py
#  NuskellCompilerProject
#
# Copyright (c) 2009-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
import logging
log = logging.getLogger(__name__)

import signal
from collections import Counter

import nuskell.verifier.crn_pathway_equivalence 
import nuskell.verifier.crn_bisimulation_equivalence 

class TimeoutError(Exception):
    pass

def handler(signum, frame):
    raise TimeoutError('Time over')

def verify(formal_crn, impl_crn, formals, interpret = None,
            method = 'bisimulation', timeout = 0):
    """Verify the equivalence of a formal CRN and its implementation CRN.

    This wrapper function for two notions of equivalence (bisimulation and
    pathway decomposition) and variations of in their implementation.

    Args:
      formal_crn (list[list[...]]): List of list data structure for formal CRNs
        as returned from the **nuskell.parser** module. The CRN must specify
        reversible reactions as two irreversible reactions.
      impl_crn (list[list[...]]): List of list data structure for implementation
        CRNs as returned from the **nuskell.parser** module. The CRN must specify
        reversible reactions as two irreversible reactions.
      formals (list[str,...]): A list of formal species.
      interpret (:obj:`dict`, optional): A dictionary storing a partial
        interpretation of the System. Typically stores a mapping between formal
        and singal species.
      method (str, optional): Choose a notion of equivalence: 'bisimulation',
        'pathway', 'integrated'.
      timeout (int, optional): Set a timeout (in seconds) for verification. Defaults
        to 0, i.e. no timeout.

    Returns:
      bool: True if equivalent, False otherwise.

    """
    v, i = None, None

    for rxn in formal_crn + impl_crn:
        if rxn.k_rev != 0:
            raise Exception('Reaction does not have the required format:', rxn)

    fcrn = [[Counter(part) for part in rxn[:2]] for rxn in formal_crn]
    ecrn = [[Counter(part) for part in rxn[:2]] for rxn in impl_crn]

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout)
    try:
        if method == 'bisimulation' or method == 'bisim-whole-graph':
            v, i = nuskell.verifier.crn_bisimulation_equivalence.test(fcrn, ecrn, formals,
                    interpretation=interpret, permissive='whole-graph')

        elif method == 'bisim-loop-search':
            v, i = nuskell.verifier.crn_bisimulation_equivalence.test(fcrn, ecrn, formals,
                    interpretation=interpret, permissive='loop-search')

        elif method == 'bisim-depth-first':
            v, i = nuskell.verifier.crn_bisimulation_equivalence.test(fcrn, ecrn, formals,
                    interpretation=interpret, permissive='depth-first')

        elif method == 'pathway':
            # NOTE: Adaptation to pathway interface
            pinter = dict()
            if interpret:
                for k, v in interpret.items():
                    pinter[k] = sorted(v.elements())
            v = nuskell.verifier.crn_pathway_equivalence.test((formal_crn, formals),
                                                              (impl_crn, pinter.keys()), 
                                                              inter = pinter,
                                                              integrated = False)
        elif method == 'integrated':
            # TODO: integrated-hybrid -> first consider some species as formal for
            # pathway decomposition, then do bisimulation. This is necessary for
            # history domains in some schemes, but it can be used for more general
            # things. E.g. if you make reversible reactions formal, then it will get
            # accepted and bisimulation can do the rest. In any case, the current
            # implementation does not cover the full hybrid theory, only some special
            # cases and some kind of bisimulation.
            pinter = dict()
            if interpret:
                for k, v in interpret.items():
                    pinter[k] = sorted(v.elements())
            v = nuskell.verifier.crn_pathway_equivalence.test((formal_crn, formals),
                                                              (impl_crn, pinter.keys()), 
                                                              inter = pinter,
                                                              integrated = True)
        else:
            raise RuntimeError('Unknown verification method.')
    except TimeoutError:
        v = None
    finally:
        signal.alarm(0)

    return v, i

def modular_bisimulation(formal_crn, impl_crn, formals, interpret = None,
                         method = 'bisimulation', timeout = 0):
    """
      Wrapper to choose from different notions of modular bisimulation equivalence.
    """

    v, i = None, None

    fcrns = [[[Counter(part) for part in rxn[:2]] for rxn in mod]
             for mod in formal_crn]
    ecrns = [[[Counter(part) for part in rxn[:2]] for rxn in mod]
             for mod in impl_crn]

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout)
    try:
        # TODO: replace second formals!!!
        if method == 'bisimulation' or method == 'bisim-whole-graph':
            v, i = nuskell.verifier.crn_bisimulation_equivalence.testModules(fcrns, ecrns, formals,
                                                            interpretation=interpret, permissive='whole-graph')

        elif method == 'bisim-loop-search':
            v, i = nuskell.verifier.crn_bisimulation_equivalence.testModules(fcrns, ecrns, formals,
                                                            interpretation=interpret, permissive='loop-search')

        elif method == 'bisim-depth-first':
            v, i = nuskell.verifier.crn_bisimulation_equivalence.testModules(fcrns, ecrns, formals,
                                                            interpretation=interpret, permissive='depth-first')

        else:
            raise RuntimeError('Unsupported verification method.')
    except TimeoutError:
        v = None
    finally:
        signal.alarm(0)

    return v, i
