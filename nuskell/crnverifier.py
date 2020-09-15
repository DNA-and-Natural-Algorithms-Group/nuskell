#!/usr/bin/env python
#
#  nuskell/verifier.py
#  NuskellCompilerProject
#
# Copyright (c) 2009-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
import logging
log = logging.getLogger(__name__)

import signal
from .crnutils import split_reversible_reactions as split
from collections import Counter
from crnverifier import (pathway_decomposition_eq,
                         crn_bisimulation_test, 
                         integrated_hybrid_test,
                         compositional_hybrid_test,
                         modular_crn_bisimulation_test)

class TimeoutError(Exception):
    pass

def handler(signum, frame):
    raise TimeoutError('Time over')

def verify(fcrn, icrn, formals, method, interpretation = None, timeout = 0):
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
    fcrn = [list(rxn[:2]) for rxn in split(fcrn)]
    icrn = [list(rxn[:2]) for rxn in split(icrn)]

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout)
    try:
        if 'crn-bisimulation' in method:
            if method[-2:] == '-ls':
                permissive_check = 'whole-graph'
            elif method[-2:] == '-rs':
                permissive_check = 'depth-first'
            else:
                permissive_check = 'whole-graph'
            v, i = crn_bisimulation_test(fcrn, 
                                         icrn, 
                                         formals,
                                         interpretation = interpretation,
                                         permissive = permissive_check) 
            if v:
                i = {k: list(v.elements()) for k, v in i.items()}
            else:
                i = {k: list(v.elements()) for k, v in i[0].items()}
        elif method == 'pathway-decomposition':
            v = pathway_decomposition_eq([fcrn, icrn], formals)
            i = None
        elif method == 'compositional-hybrid':
            v, i = compositional_hybrid_test(fcrn, 
                                             icrn, 
                                             formals,
                                             interpretation)
        elif method == 'integrated-hybrid':
            v, i = integrated_hybrid_test(fcrn, 
                                          icrn, 
                                          formals,
                                          interpretation)
        else:
            raise RuntimeError('Unknown verification method.')
    except TimeoutError:
        v, i = None, None
    finally:
        signal.alarm(0)
    return v, i

def verify_modules(fcrns, icrns, formals, method, interpretation = None, timeout = 0):
    """ Choose from different algorithms for modular CRN bisimulation. """
    fcrns = [[list(rxn[:2]) for rxn in split(mod)] for mod in fcrns]
    icrns = [[list(rxn[:2]) for rxn in split(mod)] for mod in icrns]

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout)
    try:
        if 'crn-bisimulation' in method:
            if method[-2:] == '-ls':
                permissive_check = 'whole-graph'
            elif method[-2:] == '-rs':
                permissive_check = 'depth-first'
            else:
                permissive_check = 'whole-graph'
            v, i = modular_crn_bisimulation_test(fcrns, 
                                                 icrns, 
                                                 formals,
                                                 interpretation = interpretation, 
                                                 permissive = permissive_check)
            if v:
                i = {k: list(v.elements()) for k, v in i.items()}
            else:
                i = {k: list(v.elements()) for k, v in i[2].items()}
        else:
            raise RuntimeError('Unsupported verification method.')
    except TimeoutError:
        v, i = None, None
    finally:
        signal.alarm(0)
    return v, i
