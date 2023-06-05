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
from .crnutils import split_reversible_rxns
from crnverifier import (pathway_decomposition_eq,
                         crn_bisimulation_test, 
                         integrated_hybrid_test,
                         compositional_hybrid_test,
                         modular_crn_bisimulation_test)

from threading import Thread
import functools

def timer_func(timeout):
    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            res = [Exception('function [%s] timeout [%s seconds] exceeded!' % (func.__name__, timeout))]
            def newFunc():
                try:
                   res[0] = func(*args, **kwargs)
                except Exception as e:
                    res[0] = e

            t = Thread(target=newFunc)
            t.daemon = True
            try:
                t.start()
                t.join(timeout)
            except Exception as je:
                print('error starting thread')
                raise je
            ret = res[0]
            if isinstance(ret, BaseException):
                raise ret
            return ret
        return wrapper
    return deco

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
    fcrn = [list(rxn[:2]) for rxn in split_reversible_rxns(fcrn)]
    icrn = [list(rxn[:2]) for rxn in split_reversible_rxns(icrn)]

    timer_func(timeout)

    try:
        if 'crn-bisimulation' in method:
            if method[-2:] == '-ls':
                permissive_check = 'loopsearch'
            elif method[-2:] == '-bf':
                permissive_check = 'bruteforce'
            else:
                permissive_check = 'graphsearch'
            v, i = crn_bisimulation_test(fcrn,
                                         icrn,
                                         formals,
                                         interpretation = interpretation,
                                         permissive = permissive_check)
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
    except Exception:
        v, i = None, None
    finally:
        return v, i

def verify_modules(fcrns, icrns, formals, method, interpretation = None, timeout = 0):
    """ Choose from different algorithms for modular CRN bisimulation. """
    fcrns = [[list(rxn[:2]) for rxn in split_reversible_rxns(mod)] for mod in fcrns]
    icrns = [[list(rxn[:2]) for rxn in split_reversible_rxns(mod)] for mod in icrns]

    timer_func(timeout)
    try:
        if 'crn-bisimulation' in method:
            if method[-2:] == '-ls':
                permissive_check = 'loopsearch'
            elif method[-2:] == '-bf':
                permissive_check = 'bruteforce'
            else:
                permissive_check = 'graphsearch'
            v, i = modular_crn_bisimulation_test(fcrns,
                                                 icrns,
                                                 formals,
                                                 interpretation = interpretation,
                                                 permissive = permissive_check)
        else:
            raise RuntimeError('Unsupported verification method.')
    except TimeoutError:
        v, i = None, None
    finally:
        return v, i
