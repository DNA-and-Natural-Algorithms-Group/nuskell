# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@caltech.edu)
#
# The interpreter interface for translation schemes.
#
from __future__ import absolute_import, division, print_function

import sys
import logging

from nuskell.parser import parse_ts_string
from nuskell.objects import TestTube, NuskellComplex
from nuskell.interpreter.environment import NuskellEnvironment


def ts_code_snippet():
    # Builtin funtions extending the nuskell language, loaded upon initialization
    # of the interpreter environment.
    return """
    function range(x) = if x == 0 then [] else range(x - 1) + [x - 1] ;
    function sum(x) = if len(x) == 0 then empty elseif len(x) == 1 then x[0] else x[0] + sum(tail(x)) ;
    function len(x) = if x == [] then 0 else 1 + len(tail(x)) ;
    function reverse(x) = if x == [] then [] else reverse(tail(x)) + [x[0]] ;
    function rxn_degree(x, r) = if len(x) == 0 then [] elseif len(x[0].reactants) == r then [x[0]] + rxn_degree(tail(x), r) else rxn_degree(tail(x), r) ;
    function unirxn(x) = if len(x) == 0 then [] elseif len(x[0].reactants) == 1 then [x[0]] + unirxn(tail(x)) else unirxn(tail(x)) ;
    function birxn(x) = if len(x) == 0 then [] elseif len(x[0].reactants) == 2 then [x[0]] + birxn(tail(x)) else birxn(tail(x)) ;
    function map(f, x) = if len(x) == 0 then [] else [f(x[0])] + map(f, tail(x)) ;
    function map2(f, y, x) = if len(x) == 0 then [] else [f(y, x[0])] + map2(f, y, tail(x)) """


def interpret(ts_parsed, crn_parsed, fs, modular = False, verbose = False, one = 100e-9):
    """ Interpretation of a translation scheme.

    Initializes the compiler environment, interprets the instructions of the
    translation scheme and translates the CRN.

    Args:
      ts_parsed (List[List[...]]): Low-level instructions from the translation
        scheme, as returned from the **nuskell.parser** module.
      crn_parsed (List[List[...]]): List of list data structure for CRNs as
        returned from the **nuskell.parser** module.
      fs_list ([str,...]) : A list of formal species names.

      sdlen (int, optional): Domain length of the *built-in* short() function.
      ldlen (int, optional): Domain length of the *built-in* long() function.
      verbose (bool, optional): Print logging information during compilation.
        Defaults to False.

    Returns:
      [:obj:`TestTube()`,...]: A list of TestTube objects.
      The first object contains signal and fuel species of the full DSD
      system, followed by *experimental* modular system specifications.
    """

    # Initialize the environment
    ts_env = NuskellEnvironment()

    # Parse a piece of sample code with common utilitiy functions
    header = parse_ts_string(ts_code_snippet())

    # add the code to the environment
    ts_env.interpret(header)

    # interpret the translation scheme
    ts_env.interpret(ts_parsed)

    # translate formal species list using the formal() function
    fs_result = ts_env.translate_formal_species(fs.keys())

    # translate the crn using the main() function
    cs_modules = ts_env.translate_reactions(crn_parsed, modular = modular)

    if not modular:
        assert len(cs_modules) == 1
    cs_solution = cs_modules[0]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Translation to new Complex and TestTube objects ...
    #  ... using a new naming scheme
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    solution = TestTube()
    for k, v in fs_result.items():
        v.flatten_cplx  # NusComplex-specific function.
        v = NuskellComplex(v.sequence, v.structure, name = k)
        conc = fs[k] if None in fs[k] else (fs[k][1]*one, fs[k][0] == 'constant')
        solution.add_complex(v, conc, ctype = 'signal')
        fs_result[k] = v
        logging.debug("{} - {} {}".format(type(v), v, v.kernel_string))

    num = 0
    logging.debug("SOLUTION")
    for cplx in sorted(cs_solution.complexes, key = lambda x: x.name):
        cplx.name = name = 'f' + str(num)
        assert (cs_solution.get_complex_concentration(cplx)) == (float('inf'), None)
        assert cs_solution.ReactionGraph.nodes[cplx]['ctype'] == 'fuel'
        conc = (1*one, False)
        solution.add_complex(cplx, conc, ctype = 'fuel')
        logging.debug("{} - {} {}".format(type(cplx), cplx, cplx.kernel_string))
        num += 1

    rxnmodules = []
    for e, csm in enumerate(cs_modules[1:]):
        logging.debug("MODULE {}".format(e))
        rxn = set(crn_parsed[e][0] + crn_parsed[e][1])
        module = TestTube()

        # TODO select formal species?
        for k, v in fs_result.items():
            assert k == v.name
            if module.has_complex(v):
                raise ValueError("Overwriting existing module species")
            if not solution.has_complex(v):
                raise ValueError("Cannot find formal module species in solution.")
            if k in rxn:
                conc = fs[k] if None in fs[k] else (fs[k][1]*one, fs[k][0] == 'constant')
                module.add_complex(v, conc, ctype = 'signal')
            logging.debug("{} - {} {}".format(type(v), v, v.kernel_string))

        for cplx in csm.complexes:
            if not solution.has_complex(cplx):
                raise ValueError("Cannot find constant module species in solution.")
            assert cs_solution.ReactionGraph.nodes[cplx]['ctype'] == 'fuel'
            module.add_complex(cplx, cs_solution.get_complex_concentration(cplx), ctype = 'fuel')
            logging.debug("{} - {} {}".format(type(cplx), cplx, cplx.kernel_string))
        rxnmodules.append(module)

    return solution, rxnmodules
