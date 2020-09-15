#
#  nuskell/dsdcompiler/interpreter.py
#  NuskellCompilerProject
#
# Copyright (c) 2010-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
import logging
log = logging.getLogger(__name__)

from . import parse_crn_string
from .ts_parser import parse_ts_string, parse_ts_file
from .interpreter import NuskellEnvironment 
from .objects import NuskellComplex

import os
import pkg_resources

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

def translate(input_crn, ts_file, modular = False):
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

    Returns:
      [:obj:`TestTube()`,...]: A list of TestTube objects.
      The first object contains signal and fuel species of the full DSD
      system, followed by the modular system specifications.
    """
    if not os.path.isfile(ts_file):
        builtin = '../schemes/literature/' + ts_file
        try:
            ts_file = pkg_resources.resource_filename('nuskell', builtin)
        except KeyError:
            schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
            raise InvalidSchemeError(ts_file, schemedir)

    print(ts_file)
    ts = parse_ts_file(ts_file)
    crn, fs = parse_crn_string(input_crn)
    solution, modules = compile_parsed(ts, crn, fs, modular = modular)
    return solution, modules

def ts_code_snippet():
    """ Builtin funtions for the nuskell language.

    Those functions are always automatically loaded upon initialization of the
    interpreter environment. 
    """
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

def compile_parsed(ts_parsed, crn_parsed, formals, modular = False, one = 100e-9):
    """ Translation of a CRN into a DSD system.

    Initializes the compiler environment, interprets the instructions of the
    translation scheme and then translates the CRN.

    Args:
        ts_parsed (List[List[...]]): Low-level instructions from the
            translation scheme, as returned from the **nuskell.parser** module.
        crn_parsed (List[List[...]]): List of list data structure for CRNs as
            returned from the **nuskell.parser** module.
        formals (Dict[name: trip]) : A dictionary of formal species names and their
            concentrations. Trip = (mode, value, unit)

    Returns:
        [dict,...]: Complexes and their concentrations.  If you have a modular
            representation, then the following dictionaries contain all the modules.
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
    fs_result = ts_env.translate_formal_species(list(formals.keys()))

    # translate the crn using the main() function
    cs_modules = ts_env.translate_reactions(crn_parsed, modular = modular)

    if not modular:
        assert len(cs_modules) == 1
    cs_solution = cs_modules[0]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Translation to final Complex objects using a new naming scheme.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    solution = dict()
    log.debug("Compiled species.")
    for k, v in fs_result.items():
        v.flatten_cplx  # NusComplex-specific function.
        v = NuskellComplex(v.sequence, v.structure, name = k)
        v.concentration = None #TODO
        log.debug(f"{type(v)} - {v} {v.kernel_string} {v.concentration}")
        solution[v.name] = v # flag that it is a signal?
        fs_result[k] = v

    num = 0
    for cplx in sorted(cs_solution):
        cplx.name = name = 'f' + str(num)
        assert cplx.concentration == ('constant', float('inf'), 'M')
        cplx.concentration = ('constant', 1*one, 'M')
        log.debug(f"{type(cplx)} - {cplx} {cplx.kernel_string} {cplx.concentration}")
        solution[cplx.name] = cplx # flag that it is a fuel?
        num += 1

    rxnmodules = []
    #for e, csm in enumerate(cs_modules[1:]):
    #    log.debug("MODULE {}".format(e))
    #    rxn = set(crn_parsed[e].reactants + crn_parsed[e].products)
    #    module = TestTube()

    #    # TODO select formal species?
    #    for k, v in fs_result.items():
    #        assert k == v.name
    #        if module.has_complex(v):
    #            raise ValueError("Overwriting existing module species")
    #        if not solution.has_complex(v):
    #            raise ValueError("Cannot find formal module species in solution.")
    #        if k in rxn:
    #            conc = fs[k] if None in fs[k] else (fs[k][1]*one, fs[k][0] == 'constant')
    #            module.add_complex(v, conc, ctype = 'signal')
    #        log.debug("{} - {} {}".format(type(v), v, v.kernel_string))

    #    for cplx in csm.complexes:
    #        if not solution.has_complex(cplx):
    #            raise ValueError("Cannot find constant module species in solution.")
    #        assert cs_solution.ReactionGraph.nodes[cplx]['ctype'] == 'fuel'
    #        module.add_complex(cplx, cs_solution.get_complex_concentration(cplx), ctype = 'fuel')
    #        log.debug("{} - {} {}".format(type(cplx), cplx, cplx.kernel_string))
    #    rxnmodules.append(module)
    return solution, rxnmodules



