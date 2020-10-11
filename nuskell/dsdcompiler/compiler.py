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

import os
import pkg_resources

from .crn_parser import parse_crn_file, parse_crn_string
from .ts_parser import parse_ts_string, parse_ts_file
from .interpreter import NuskellEnvironment 
from .objects import NuskellComplex
SCHEME_DIRS = ['schemes/literature/', 'schemes/variants/'] 

class NuskellInterpreterError(Exception):
    pass

class InvalidSchemeError(Exception):
    def __init__(self, ts_file):
        self.message = f"Cannot find translation scheme: {ts_file}\n"
        schemes = get_builtin_schemes()
        for k, v in schemes.items():
            log.error(f'Listing installed schemes: "{k}"')
            for s in v:
                log.error(f"   --ts {s}")
        super(InvalidSchemeError, self).__init__(self.message)

def get_builtin_schemes():
    global SCHEME_DIRS
    schemes = dict()
    for d in SCHEME_DIRS:
        builtin = pkg_resources.resource_filename(__name__, d)
        schemes[builtin] = sorted(os.listdir(builtin))
    return schemes

def get_canonical_schemes():
    # Candidates for "default" translation schemes.
    return {'canonical': ['soloveichik2010.ts',
                          'qian2011_3D_var1.ts', 
                          'cardelli2011_NM.ts',
                          'cardelli2011_NM_noGC.ts',
                          'lakin2012_3D_var1.ts',
                          'cardelli2013_2D_3I_noGC.ts',
                          'chen2013_2D_JF_var1.ts',
                          'srinivas2015.ts']}

def find_scheme_file(ts):
    global SCHEME_DIRS
    if not os.path.isfile(ts):
        for d in SCHEME_DIRS:
            builtin = d + ts
            try:
                tsx = pkg_resources.resource_filename(__name__, builtin)
                if os.path.exists(tsx):
                    ts = tsx
                    break
            except KeyError as err:
                pass
        else:
            raise InvalidSchemeError(ts)
    return ts

def translate(input_crn, ts_file, modular = False):
    """ CRN-to-DSD translation wrapper function.

    A formal chemical reaction network (CRN) is translated into a domain-level
    strand displacement (DSD) system. The translation-scheme and the CRN are
    parsed into low-level instructions, passed on to the **interpreter** and
    returned in form of a list of complexes.

    Args:
      input_crn (str): An input string representation of the formal CRN.
      ts_file (str): The input file name of a translation scheme.
      modular (bool, optional): Split CRN into modules.

    Returns:
      [:obj:`TestTube()`,...]: A list of TestTube objects.
      The first object contains signal and fuel species of the full DSD
      system, followed by the modular system specifications.
    """
    ts_file = find_scheme_file(ts_file)
    ts = parse_ts_file(ts_file)
    crn, fs = parse_crn_string(input_crn)
    solution, modules = interpret(ts, crn, fs, modular = modular)
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

def interpret(ts_parsed, crn_parsed, formals, modular = False, one = 100):
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

    for name in formals.keys():
        if name[0] == 'f': # fuels
            raise NuskellInterpreterError(f'Formal species name must not start with "f": {name}.')
        elif name[0] == 'i': # intermediates
            raise NuskellInterpreterError(f'Formal species name must not start with "i": {name}.')
        elif name[0] == 'w': # wastes
            raise NuskellInterpreterError(f'Formal species name must not start with "w": {name}.')

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
    log.debug("Compiled signal species:")
    for k, v in fs_result.items():
        v.flatten_cplx  # NusComplex-specific function.
        cplx = NuskellComplex(v.sequence, v.structure, name = k)
        cplx.concentration = None if None in formals[k] else (*formals[k], 'nM')
        log.debug(f"{type(cplx)} - {cplx} {cplx.kernel_string} {cplx.concentration}")
        solution[cplx.name] = cplx
        fs_result[k] = cplx

    num = 1
    log.debug("Compiled fuel species:")
    for cplx in sorted(cs_solution):
        cplx.name = 'f' + str(num)
        assert cplx.concentration == ('constant', float('inf'), 'nM')
        cplx.concentration = ('constant', 1*one, 'nM')
        log.debug(f"{type(cplx)} - {cplx} {cplx.kernel_string} {cplx.concentration}")
        solution[cplx.name] = cplx
        num += 1

    rxnmodules = []
    for e, csm in enumerate(cs_modules[1:]):
        module = dict()
        log.debug(f"Compiled module {e}")
        rxn = set(crn_parsed[e].reactants + crn_parsed[e].products)
        for k, cplx in fs_result.items():
            assert k == cplx.name
            if cplx.name in module:
                raise NuskellInterpreterError("Overwriting existing signal species in module.")
            if cplx.name not in solution:
                raise NuskellInterpreterError("Cannot find signal species of module in solution.")
            if k in rxn:
                module[cplx.name] = cplx
            log.debug(f"{type(cplx)} - {cplx} {cplx.kernel_string} {cplx.concentration}")

        for cplx in sorted(csm):
            if cplx.name not in solution:
                raise NuskellInterpreterError("Cannot find fuel species of module in solution.")
            module[cplx.name] = cplx
            log.debug(f"{type(cplx)} - {cplx} {cplx.kernel_string} {cplx.concentration}")
        rxnmodules.append(module)
    return solution, rxnmodules

if __name__ == '__main__':
   translate()

