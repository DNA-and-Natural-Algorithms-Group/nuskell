# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@caltech.edu)
#
# The interpreter interface for translation schemes.
#

import sys
from copy import copy

from nuskell.parser import parse_ts_string
from nuskell.objects import TestTube, Complex, Domain
from nuskell.interpreter.environment import Environment 

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

def interpret(ts_parsed, crn_parsed, fs_list,
    sdlen=6, ldlen=15, verbose=False):
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
  ts_env = Environment('', sdlen=sdlen, ldlen=ldlen)

  # Parse a piece of sample code with common utilitiy functions
  header = parse_ts_string(ts_code_snippet())

  # add the code to the environment
  ts_env.interpret(header)

  # interpret the translation scheme
  ts_env.interpret(ts_parsed)

  # translate formal species list using the formal() function 
  fs_result = ts_env.translate_formal_species(fs_list)

  # translate the crn using the main() function 
  cs_solution = ts_env.translate_reactions(crn_parsed)

  solution = TestTube()
  for k,v in fs_result.items():
    v.flatten_cplx
    #print type(v), k, map(str, v.sequence), v.structure
    if k in map(str, solution.complexes) :
      raise ValueError("Overwriting existing name")
    new = Complex(name = k, 
        sequence = v.sequence, 
        structure = v.structure)
    solution.add_complex(new, (None, None), sanitycheck=True)

  #for k,v in solution.complexes.items():
  #  print type(v), v, map(str, v.sequence), v.structure

  num=1
  for cplx in sorted(cs_solution.complexes, key=lambda x : x.name):
    rename = 'f{}_'.format(str(num))
    cplx.name = rename
    ## Make sure nobody calls calls formal species like fuel species
    #if rename in map(str, solution.complexes) :
    #  raise ValueError("Duplicate fuel species name!")
    new = Complex(sequence = cplx.sequence, structure = cplx.structure, name = rename)
    solution.add_complex(new, cs_solution.get_complex_concentration(cplx), sanitycheck=True)
    num += 1

  return solution, cs_solution

