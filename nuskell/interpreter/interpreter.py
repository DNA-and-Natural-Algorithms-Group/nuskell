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
  """ Sample code of the ts language, which is parsed upon initialization of
  the ts environment.

  :return: a code snippet in the ts language format
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

def interpret(ts_parsed, crn_parsed, fs_list, cs_list,
    name='ts_author', sdlen=6, ldlen=15, verbose=False):
  """Interface to the nuskell environment.

  The translation scheme is executed in the nuskell environment.

  Args:
    ts_parsed (List[List[...]]) : nuskell code in a complex data structure as
      it is returned after parsing it with the **nuskell.parser** module.
    crn_parsed (Llist[List[...]]) : A data structure that describes a formal
      crn. See **nuskell.parser** for details.
    fs_list () : A list of formal species. 
    name (str) : A name for the nuskell environment, usually referring to the
      implemted translation scheme.
    sdlen (Optional[int]) : Domain length of the *built-in* short() function
    ldlen (Optional[int]) : Domain length of the *built-in* long() function
    .. verbose (Optional[bool]) : Toggle a verbose mode for additional feedback 
      during computations.

  Returns:
    The given CRN translated with the given translation scheme into a DSD circuit.
  """

  # Initialize the environment
  ts_env = Environment(name, sdlen=sdlen, ldlen=ldlen)

  # Parse a piece of sample code with common utilitiy functions
  header = parse_ts_string(ts_code_snippet())

  # add the code to the environment
  ts_env.interpret(header)

  # interpret the translation scheme
  ts_env.interpret(ts_parsed)

  # translate formal species list using the formal() function 
  fs_result = ts_env.translate_formal_species(fs_list)

  if cs_list :
    # translate a given constant species list using the constant() function 
    cs_solution = ts_env.translate_constant_species(cs_list, crn_parsed)
    cs_result = ts_env.constant_species_dict

  else :
    # translate the crn using the main() function 
    cs_solution = ts_env.translate_reactions(crn_parsed)
    cs_result = {}

  # # Alternative way to extract data at the end.
  # # get the results in form of a dictionary d={fs:Object}
  # fs_result = ts_env.formal_species_dict
  # # get the final solution object
  # cs_solution = ts_env.constant_species_solution

  solution = TestTube()
  for k,v in fs_result.items():
    v.flatten_cplx
    #print type(v), k, map(str, v.sequence), v.structure
    if k in map(str, solution.complexes) :
      raise ValueError("Overwriting existing name")
    new = Complex(name = k, 
        sequence = v.sequence, 
        structure = v.structure)
    solution.add_complex(new, (10, False), sanitycheck=True)

  #for k,v in solution.complexes.items():
  #  print type(v), v, map(str, v.sequence), v.structure

  num=1
  for cplx in cs_solution.complexes:
    rename = 'f{}_'.format(str(num))
    cplx.name = rename
    ## Make sure nobody calls calls formal species like fuel species
    #if rename in map(str, solution.complexes) :
    #  raise ValueError("Duplicate fuel species name!")
    new = Complex(sequence = cplx.sequence, structure = cplx.structure, name = rename)
    solution.add_complex(new, cs_solution.get_complex_concentration(cplx), sanitycheck=True)
    num += 1

  ##############################################################################
  # The following code is a greedy replacement of history-domains *before*
  # reaction enumeration.  This seems to work just fine, but there might be an
  # example where it makes a difference which produce domain is used to replace
  # the history domain?  
  # The primary reason why this is code is commented out is that the
  # verification preprocessing for history schemes relies on finding history
  # domains in the formal species, ... :-(

  #for nx, x in solution.complexes.items():
  #  hd = filter(lambda y:str(y)[0]=='h', x.sequence)
  #  if len(hd) == 1:
  #    if '+' in x.sequence:
  #      raise NotImplementedError("history domain in multistranded complex!")
  #    hd = hd[0]
  #    for ns, s in solution.strands.items():
  #      if x.sequence == s : 
  #        continue
  #      if patternMatch(x.sequence, s, ignore=str(hd)):
  #        dummy = Complex(s, x.structure, name=x.name)
  #        solution.rm_complex(x)
  #        solution.add_complex(dummy)
  #        break
  #  elif len(hd) > 1:
  #    raise NotImplementedError("multiple history domains in one complex!")

  return solution, cs_solution

def patternMatch(x, y, ignore = '?'):
  """Matches two complexes if they are the same, ignoring history domains. 

  Note: The strand order of the second complex changes to the strand order of
  the first complex, if there is a rotation under which both complexes are
  patternMatched.

  Args: 
    x (Complex()) : A nuskell Complex() object.
    y (Complex()) : A nuskell Complex() object.

  Returns: True/False
  """
  if len(x) != len(y) :
    return False

  def pM_check(pMx, pMy):
    """Recursively parse the current sequences and structures. 

    Args: 
      pMx [seqX,strX]: A list of two lists (sequence, structrure)
      pMy [seqY,strY]: A list of two lists (sequence, structrure)

    Returns: True/False
    """
    if len(pMx) == 0 :
      return True

    if pMx[0] != ignore and pMy[0] != ignore and pMx[0] != pMy[0]:
      return False
    return pM_check(pMx[1:], pMy[1:])

  pMx = map(str, x)
  pMy = map(str, y)
  if pM_check(pMx,pMy) :
    return True
  return False

