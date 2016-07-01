#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# The interpreter for translation schemes.
#

import sys
from copy import copy

import dnaobjects as DNAObjects

from nuskell.parser import parse_ts_string
from nuskell.interpreter.environment import Environment, Structure
from nuskell.interpreter.environment import builtin_functions


def flatten(x):
  """Takes a Structure instance and convert it to the following format:
     [(a, "("), (b, "("), ..., (-a, ")")]
  """
  if isinstance(x, Structure):
    l = builtin_functions.flip([[x.domains, x.dotparens], len(x.domains)])
    l = map(lambda x: (x[0], x[1]), l)
    return flatten(l)
  if type(x) == list:
    def concat(l):
      res = []
      for x in l:
        res += x
      return res
    return concat(map(flatten, x))
  else:
    x, y = x
    if y != "~":
      return [(x, y)]
    else:
      return flatten(x)

def strip_consecutive_strandbreaks(l):
  """Make sure that there are no empty sequences.  """
  #print 'l', l
  flag = 1
  res = []
  for (x, y) in l:
    if x == "+" and flag == 1: continue
    if x == "+":
      flag = 1
    else:
      flag = 0
    res.append((x, y))
  (x, y) = res[-1]
  if x == "+": res = res[:-1]
  return res

def rotate(complex):
  """Take a complex in form of a list of Strands and their secondary structure,
  pop() the first element and append it at the end. 
  
  Example:
    [[S1, S2, S3],['(', '+', '.', ')', '+', '.']]
    [[S2, S3, S1],['.', '(', '+', '.', '+', ')']]

  :param complex: a list of lists of Strand Objects and their structure 

  :return: a list of lists [[Strand Objects],[dotparens]]
  """
  def find(l, key):
    for i in range(len(l)):
      if l[i] == key:
        return i
    return None

  def hardcopyList(l):
    if type(l) != list:
      return l
    return map(hardcopyList, l)

  complex = hardcopyList(complex)

  if "+" not in complex[1]:
    return complex        
  else:
    dom = complex[0][1:] + complex[0][:1]

    # change parentheses appropriately
    p = find(complex[1], "+")
    dpr = complex[1]
    stack = []
    for i in range(p):
      if dpr[i] == "(": stack.append(i)
      elif dpr[i] == ")": stack.pop()
    for i in stack:
      dpr[i] = ")"
    stack = []
    for i in reversed(range(p + 1, len(dpr))):
      if dpr[i] == ")": stack.append(i)
      elif dpr[i] == "(": stack.pop()
    for i in stack:
      dpr[i] = "("

    dpr = dpr[p + 1:] + ["+"] + dpr[:p]
    return [dom, dpr]

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

def interpret(ts_parsed, crn_parsed, fs_list, 
    name='ts_author', sdlen=6, ldlen=15, verbose=True):
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
  ts_env.translate_formal_species(fs_list)
  # translate the crn using the main() function 
  ts_env.translate_reactions(crn_parsed)

  # get the results in form of a dictionary d={fs:Object}
  fs_result = ts_env.formal_species_dict
  # get the final solution object
  cs_result = ts_env.constant_species_solution

  return post_process(fs_result, cs_result)

def post_process(fs_result, solution):
  """Convert output to DNA Objects format.

  Args:
    fs_result (Dict{'FS':Species()}) : Compiled formal species
    solution (Solution()): Compiled Solution of the DSD S
  """

  # flatten outputs
  for k in fs_result.keys():
    fs_result[k] = \
        strip_consecutive_strandbreaks(flatten(fs_result[k]))

  solution_as_list = list()
  for m in solution.molecules:
    # Hack to enable empty solution objects.. 
    # not exactly sure if we want this
    if flatten(m) == [] : continue
    solution_as_list.append(strip_consecutive_strandbreaks(flatten(m)))

  fs_list = []
  cs_list = []

  # list of (DNAObject.Domain) domains in the system, modified by get_domain
  domains = []
  def get_domain(x):
    """Convert Domain format from Nuskell into DNAObjects """
    domain_name = str(x)
    starred = False
    if domain_name[-1] == "*":
      domain_name = domain_name[:-1]
      starred = True

    # Check if we have seen the domain previously
    for d in domains:
      if d.name == domain_name:
        return d.complement if starred else d

    # Otherwise initialize a new domain object
    constraint = 'N' * x.length
    new_dom = DNAObjects.Domain(name = domain_name, constraints = constraint)
    domains.append(new_dom)
    return new_dom.complement if starred else new_dom

  # list of (DNAObject.Strand) strands in the system, modified by get_strand
  strands = []
  def get_strand(strand):
    for x in strands:
      if x.domains == strand:
        return x
    new_strand = DNAObjects.Strand(domains = strand)
    strands.append(new_strand)
    return new_strand

  wildcard_domain = DNAObjects.Domain(name = "?", constraints = 'N' * 15)

  # convert formal species
  for fs_name in fs_result:
    # add a strand break at the end for convenience
    complex_old_format = fs_result[fs_name] + [("+", "+")]
    complex = []
    strand = []
    structure = ""
    for (x, y) in complex_old_format:
      structure += y
      if x == "+":
        strand = get_strand(strand)
        complex.append(strand)
        strand = []
        continue
      if x == "?":
        x = wildcard_domain
      else:
        x = get_domain(x)
      strand.append(x)

    # remove the strand break that was added at the beginning
    structure = structure[:-1]
    fs_list.append(
        DNAObjects.Complex(name = fs_name, 
          strands = complex, structure = structure))

  previous = []
  # convert constant species
  for cs in solution_as_list:
    # add a strand break at the end for convenience
    complex_old_format = cs + [("+", "+")]
    complex = []
    strand = []
    structure = ""
    for (x, y) in complex_old_format:
      structure += y
      if x == "+":
        strand = get_strand(strand)
        complex.append(strand)
        strand = []
        continue
      if x == "?":
        x = wildcard_domain
      else:
        x = get_domain(x)
      strand.append(x)
    # remove the strand break that was added at the beginning
    structure = structure[:-1]

    # has it been already added?
    seen = False
    for y in previous:
      x = [complex, list(structure)]
      p = rotate(x)
      while True:
        if p == y:
          seen = True
        if p == x:
          break
        p = rotate(p)

    if not seen:
      previous.append([complex, list(structure)])
      cs_list.append(
          DNAObjects.Complex(strands = complex, structure = structure))

  return domains, strands, fs_list, cs_list

