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

from nuskell.parser import parse_ts_string
import nuskell.include.DNAObjects as DNAObjects

import nuskell.interpreter.environment as tlsenv
from nuskell.interpreter.environment import Structure
from nuskell.interpreter.environment import Species, Reaction

def flatten(x):
  """Takes a Structure instance and convert it to the following format:
     [(a, "("), (b, "("), ..., (-a, ")")]
  """
  if isinstance(x, Structure):
    l = tlsenv._builtin()._flip([[x.domains, x.dotparens], len(x.domains)])
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
  """Make sure that there are no empty sequences. 
  """
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

def tls_code_snippet():
  """ Sample code of the tls language, which is parsed upon initialization of
  the tls environment.

  :return: a code snippet in the tls language format
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

def interpret(tls_parsed, crn_parsed, fs_list, 
    name='tls_author', sdlen=6, ldlen=15, verbose=True):

  # Setup the environment
  tls_env = tlsenv.Environment(name, sdlen=sdlen, ldlen=ldlen)

 # if verbose : tls_env.print_environment()

  # Parse a piece of sample code with utilities
  header = parse_ts_string(tls_code_snippet())

  # add the code to the environment
  tls_env.interpret(header)

  tls_env.interpret(tls_parsed)
  #raise Exception, "stop here"

  # create formal species
  formal_species_objects = map(Species, fs_list)

  ####### include into evironment? ##########

  # compile the formal species
  tls_env._create_binding("__formalspecies__", formal_species_objects)

  tls_env._env, formal_species_result = tls_env.interpret_expr(
      ["trailer",
        ["id", "map"], 
        ["apply", 
          ["id", "formal"], 
          ["id", "__formalspecies__"]]])

  # final tranlsation of formal species
  formal_species_dict = {}
  for i in range(len(fs_list)):
    #print fs_list, i, formal_species_result[i]
    formal_species_dict[fs_list[i]] = formal_species_result[i]

  # compile reactions -- scary!
  crn_remap = map(lambda x: [x[0]] + map(
    lambda y: map(lambda z: formal_species_dict[z], y), x[1:]), crn_parsed)
  crn_object = map(lambda x: Reaction(x[1], x[2], x[0] == "reversible"), crn_remap)

  tls_env._create_binding("__crn__", crn_object)
  tls_env._env, solution = tls_env.interpret_expr( 
      ["trailer", 
        ["id", "main"], 
        ["apply", 
          ["id", "__crn__"]]])

  # flatten outputs
  for k in formal_species_dict.keys():
    #print formal_species_dict[k]
    #print flatten(formal_species_dict[k])
    formal_species_dict[k] = \
        strip_consecutive_strandbreaks(flatten(formal_species_dict[k]))

  solution_as_list = list()
  for m in solution.molecules:
    solution_as_list.append(strip_consecutive_strandbreaks(flatten(m)))

  #print 'fsd',formal_species_dict 
  #for k in formal_species_dict:
  #  print k, formal_species_dict[k]

  #print 'sal', solution_as_list
  #for s in solution_as_list:
  #  print 's', s


  #################################
  # convert output to DNA Objects #
  #################################
  domains = []
  strands = []
  fs_list = []
  constant_species = []

  def get_domain(x):
      domain_name = str(x)
      starred = False
      if domain_name[-1] == "*":
          domain_name = domain_name[:-1]
          starred = True
      for d in domains:
          if d.name == domain_name:
              if starred: return d.C
              return d
      new_dom = DNAObjects.Domain(name = domain_name, length = x.length)
      domains.append(new_dom)
      if starred: new_dom = new_dom.C
      return new_dom

  def get_strand(strand):
      for x in strands:
          if x.domain_list == strand:
              return x
      new_strand = DNAObjects.Strand(domains = strand)
      strands.append(new_strand)
      return new_strand

  wildcard_domain = DNAObjects.Domain(name = "?", length = 1)

  # convert formal species
  for fs_name in formal_species_dict:
      # add a strand break at the end for convenience
      complex_old_format = formal_species_dict[fs_name] + [("+", "+")]
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
      fs_list.append( \
          DNAObjects.Complex(name = fs_name, \
                             strands = complex, \
                             structure = structure))

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
      flag = False
      for y in previous:
          x = [complex, list(structure)]
          p = rotate(x)
          while True:
              if p == y:
                  flag = True
              if p == x:
                  break
              p = rotate(p)
      if not flag:
          previous.append([complex, list(structure)])
          constant_species.append( \
              DNAObjects.Complex(strands = complex, \
                                 structure = structure))

  return domains, strands, fs_list, constant_species


  # print properties of the current environment
  # print tls_env.report_env()
  # print tls_env.report_tls()
  # print tls_env.report_crn()
  # print tls_env.report_dom()

  # do the translation
  # tls_env.translate(crn)

