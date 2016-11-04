# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@caltech.edu)
#
# Preprocessing and interface to peppercorn enumerator
#

from nuskell.parser import parse_dom_file
from nuskell.objects import Domain, Complex, TestTube

import nuskell.include.peppercorn.utils as peputils
from nuskell.include.peppercorn.enumerator import Enumerator
import nuskell.include.peppercorn.reactions as reactions
from nuskell.include.peppercorn.condense import condense_resting_states

# Once peppercorn is its own package, we import it as depencency.
# try :
#   import peppercorn
# except ImportError:
#   sys.exit("""nuskell depends on the peppercorn package
#   -- download at http://dna.caltech.edu/peppercorn """)

def peppercorn_enumerate(args, solution, condense = True, verbose=False):
  """Nuskell interface to the enumerator. 

  Args:
    args (Argparse Object): Arguments for the peppercorn enumerator
    solution (TestTube()):  A set of complexes for enumeration
  
  Returns:
    enum_crn ([[[re],[pr]],..]: A CRN of enumerated species
    enum_solution (TestTube()): Enumerated complexes
  """
  enum = initialize_enumerator(solution)
  set_enumargs(enum, args)

  enum.enumerate()

  if condense :
    enum_crn = get_enum_crn_condense(enum, verbose)
  else :
    enum_crn = get_enum_crn(enum, verbose)
  enum_solution = get_enum_solution(enum, solution, verbose)

  return enum_crn, enum_solution

def initialize_enumerator(solution) :
  """Initialize the peppercorn enumerator object.

  Args:
    solution (nuskell.objects.TestTube()): A set of complexes

  Returns:
    Enumerator (peppercorn.enumerator.Enumerator())
  """

  # Translate to peppercorn domains
  domains = {}
  for n, d in solution.domains.items() :
    if n[-1] == '*' :
      new_dom = peputils.Domain(n[:-1], d.length, sequence=''.join(d.sequence), is_complement=True)
    else :
      new_dom = peputils.Domain(n, d.length, sequence=''.join(d.sequence))
    domains[n] = new_dom
  #print domains.values()

  # Translate to peppercorn strands
  strands = {}
  dom_to_strand = {}
  for n, s in solution.strands.items() :
    dom_to_strand[tuple(map(str,s))] = n
    doms = []
    for d in map(str, s) :
      doms.append(domains[d])
    strands[n] = peputils.Strand(n, doms)
  #print strands.values()

  # Translate to peppercorn complexes
  complexes = {}
  for n, c in solution.complexes.items():
    cplx_strands = []
    for s in c.lol_sequence :
      ns = dom_to_strand[tuple(map(str,s))]
      cplx_strands.append(strands[ns])
    complex_structure = peputils.parse_dot_paren(''.join(c.structure))
    complex = peputils.Complex(n, cplx_strands, complex_structure)
    complex.check_structure()
    complexes[n] = complex		
  #print complexes.values()

  domains = domains.values()
  strands = strands.values()
  complexes = complexes.values()

  return Enumerator(domains, strands, complexes)

def enum_cplx_rename(x) :
  """ A function to rename enumerator species names to a format 
  which is compatible with nuskell.objects
  """
  x = 'e'+x if x[0].isdigit() else x
  return x +'_' if x[-1].isdigit() else x

def enum_rs_rename(x) :
  """ A function to rename enumerator species names to a format 
  which is compatible with nuskell.objects
  """
  x = 'r'+x if x[0].isdigit() else x
  return x +'_' if x[-1].isdigit() else x


def enum_domain_rename(x) :
  # NOTE: This function is obsolete, because domains are initialized using
  # the original solution object. However, it translates domain names to be
  # compatible with nuskell.objects
  x = 'e'+x if x[0].isdigit() else x
  if x[-1]=='*':
    return x[:-1] + 'p' + '*'
  elif x[-1].isdigit() :
    return x + 'p'
  else :
    return x

def get_enum_crn_condense(enum, verbose=False):
  # Extract condensed reactions
  condense_options={}
  condensed = condense_resting_states(enum, **condense_options)
  reactions = condensed['reactions']

  enum_crn = []
  for r in sorted(reactions):
    #print map(lambda x: map(str, x.complexes), r.reactants), '->', 
    #print map(lambda x: map(str, x.complexes), r.products)
    #print r.kernel_string()

    react = []
    for rs in r.reactants:
      if len(rs.complexes) == 1: 
        react.append(enum_cplx_rename(str(rs.complexes[0])))
      else :
        react.append(enum_rs_rename(rs.name))

    prod = []
    for rs in r.products:
      if len(rs.complexes) == 1: 
        prod.append(enum_cplx_rename(str(rs.complexes[0])))
      else :
        prod.append(enum_rs_rename(rs.name))

    enum_crn.append([react, prod])
  return enum_crn

def get_enum_crn(enum, verbose=False):
  reactions = enum.reactions
  enum_crn = []
  for r in sorted(reactions):
    react = map(enum_cplx_rename, map(str, r.reactants))
    prod = map(enum_cplx_rename, map(str, r.products))
    #print r.kernel_string()
    #print [react, prod]
    enum_crn.append([react, prod])
  return enum_crn

def get_enum_solution(enum, solution, verbose=False):
  # Extraction of complexes in the enumerated CRN:
  enum_cplxs = TestTube()
  enum_cplxs += solution
  for rs in enum.resting_states:
    enum_cplxs.load_enum_cplxs(rs.complexes, enum_cplx_rename)

  return enum_cplxs

def set_enumargs(enum, args):
  """Transfer options to Enumerator Object. 
  
  Do NOT set Nuskell-defaults here. Defaults are set with the argparse object
  of peppercorn: nuskell/include/peppercorn/enumerator.py 
  
  """

  #enum.MAX_REACTION_COUNT = 5000000
  #enum.MAX_COMPLEX_COUNT  = 100000
  #enum.MAX_COMPLEX_SIZE   = 1000
  #enum.k_fast = 2.0
  #enum.REJECT_REMOTE = True
  if args.verbose is not None:
    import logging
    logger = logging.getLogger()
    if args.verbose == 1:
      logger.setLevel(logging.INFO)
    elif args.verbose == 2:
      logger.setLevel(logging.DEBUG)
    elif args.verbose >= 3:
      logger.setLevel(logging.NOTSET)

  if args.k_slow is not None:
    enum.k_slow = args.k_slow
  if args.k_fast is not None:
    enum.k_fast = args.k_fast

  if args.MAX_REACTION_COUNT is not None:
    enum.MAX_REACTION_COUNT = args.MAX_REACTION_COUNT

  if args.MAX_COMPLEX_COUNT is not None:
    enum.MAX_COMPLEX_COUNT = args.MAX_COMPLEX_COUNT

  if args.MAX_COMPLEX_SIZE is not None:
    enum.MAX_COMPLEX_SIZE = args.MAX_COMPLEX_SIZE

  if args.RELEASE_CUTOFF is not None:
    enum.RELEASE_CUTOFF = args.RELEASE_CUTOFF

  if args.RELEASE_CUTOFF_1_1 is not None:
    enum.RELEASE_CUTOFF_1_1 = args.RELEASE_CUTOFF_1_1

  if args.RELEASE_CUTOFF_1_N is not None:
    enum.RELEASE_CUTOFF_1_N = args.RELEASE_CUTOFF_1_N

  if args.REJECT_REMOTE is not None:
    enum.REJECT_REMOTE = args.REJECT_REMOTE

  if args.UNZIP is not None:
    enum.UNZIP = args.UNZIP

  if args.LEGACY_UNZIP is not None:
    enum.LEGACY_UNZIP = args.LEGACY_UNZIP

  # TODO: what is this?
  # enum.DFS = not args.bfs

  # Modify enumeration events based on command line options.
  if args.ignore_branch_3way:
    if reactions.branch_3way in enum.FAST_REACTIONS:
      enum.FAST_REACTIONS.remove(reactions.branch_3way)

  if args.ignore_branch_4way:
    if reactions.branch_4way in enum.FAST_REACTIONS:
      enum.FAST_REACTIONS.remove(reactions.branch_4way)

  return

