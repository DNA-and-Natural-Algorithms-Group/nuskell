#!/usr/bin/env python
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@dna.caltech.edu)
#
# Compiler module.
#

import os
import sys
import argparse

from nuskell.parser import parse_crn_string, parse_ts_file
from nuskell.interpreter.interpreter import interpret
from nuskell.enumeration import enumerate_crn, enumerate_crn_old
from nuskell.verifier.verifier import verify

from nuskell.include.peppercorn.enumerator import get_peppercorn_args

def translate(input_crn, ts_file, pilfile=None, domfile=None, sdlen=6, ldlen=15):
  """A formal chemical reaction network (CRN) is translated into domain level
  representation (DOM) of a DNA strand displacement circuit (DSD).  The
  translation-scheme has to be formulated using the nuskell programming
  language. 

  Args: 
    input_crn: An input string representation of the formal CRN
    is_file: The input file name of a translation scheme
    pilfile: The output file name of a DOM-level cirucit in .pil format
    domfile: The output file name of a DOM-level cirucit in .dom format

  Returns:
    domains: A list of Domain objects (see DNAObjects)
    strands: A list of Strand objects (see DNAObjects)
    formal_species: A list of Complex objects (see DNAObjects)
    constant_species: A list of Complex objects (see DNAObjects)

  """

  ts = parse_ts_file(ts_file)
  (crn, formal_species, const_species) = parse_crn_string(input_crn)

  domains, strands, formal_species, constant_species = \
      interpret(ts, crn, formal_species, const_species, sdlen=sdlen, ldlen=ldlen)

  if pilfile :
    print_as_PIL(pilfile, domains, strands, formal_species, constant_species)
  if domfile :
    print_as_DOM(domfile, domains, strands, formal_species, constant_species)

  return domains, strands, formal_species, constant_species

def print_as_PIL(outfile, domains, strands, formal_species, constant_species):
  """ Print an implementation CRN in the pepper internal language (PIL) file format. 
    
  Args:
    outfile: filename of the PIL file
    domains: array of Domain objects (see DNAObjects)
    strands: array of Strand objects (see DNAObjects)
    formal_species: A list of Complex objects (see DNAObjects)
    constant_species: A list of Complex objects (see DNAObjects)
  """
  F = open(outfile, 'w')
  for x in domains:
    if x.name != "?":
      print >> F, "sequence " + x.name + " = " + "N" * x.length + " : " + str(x.length)
  for x in strands:
    total_length = 0
    print >> F, "strand " + x.name + " =",
    for d in x.domains:
      total_length += d.length
      print >> F, d.name,
    print >> F, ": " + str(total_length)

  print >> F, "# Formal species"
  for x in formal_species:
    print >> F, "structure " + x.name + " =",
    first = True
    for s in x.strands:
      if first: first = False
      else: print >> F, "+",
      print >> F, s.name,
    print >> F, ": "  + x.structure.to_dotparen()

  print >> F, "# Constant species"
  for x in constant_species:
    print >> F, "structure " + x.name + " =",
    first = True
    for s in x.strands:
      if first: first = False
      else: print >> F, "+",
      print >> F, s.name,
    print >> F, ": "  + x.structure.to_dotparen()
  # TODO : figure out a way to deal with the wildcards.
  return

def print_as_DOM(outfile, domains, strands, formal_species, constant_species):
  """ Print an implementation CRN in the domain (DOM) file format. 

  Args:
    outfile: filename of the DOM file
    domains: array of Domain objects (see DNAObjects)
    strands: array of Strand objects (see DNAObjects)
    formal_species: A list of Complex objects (see DNAObjects)
    constant_species: A list of Complex objects (see DNAObjects)
  """
  F = open(outfile, 'w')
  for x in domains:
    if x.name != "?":
      print >> F, "sequence " + x.name + " : " + str(x.length)

  print >> F, "# Formal species"
  for x in formal_species:
    print >> F, x.name + " :"
    first = True
    t = []
    for s in x.strands:
      if first: first = False
      else: print >> F, "+",
      for d in s.domains:
        print >> F, d.name,
        t.append(d)
    print >> F
    i = 0
    while i < len(x.structure.to_dotparen()):
      if x.structure.to_dotparen()[i] == "+":
        print >> F, "+",
        i += 1
      else:
        if len(t) == 0:
          break
        elif t[0].length == 0:
          t = t[1:]
        else:
          print >> F, x.structure.to_dotparen()[i],
          i += t[0].length
          t = t[1:]
    print >> F
              
  print >> F, "# Constant species"
  for x in constant_species:
    print >> F, x.name + " :"
    first = True
    t = []
    for s in x.strands:
      if first: first = False
      else: print >> F, "+",
      for d in s.domains:
        print >> F, d.name,
        t.append(d)
    print >> F
    i = 0
    while i < len(x.structure.to_dotparen()):
      if x.structure.to_dotparen()[i] == "+":
        print >> F, "+",
        i += 1
      else:
        if len(t) == 0:
          break
        elif t[0].length == 0:
          t = t[1:]
        else:
          print >> F, x.structure.to_dotparen()[i],
          i += t[0].length
          t = t[1:]
    print >> F
  return

def get_nuskell_args(parser) :
  """ A collection of arguments for Nuskell 
  nuskell reads and processes chemical reaction networks (CRNs). Three
  different modes are supported to process this CRN:
    
    Translate a formal CRN to a domain-level (DOM) implementation CRN 
      - requires a translation scheme file

    Compare a formal CRN to an implementation CRN using bisimulation
      - requires an implementation CRN file 
      - optionally reads an interpretation CRN 

    .. Complain if both a translation scheme and an iCRN are specified.
  """
  # Enter the translation mode of Nuskell
  parser.add_argument("--ts", required=True, action = 'store',
      help="Specify path to the translation scheme")

  # Enter the verification-only mode of Nuskell (Robert's mode)
  parser.add_argument("--compare", nargs="+", action = 'store',
      help="Specify path to an implementation CRN and (optionally) \
      also to an interpretation CRN.")

  # Choose a verification method.
  parser.add_argument("--verify", default='', action = 'store',
      help="Specify name a verification method: \
          (bisimulation, pathway, integrated, bisim-loopsearch, bisim-wholegraph)") 

  # Convenience options 
  parser.add_argument("-o", "--output", default='', action = 'store',
      help="Specify name of output file")
  parser.add_argument("--dom_short", type=int, default=6,
      help="Length of short domains when using the short() built-in function")
  parser.add_argument("--dom_long", type=int, default=15,
      help="Length of long domains when using the long() built-in function")
  return parser

def main() :
  """ The Nuskell compiler.  
  
  Commandline-parameters are collected in order to 
    - translate formal CRNs into implementation CRNs. 
    - verify the equivalence between fCRN and iCRN.

  Output:
    - Domain-level DSD circuits printed into .pil and/or .dom files
    - Verbose information during verification
    - Result of verification
    
  """
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser = get_nuskell_args(parser)
  parser = get_peppercorn_args(parser)
  args = parser.parse_args()

  # Parse the input CRN 
  input_crn = sys.stdin.readlines()
  input_crn = "".join(input_crn)

  if args.output == '' :
    args.output = 'implementation'

  pilfile = args.output + '.pil'
  domfile = args.output + '.dom'
  d, s, f, c = translate(input_crn, args.ts,
          pilfile=pilfile, 
          domfile=domfile, 
          sdlen = args.dom_short,
          ldlen = args.dom_long)

  # Determine if we need the enumerator, at this point, we only need it for
  # verification, but in the future, we need it also for simulation,
  # optimization, etc.
  if True and args.verify :
    print "Translation done. Enumerating species of the implementation CRN:"
    
    # TODO: need to check *ignoring* of history domains in the pil-file before 
    # switching to the new method. 
    old = True
    if old :
      enum_crn, init_cplxs, enum_cplxs = enumerate_crn_old(args, domfile)
    else :
      enum_crn, init_cplxs, enum_cplxs = enumerate_crn(args, pilfile)

  if args.verify :
    print "Verification using:", args.verify
    # --pathway --bisimulation --integrated

    v = verify(input_crn, enum_crn, init_cplxs, enum_cplxs, 
        method = args.verify, verbose = True)

    if v:
      print "verify: compilation was correct."
    else:
      print "verify: compilation was incorrect." 

if __name__ == '__main__':
  main()

