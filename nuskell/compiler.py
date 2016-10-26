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
from nuskell.enumeration import peppercorn_enumerate
from nuskell.verifier.verifier import verify
from nuskell.objects import TestTube
from nuskell.include.peppercorn.enumerator import get_peppercorn_args

def translate(input_crn, ts_file, pilfile=None, verbose = False):
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

  solution, constant_solution = interpret(ts, crn, 
      formal_species, const_species)

  if pilfile :
    solution.write_pilfile(pilfile)

  return solution, constant_solution

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

  # Choose a verification method.
  parser.add_argument("--verify", default='', action = 'store',
      help="Specify name a verification method: \
          (bisimulation, pathway, integrated, bisim-loop-search,\
          bisim-depth-first, bisim-whole-graph)") 

  parser.add_argument("--enumerate", action = 'store_true',
      help="Enumerate the implementation CRN.")

  parser.add_argument("--simulate", action = 'store_true',
      help="Simulate the CRNs.")

  parser.add_argument("-o", "--output", default='', action = 'store',
      help="Specify name of output file")

  ## Enter the verification-only mode of Nuskell 
  #parser.add_argument("--compare", nargs="+", action = 'store',
  #    help="Specify path to an implementation CRN and (optionally) \
  #    also to an interpretation CRN.")
  return parser

def main() :
  """ The Nuskell compiler.  
  
  Commandline-parameters are collected in order to 
    - translate formal CRNs into implementation CRNs. 
    - verify the equivalence between formal CRN and enumerated CRN.

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

  # ~~~~~~~~~~~~~~~~~~~~~~~~
  # Prepare CRN and TestTube
  # ~~~~~~~~~~~~~~~~~~~~~~~~
  print "Translating..."
  if args.ts : # Translate CRN using a translation scheme
    solution, constant_solution = translate(input_crn, args.ts, 
        verbose = args.verbose)
    if pilfile : solution.write_pilfile(pilfile)
  elif args.pilfile : # Parse implementation species from a PIL file
    solution = TestTube()
    raise NotImplementedError
    solution.load_pilfile(args.pilfile)
  print "...done."

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Prepare enumerated CRN and TestTube
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if args.verify or args.enumerate :
    print "Enumerating..."
    #TODO: remove pilfile
    enum_crn, enum_solution = peppercorn_enumerate(args, pilfile, solution, 
        verbose = args.verbose)
    print "...done."

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Verify equivalence of CRNs
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  if args.verify :
    print "Verification using:", args.verify
    v = verify(input_crn, enum_crn, solution, enum_solution, 
        method = args.verify, verbose = args.verbose)
    if v:
      print "compilation was correct."
    else:
      print "compilation was incorrect." 

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Simulate CRNs in a TestTube
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if args.simulate :
    raise NotImplementedError

if __name__ == '__main__':
  main()

