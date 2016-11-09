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
from nuskell.parser import split_reversible_reactions
from nuskell.parser import combine_reversible_reactions

from nuskell.interpreter import interpret
from nuskell.enumeration import TestTubePeppercornIO
from nuskell.verifier import preprocess, verify
from nuskell.objects import TestTube, TestTubeIO

from nuskell.include.peppercorn.enumerator import get_peppercorn_args

def translate(input_crn, ts_file, pilfile=None, domfile=None, verbose = False):
  """A formal chemical reaction network (CRN) is translated into domain level
  representation (DOM) of a DNA strand displacement circuit (DSD).  The
  translation-scheme has to be formulated using the nuskell programming
  language. 

  Args: 
    input_crn: An input string representation of the formal CRN
    ts_file: The input file name of a translation scheme
    pilfile: The output file name of a domain-level cirucit in .pil format

  Returns:
    solution: A TestTube object 
    constant_soluiton: A TestTube object that contains only constant species
  """

  ts = parse_ts_file(ts_file)
  (crn, formal_species, const_species) = parse_crn_string(input_crn)

  solution, constant_solution = interpret(ts, crn, 
      formal_species, const_species)

  if pilfile :
    TestTubeIO(solution).write_pilfile(pilfile)
  if domfile :
    TestTubeIO(solution).write_domfile(pilfile)

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
  parser.add_argument("--ts", action = 'store',
      help="Specify path to the translation scheme")

  parser.add_argument("--pilfile", action = 'store',
      help="Specify path to a *.pil file.")

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

  # ~~~~~~~~~~~~~~~~~~~~~~~~
  # Prepare CRN and TestTube
  # ~~~~~~~~~~~~~~~~~~~~~~~~
  if args.ts : # Translate CRN using a translation scheme
    print "\nTranslating..."
    solution, constant_solution = translate(input_crn, args.ts, 
        verbose = (args.verbose > 1)) 
    pilfile = args.output + '.pil'
    domfile = args.output + '.dom'
    TestTubeIO(solution).write_pilfile(pilfile)
    TestTubeIO(solution).write_domfile(domfile)
    print "Wrote file:", pilfile, domfile
  elif args.pilfile : # Parse implementation species from a PIL file
    print "Parsing PIL file..."
    solution = TestTube()
    raise NotImplementedError
    solution.load_pilfile(args.pilfile)
  else :
    raise NotImplementedError("Automated choice of translation scheme is not implemented")
    solution = None

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Prepare enumerated CRN and TestTube
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if args.verify or args.simulate or args.enumerate :
    print "\nEnumerating reaction pathways..."

    peppercorn = TestTubePeppercornIO(solution, args)
    peppercorn.enumerate()
    enum_solution = peppercorn.testtube
    enum_crn = peppercorn.condense_reactions
    #enum_crn = peppercorn.all_reactions

    if args.verbose :
      rev_crn = combine_reversible_reactions(enum_crn)
      print "Enumerated CRN:"
      for rxn in rev_crn :
        if rxn[2] == 'reversible' :
          print ' + '.join(rxn[0]), '<=>', ' + '.join(rxn[1])
        else :
          print ' + '.join(rxn[0]), '->', ' + '.join(rxn[1])

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Verify equivalence of CRNs
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  if args.verify :
    print "\nVerification preprocessing..."
    (fcrn, fs, cs) = parse_crn_string(input_crn) 
    fcrn = split_reversible_reactions(fcrn)

    if True :
      # Reduce the enumerated CRN and find an interpretation
      icrn, interpret = preprocess(fcrn, enum_crn, fs,
          solution, enum_solution, verbose=(args.verbose>1))
    else :
      # Reduce the enumerated CRN without finding an interpretation
      import nuskell.verifier.verifier
      cs = filter(lambda x: x not in fs, solution.complexes)
      icrn = nuskell.verifier.verifier.removeSpecies(enum_crn, cs)
      icrn = nuskell.verifier.verifier.removeDuplicates(icrn)
      interpret = None

    if args.verbose :
      print "Implementation CRN:"
      rev_crn = combine_reversible_reactions(icrn)
      for rxn in rev_crn :
        if rxn[2] == 'reversible' :
          print ' + '.join(rxn[0]), '<=>', ' + '.join(rxn[1])
        else :
          print ' + '.join(rxn[0]), '->', ' + '.join(rxn[1])

    print "\nVerification using:", args.verify
    if args.verbose :
      print "Partial interpretation:"
      for impl, formal in sorted(interpret.items()) :
        print "  {} => {}".format(impl, ', '.join([x for x in formal.elements()]))

    v, i = verify(fcrn, icrn, fs, interpret=interpret, 
        method=args.verify, verbose=(args.verbose>1))

    if i and args.verbose :
      if not v : i = i[0]
      print "Returned interpretation:"
      for impl, formal in sorted(i.items()) :
        print "  {} => {}".format(impl, ', '.join([x for x in formal.elements()]))


    if v:
      print "\nCompilation was correct."
    else:
      print "\nCompilation was incorrect." 

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Simulate CRNs in a TestTube
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if args.simulate :
    raise NotImplementedError

if __name__ == '__main__':
  main()

