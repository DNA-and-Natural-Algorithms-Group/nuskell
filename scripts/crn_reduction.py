#!/usr/bin/env python
#
# Written by Stefan Badelt (badelt@dna.caltech.edu)
#
# Condense a strongly connected hypergraph to minimal equivalent CRNs
#

import os
import sys
import argparse

from nuskell.parser import parse_crn_string, parse_crn_file
from nuskell.parser import split_reversible_reactions
from nuskell.parser import combine_reversible_reactions
from nuskell.verifier import verify

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
  # Choose a verification method.
  parser.add_argument("--verify", default='bisimulation', action = 'store',
      help="Specify name a verification method: \
          (bisimulation, pathway, integrated, bisim-loop-search,\
          bisim-depth-first, bisim-whole-graph)") 

  parser.add_argument("-o", "--output", default='', action = 'store',
      help="Specify name of output file")

  parser.add_argument("-v", "--verbose", default=0, action = 'count',
      help="Verbosity of output")

  ## Enter the verification-only mode of Nuskell 
  parser.add_argument("--compare", action = 'store', required=True,
      help="Specify path to an implementation CRN.")
        #and (optionally) also to an interpretation CRN.")
  return parser

def main() :
  """ CRN condensation.
  
  Takes a formal and an implementation CRN and verifies the equivalence
  according to pathway decomposition or bisimulation.

  Throws out edges from the graph to see how equivalence notion changes.

  """
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser = get_nuskell_args(parser)
  args = parser.parse_args()

  # Parse the input CRN 
  input_crn = sys.stdin.readlines()
  input_crn = "".join(input_crn)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Verify equivalence of CRNs
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  print "\nVerification preprocessing..."
  (fcrn, fs, cs) = parse_crn_string(input_crn) 
  fcrn = split_reversible_reactions(fcrn)

  if args.verbose :
    rev_crn = combine_reversible_reactions(fcrn)
    print "Formal CRN:"
    for rxn in rev_crn :
      if rxn[2] == 'reversible' :
        print ' + '.join(rxn[0]), '<=>', ' + '.join(rxn[1])
      else :
        print ' + '.join(rxn[0]), '->', ' + '.join(rxn[1])

  (ccrn, _, _) = parse_crn_file(args.compare)
  ccrn = split_reversible_reactions(ccrn)

  if args.verbose :
    rev_crn = combine_reversible_reactions(ccrn)
    print "Implementation CRN:"
    for rxn in rev_crn :
      if rxn[2] == 'reversible' :
        print ' + '.join(rxn[0]), '<=>', ' + '.join(rxn[1])
      else :
        print ' + '.join(rxn[0]), '->', ' + '.join(rxn[1])

  interpret = dict()
  #if True :
  #  # Reduce the enumerated CRN and find an interpretation
  #  icrn, interpret = preprocess(fcrn, enum_crn, fs,
  #      solution, enum_solution, verbose=(args.verbose>1))
  #else :
  #  # Reduce the enumerated CRN without finding an interpretation
  #  import nuskell.verifier.verifier
  #  cs = filter(lambda x: x not in fs, solution.complexes)
  #  icrn = nuskell.verifier.verifier.removeSpecies(enum_crn, cs)
  #  icrn = nuskell.verifier.verifier.removeDuplicates(icrn)

  print "\nVerification using:", args.verify
  if args.verbose :
    print "Partial interpretation:"
    for impl, formal in sorted(interpret.items()) :
      print "  {} => {}".format(impl, ', '.join([x for x in formal.elements()]))

  v, i = verify(fcrn, ccrn, fs, interpret=interpret, 
      method=args.verify, verbose=(args.verbose>1))

  if i and args.verbose :
    if not v : i = i[0]
    print "Returned interpretation:"
    for impl, formal in sorted(i.items()) :
      print "  {} => {}".format(impl, ', '.join([x for x in formal.elements()]))

  if v:
    print "\nCRNs are equivalent."
  else:
    print "\nCRNs are not equivalent." 

if __name__ == '__main__':
  main()

