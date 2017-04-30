#!/usr/bin/env python
#
# Written by Stefan Badelt (badelt@dna.caltech.edu)
#

import os
import sys
import argparse
from collections import Counter

from nuskell import printCRN
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
  parser.add_argument("--verify", nargs='+', default='', action = 'store', 
      metavar = '<str>', help="Specify verification methods: \
          (bisimulation, pathway, integrated, bisim-loop-search,\
          bisim-depth-first, bisim-whole-graph)") 

  parser.add_argument("--verify-timeout", type=int, default=30, metavar='<int>',
      help="Specify time [seconds] to wait for verification to complete.")

  parser.add_argument("-o", "--output", default='', action = 'store',
      help="Specify name of output file")

  parser.add_argument("-v", "--verbose", default=0, action = 'count',
      help="Verbosity of output")

  parser.add_argument("--independent", action = 'store_true',
      help="Do not assume ad-hoc interpretation from species names.")

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

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Parse and process input CRN 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  input_crn = sys.stdin.readlines()
  input_crn = "".join(input_crn)
  crn1, crn1s, _, _ = parse_crn_string(input_crn) 

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Verify equivalence of CRNs
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  print "\nVerification preprocessing..."
  crn1 = split_reversible_reactions(crn1)

  if args.verbose :
    print "Formal CRN:"
    printCRN(crn1, reversible = True, rates=False)

  crn2, crn2s, _, _ = parse_crn_file(args.compare)
  crn2 = split_reversible_reactions(crn2)

  if args.verbose :
    print "Implementation CRN:"
    printCRN(crn2, reversible = True, rates=False)

  interpret = dict()
  if not args.independent:
    print 'Ad-hoc partial interpretation:'
    formals = set(crn1s)
    for sp in crn2s:
      if sp in formals:
        interpret[sp] = Counter([sp])

    if interpret and args.verbose :
      for impl, formal in sorted(interpret.items()) :
        print "  {} => {}".format(impl, ', '.join([x for x in formal.elements()]))

  print "\nVerification using:", args.verify
  for meth in args.verify :
    v, i = verify(crn1, crn2, crn1s, method=meth, interpret=interpret,
        verbose=(args.verbose>1),
        timeout=args.verify_timeout)

    if i and args.verbose :
      if not v : i = i[0]
      print "Returned interpretation ['{}']:".format(meth)
      for impl, formal in sorted(i.items()) :
        print "  {} => {}".format(impl, ', '.join([x for x in formal.elements()]))

    if v is True:
      print " {}: CRNs are {} equivalent.".format(v, meth)
    elif v is False:
      print " {}: CRNs are not {} equivalent.".format(v, meth)
    elif v is None:
      print " {}: {} verification did not terminate within {} seconds.".format(v, 
          meth, args.verify_timeout)

if __name__ == '__main__':
  main()

