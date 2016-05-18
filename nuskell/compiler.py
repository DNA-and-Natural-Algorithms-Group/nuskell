#!/usr/bin/env python
#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Compiler module.
#

import os
import sys

import nuskell.parser.ts_parser as ts_parser
import nuskell.parser.crn_parser as crn_parser

from nuskell.interpreter.interpreter import interpret
from nuskell.verifier.verifier import verify

def compile(ts_file, input_crn, pilfile=None, domfile=None, 
    sdlen=6, ldlen=15):
  """ A formal chemical reaction network (CRN) is compiled into a DNA strand
  displacement circuit in domain level representation (DOM). The translation
  scheme must be formulated in the nuskell programming language. 

  :param ts_file: The input file name of a translation scheme
  :param input_crn: An input string representation of the formal CRN
  :param pilfile: The output file name of a DOM-level cirucit in .pil format
  :param pilfile: The output file name of a DOM-level cirucit in .dom format
  """

  ts = ts_parser.parse_file(ts_file)
  (crn, formal_species, const_species) = crn_parser.parse_string(input_crn)

  domains, strands, formal_species, constant_species = \
      interpret(ts, crn, formal_species, sdlen=sdlen, ldlen=ldlen)

  if pilfile :
    print_as_PIL(pilfile, domains, strands, formal_species, constant_species)
  if domfile :
    print_as_DOM(domfile, domains, strands, formal_species, constant_species)
  return 

def print_as_PIL(outfile, domains, strands, formal_species, constant_species):
  """ Print a compiled circuit in the pepper internal language (PIL) file format. 

  :param outfile: filename of the PIL file
  :param domains: array of *domain* objects (see DNAObjects)
  :param strands: array of *strand* objects (see DNAObjects)
  :param formal_species: array of *complex* objects (see DNAObjects)
  :param constant_species: array of *complex* objects (see DNAObjects)
  """
  F = open(outfile, 'w')
  for x in domains:
    if x.name != "?":
      print >> F, "sequence " + x.name + " = " + "N" * x.length + " : " + str(x.length)
  for x in strands:
    total_length = 0
    print >> F, "strand " + x.name + " =",
    for d in x.domain_list:
      total_length += d.length
      print >> F, d.name,
    print >> F, ": " + str(total_length)

  print >> F, "# Formal species"
  for x in formal_species:
    print >> F, "structure " + x.name + " =",
    first = True
    for s in x.strand_list:
      if first: first = False
      else: print >> F, "+",
      print >> F, s.name,
    print >> F, ": "  + x.structure

  print >> F, "# Constant species"
  for x in constant_species:
    print >> F, "structure " + x.name + " =",
    first = True
    for s in x.strand_list:
      if first: first = False
      else: print >> F, "+",
      print >> F, s.name,
    print >> F, ": "  + x.structure
  # TODO : figure out a way to deal with the wildcards.
  return

def print_as_DOM(outfile, domains, strands, formal_species, constant_species):
  """ Print a compiled circuit in the domain (DOM) file format. 

  :param outfile: filename of the DOM file
  :param domains: array of *domain* objects (see DNAObjects)
  :param strands: array of *strand* objects (see DNAObjects)
  :param formal_species: array of *complex* objects (see DNAObjects)
  :param constant_species: array of *complex* objects (see DNAObjects)
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
    for s in x.strand_list:
      if first: first = False
      else: print >> F, "+",
      for d in s.domain_list:
        print >> F, d.name,
        t.append(d)
    print >> F
    i = 0
    while i < len(x.structure):
      if x.structure[i] == "+":
        print >> F, "+",
        i += 1
      else:
        if len(t) == 0:
          break
        elif t[0].length == 0:
          t = t[1:]
        else:
          print >> F, x.structure[i],
          i += t[0].length
          t = t[1:]
    print >> F
              
  print >> F, "# Constant species"
  for x in constant_species:
    print >> F, x.name + " :"
    first = True
    t = []
    for s in x.strand_list:
      if first: first = False
      else: print >> F, "+",
      for d in s.domain_list:
        print >> F, d.name,
        t.append(d)
    print >> F
    i = 0
    while i < len(x.structure):
      if x.structure[i] == "+":
        print >> F, "+",
        i += 1
      else:
        if len(t) == 0:
          break
        elif t[0].length == 0:
          t = t[1:]
        else:
          print >> F, x.structure[i],
          i += t[0].length
          t = t[1:]
    print >> F
  return

def main() :
  """ Main interface to the nuskell compiler. This function uses arparse to
  collect non-standard system parameters, compiles the formal CRN to a
  DOM-level circuit and prints the results as PIL or DOM files. The function
  also allows to verify if a given DOM representation implements the input CRN
  (see nuskell.verifier).
  """
  import sys
  import argparse
  def get_nuskell_args() :
    """ A collection of arguments for nuskell 
    nuskell reads and processes chemical reaction networks (CRNs). Three
    different modes are supported to process this CRN:
      
      Compile a formal CRN to a domain-level (DOM) implementation CRN 
        - requires a translation scheme file

      Compare a formal CRN to an implementation CRN using bisimulation
        - requires an implementation CRN file 
        - optionally reads an interpretation CRN 

      Verify that the formal CRN is implemented by an implementation CRN using
      pathway decomposition
        - requires a translation scheme file

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--compile", action = 'store',
        help="Specify path to the translation scheme")
  
    parser.add_argument("--compare", action = 'store',
        help="Specify path to the translation scheme")

    parser.add_argument("--ts", required=True, action = 'store',
        help="Specify path to the translation scheme")
    parser.add_argument("-o", "--output", default='', action = 'store',
        help="Specify name of output file")

    parser.add_argument("--dom_short", type=int, default=6,
        help="Length of short domains when using the short() built-in function")
    parser.add_argument("--dom_long", type=int, default=15,
        help="Length of long domains when using the long() built-in function")

    # could make it an array of methods: 
    #   [bisimulation, pathway-equivalence]
    parser.add_argument("--verify", default='', action = 'store',
        help="Specify name a verification method: \
            (standard, wolfe, ...)")
    #parser.add_argument("--verfy_condense", action = 'store_true',
    #    help="Verify the equivalence between detailed semantics and " +
    #    "condensed semantics for DOMFILE or for CRNFILE compiled " + 
    #    "using TSFILE."
  
    return parser.parse_args()
  args = get_nuskell_args()

  # Parse the input CRN 
  input_crn = sys.stdin.readlines()
  input_crn = "".join(input_crn)

  if args.output == '' :
    args.output = 'nuskell_file'

  pilfile = args.output + '.pil'
  domfile = args.output + '.dom'
  compile(args.ts, input_crn, 
      pilfile=pilfile, 
      domfile=domfile, 
      sdlen = args.dom_short,
      ldlen = args.dom_long)
 
  if args.verify :
    print "Compilation done. enumerating pathways ... "

    # if "--pathway" in options: method = "--pathway"
    # if "--bisimulation" in options: method = "--bisimulation"
    # if "--integrated" in options: method = "--integrated"
    # if "--interactive" in options: interactive = True

    print "Enumeration done. Verification using:", args.verify
    v = verify(input_crn, domfile, method = args.verify, verbose = True) # condense = True, interactive = False

    if v:
      print "verify: compilation was correct."
    else:
      print "verify: compilation was incorrect." 

if __name__ == '__main__':
  main()

