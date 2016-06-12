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

from nuskell.parser import parse_crn_string, parse_ts_file
from nuskell.parser import parse_dom_file

from nuskell.interpreter.interpreter import interpret
from nuskell.verifier.verifier import verify

def compile(input_crn, ts_file, pilfile=None, domfile=None, sdlen=6, ldlen=15):
  """A formal chemical reaction network (CRN) is compiled into domain level
  representation (DOM) of a DNA strand displacement circuit (DSD).  The
  translation-scheme has to be formulated using the nuskell programming
  language. 

  :param input_crn: An input string representation of the formal CRN
  :param ts_file: The input file name of a translation scheme
  :param pilfile: The output file name of a DOM-level cirucit in .pil format
  :param pilfile: The output file name of a DOM-level cirucit in .dom format
  """

  ts = parse_ts_file(ts_file)
  (crn, formal_species, const_species) = parse_crn_string(input_crn)

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

def main() :
  """Standard interface to the nuskell compiler.  Commandline-parameters are
  collected in order to compile and/or verify CRNs using DNA strand
  displacement translation-schemes.  Domain-level DSD circuits are printed in
  the .pil or .dom fileformat, verbose information, as well as results for
  verification are printed to STDOUT.
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
  compile(input_crn, args.ts,
      pilfile=pilfile, 
      domfile=domfile, 
      sdlen = args.dom_short,
      ldlen = args.dom_long)


 
  if args.verify :
    print "Compilation done. enumerating pathways ... "
    import peppercorn.input as pepin
    import peppercorn.output as pepout

    #TODO: call enumerator_argparse, pass parser_object on
    enum = pepin.input_pil(pilfile)
    enum.MAX_COMPLEX_COUNT  = 10000
    enum.MAX_REACTION_COUNT = 50000
    enum.MAX_COMPLEX_SIZE   = 100
    enum.enumerate()

    # Write the output of peppercorn into the enumfile
    enumfile = args.output + '.enum'
    pepout.output_crn(enum, enumfile, output_condensed = True)

    # Post-process enumerator results to extract the condensed crn
    enum_crn = []
    with open(enumfile, 'r') as enu :
      for line in enu :
        react = line.split('->')
        react[0] = sorted([x.strip() for x in react[0].split('+')])
        react[1] = sorted([x.strip() for x in react[1].split('+')])
        enum_crn.append(react)

    slow_cplxs = []
    for rs in enum.resting_states:
      name = str(rs)
      for cx in rs.complexes:
        cxs = []
        #TODO: this should be easier in the peppercorn interface
        for sd in cx.strands:
          cxs.append('+')
          for ds in map(str, sd.domains):
            if ds[-1] == '*':
              cxs.append([ds[:-1], '*'])
            else:
              cxs.append([ds])
        # remove the first '+' again
        if len(cxs)>0: cxs = cxs[1:]
        cx = [name, cxs, list(cx.dot_paren_string())]
        slow_cplxs.append(cx)

    dom = parse_dom_file(domfile)
    cplxs = dom[1] if len(dom)==2 else dom[0] # else: no sequence information

    #print 'ec', enum_crn
    #print 'sc', slow_cplxs
    #print 'ic', input_crn
    #print 'cp', cplxs

    print "Enumeration done. Verification using:", args.verify
    # --pathway --bisimulation --integrated
    v = verify(input_crn, enum_crn, cplxs, slow_cplxs, 
        method = args.verify, verbose = True)

    if v:
      print "verify: compilation was correct."
    else:
      print "verify: compilation was incorrect." 

if __name__ == '__main__':
  main()

