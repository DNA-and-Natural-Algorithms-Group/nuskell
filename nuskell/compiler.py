#!/usr/bin/env python
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@dna.caltech.edu)
#
# Compiler module.
#

""" Wrapper functions used by the nuskell compiler script. """

import os
import pkg_resources

from nuskell.parser import parse_crn_string, parse_ts_file
from nuskell.parser import split_reversible_reactions
from nuskell.parser import combine_reversible_reactions

from nuskell.interpreter import interpret
from nuskell.enumeration import TestTubePeppercornIO
from nuskell.objects import TestTube, TestTubeIO

class InvalidSchemeError(Exception):
  """Raise Error: Cannot find scheme."""

  def __init__(self, ts_file, builtin=None):
    self.message = "Cannot find translation scheme: {}\n".format(ts_file)

    if builtin:
      self.message += "You may want to use one of the built-in schemes instead:\n"
      for s in os.listdir(builtin) :
        self.message += " * {}\n".format(s) 
 
    super(InvalidSchemeError, self).__init__(self.message) 

def simulateTT():
  raise NotImplementedError

def printCRN(crn, reversible=True):
  """A wrapper function to print CRNs. """
  if reversible :
    pcrn = combine_reversible_reactions(crn)
  else :
    pcrn = split_reversible_reactions(crn)

  for rxn in pcrn :
    if len(rxn) > 2 and rxn[2] == 'reversible' :
      print ' + '.join(rxn[0]), '<=>', ' + '.join(rxn[1])
    else :
      print ' + '.join(rxn[0]), '->', ' + '.join(rxn[1])

def enumerateTT(testtube, condensed=True, args=None):
  """A wrapper function to enumerate species in a nuskell.TestTube object. """
  peppercorn = TestTubePeppercornIO(testtube, args)
  peppercorn.enumerate()
  enum_testtube = peppercorn.testtube

  if condensed :
    enum_crn = peppercorn.condense_reactions
  else :
    enum_crn = peppercorn.all_reactions

  return enum_testtube, enum_crn

def translate(input_crn, ts_file, pilfile=None, domfile=None, dnafile=None, verbose = False):
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

  if not os.path.isfile(ts_file):
    builtin = 'schemes/' + ts_file

    try :
      ts_file = pkg_resources.resource_filename('nuskell', builtin)
      print "Using scheme:", ts_file
    except KeyError :
      schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
      raise InvalidSchemeError(ts_file, schemedir)


  ts = parse_ts_file(ts_file)
  (crn, formal_species, const_species) = parse_crn_string(input_crn)

  solution, constant_solution = interpret(ts, crn, 
      formal_species, const_species)

  if pilfile :
    with open(pilfile, 'w') as pil:
      TestTubeIO(solution).write_pilfile(pil)
  if domfile :
    with open(domfile, 'w') as dom:
      TestTubeIO(solution).write_domfile(dom)
  if dnafile :
    with open(dnafile, 'w') as dna:
      TestTubeIO(solution).write_dnafile(dna, formal=formal_species, crn=crn, ts=os.path.basename(ts_file))

  return solution, constant_solution

