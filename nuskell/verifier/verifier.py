#!/usr/bin/env python
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@dna.caltech.edu)
#
# Verification interface
#

from collections import Counter

from nuskell.parser import combine_reversible_reactions
import crn_bisimulation_equivalence
import crn_pathway_equivalence

def printRxn(rxn, inter={}):
  # First, replace every Instance of a variable by its interpretation
  # Second, print it
  react = map(lambda x: ' + '.join(
    sorted(inter[x].elements())) if x in inter else x, rxn[0])
  produ = map(lambda x: ' + '.join(
    sorted(inter[x].elements())) if x in inter else x, rxn[1])
  if len(rxn) == 3 and rxn[2] == 'reversible':
    print ' + '.join(react), '<=>', ' + '.join(produ)
  else :
    print ' + '.join(react), '->', ' + '.join(produ)

def removeSpecies(crn, fuel):
  crn = [[filter(lambda s: s not in fuel, rxn[0]),
          filter(lambda s: s not in fuel, rxn[1])]
         for rxn in crn]
  return crn

def removeDuplicates(l):
  def helper(l) :
    r = []
    if len(l) == 0: return []
    l.sort()
    while len(l) > 1:
      if l[0] != l[1]:
        r.append(l[0])
      l = l[1:]
    r.append(l[0])
    return r
  l = sorted(map(lambda x: [sorted(x[0]), sorted(x[1])], l))
  return helper(l)

def verify(irrev_crn, enum_crn, input_fs, interpret = None, 
    method = 'bisimulation', verbose = False):
  """Wrapper to choose from different notions of equivalence.
  """
  interactive = False
  v, i = None, None

  fcrn = [[Counter(part) for part in rxn] for rxn in irrev_crn]
  ecrn = [[Counter(part) for part in rxn] for rxn in enum_crn]

  if method == 'bisimulation' or method == 'bisim-whole-graph' :
    v, i = crn_bisimulation_equivalence.test(fcrn, ecrn, input_fs, 
        interpretation = interpret, permissive='whole-graph', verbose=verbose)

  elif method == 'bisim-loop-search':
    v, i = crn_bisimulation_equivalence.test(fcrn, ecrn, input_fs,
        interpretation = interpret, permissive='loop-search', verbose=verbose)

  elif method == 'bisim-depth-first':
    v, i = crn_bisimulation_equivalence.test(fcrn, ecrn, input_fs,
        interpretation = interpret, permissive='depth-first', verbose=verbose)

  elif method == 'pathway':
    # NOTE: Adaptation to pathway interface
    pinter = dict()
    if interpret :
      for k,v in interpret.items() :
        if v :
          v = sorted(v.elements())[0]
          pinter[k]=[v]
        else :
          pinter[k]=[]
    v = crn_pathway_equivalence.test((irrev_crn, input_fs), 
        (enum_crn, pinter.keys()), pinter, False, interactive, verbose)
  elif method == 'integrated':
    # TODO: integrated-hybrid -> first consider some species as formal for
    # pathway decomposition, then do bisimulation. This is necessary for
    # history domains in some schemes, but it can be used for more general
    # things. E.g. if you make reversible reactions formal, then it will get
    # accepted and bisimulation can do the rest. In any case, the current
    # implementation does not cover the full hybrid theory, only some special
    # cases and some kind of bisimulation.
    pinter = dict()
    if interpret :
      for k,v in interpret.items() :
        if v :
          v = sorted(v.elements())[0]
          pinter[k]=[v]
        else :
          pinter[k]=[]
    v = crn_pathway_equivalence.test((irrev_crn, input_fs), 
        (enum_crn, pinter.keys()), pinter, True, interactive, verbose)
  else:
    raise RuntimeError('Unknown verification method.')

  return v, i

