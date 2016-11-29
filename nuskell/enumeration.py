# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@caltech.edu)
#
# Preprocessing and interface to peppercorn enumerator
#

from nuskell.objects import TestTube, Complex

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

class TestTubePeppercornIO(object):
  """ A TestTube wrapperclass to communicate with the peppercorn enumerator.

  This class reads a Nuskell TestTube object and translates it into a
  Peppercorn Enumerator object. Both objects refer to the same complexes and
  domains.
  """
  def __init__(self, ttube, args=None):
    # Input
    self._testtube = TestTube()
    self._testtube += ttube
    self._enumerator = self._initialize_peppercorn(ttube)
    if args :
      set_enum_args(self._enumerator, args)

    # If True, the enumerator has been called and the TestTube is up-to-date.
    # only flips if used by the enumerate() function in *this* instance.
    self._processed = False 
    self._enum_to_ttube = dict()

    # Output: Autoname will reset the naming of complexes of the peppercorn
    # enumerator format and replace it with autmatic TestTube names. This is
    # nice for readability, but it becomes harder to debug code.
    self._autoname = False
    self._enum_prefix = 'e'
    self._enum_append = '' if self._autoname else '_'

  @property
  def testtube(self):
    if not self._processed and self._enumerator._resting_states != None:
      raise ValueError('Enumerator has to be called within the IO object.')
    return self._testtube

  @property
  def enumerator(self):
    return self._enumerator

  @property
  def condense_reactions(self):
    condensed = condense_resting_states(self.enumerator, compute_rates=True, k_fast=0.)
    reactions = condensed['reactions']

    enum_crn = []
    for r in sorted(reactions):
      # print map(lambda x: map(str, x.complexes), r.reactants), '->', 
      # print map(lambda x: map(str, x.complexes), r.products)
      # print r.kernel_string()

      react = []
      for rs in r.reactants:
        react.append(self._enum_to_ttube[rs.complexes[0].name])
        #if len(rs.complexes) > 1 :
        #  print 'resting state', rs.complexes, react[-1]

      prod = []
      for rs in r.products:
        prod.append(self._enum_to_ttube[rs.complexes[0].name])
        #if len(rs.complexes) > 1 :
        #  print 'resting state', rs.complexes

      enum_crn.append([react, prod])
    return enum_crn

  @property
  def all_reactions(self):
    enum_crn = []
    for r in sorted(self.enumerator.reactions):
      react = []
      for re in map(str, r.reactants) :
        if re in self._enum_to_ttube:
          re = self._enum_to_ttube[re]
        react.append(re)
      prod = []
      for pr in map(str, r.products) :
        if pr in self._enum_to_ttube:
          pr = self._enum_to_ttube[pr]
        prod.append(pr)
      #print r.kernel_string()
      print [react, prod]
      enum_crn.append([react, prod])
    return enum_crn

  def enumerate(self):
    # Dangerous, if you were to call self._enumerator.enumerate() directly,
    self._enumerator.enumerate()
    self._get_resting_cplxs()
    self._processed = True

  def _initialize_peppercorn(self, testtube):
    """Initialize the peppercorn enumerator object.

    Args:
      testtube <nuskell.objects.TestTube()>: complexes, strands and domains
      args <arparse.parser()>: arguments for the enumerator

    Returns:
      Enumerator <peppercorn.enumerator.Enumerator()>
    """

    # Translate to peppercorn domains
    domains = {}
    for n, d in testtube.domains.items() :
      if n[-1] == '*' :
        new_dom = peputils.Domain(n[:-1], d.length, 
            sequence=''.join(d.sequence), is_complement=True)
      else :
        new_dom = peputils.Domain(n, d.length, sequence=''.join(d.sequence))
      domains[n] = new_dom
    #print domains.values()

    # Translate to peppercorn strands
    strands = {}
    dom_to_strand = {}
    for n, s in testtube.strands.items() :
      dom_to_strand[tuple(map(str,s))] = n
      doms = []
      for d in map(str, s) :
        doms.append(domains[d])
      strands[n] = peputils.Strand(n, doms)
    #print strands.values()

    # Translate to peppercorn complexes
    complexes = {}
    for n, c in testtube.complexes.items():
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

  def _ecplx_rename(self, x):
    """ A function to rename enumerator species names to a format 
    which is compatible with nuskell.objects
    """
    assert x[0].isdigit()
    assert x[-1].isdigit()
    return self._enum_prefix + x + self._enum_append

  def _get_resting_cplxs(self):
    """Merge self.enumerator complexes into the self.testtube object. 

    Note that the enumerator cannot generate new domains or new strands, it
    only finds new complexes. This function adds all complexes that are
    present in resting states in the system. Transient states are ignored.
    """
    oldcplxs = self._testtube.complexes.keys()
    if not self._autoname and any(
        map(lambda x:(x[0]==self._enum_prefix), oldcplxs)) :
      raise Exception(self._enum_prefix, 
          "Namespace prefix reserved for enumerated complexes")

    domains = self._testtube.domains
    for rs in self._enumerator.resting_states:
      for cx in rs.complexes:
        if cx.name in oldcplxs: 
          self._enum_to_ttube[cx.name] = cx.name
          continue
        domseq = []
        for sd in cx.strands:
          domseq.append('+')
          for do in sd.domains:
            assert (do.name in domains)
            dom = domains[do.name]
            domseq.append(dom)
        # remove the first '+' again
        if len(domseq)>0: domseq = domseq[1:]
        domstr = list(cx.dot_paren_string())
        if self._autoname :
          mycplx = Complex(sequence=domseq, structure=domstr, 
              prefix=self._enum_prefix)
        else :
          cplxname = self._ecplx_rename(cx.name)
          if cplxname in self._testtube.complexes:
            raise ValueError("Complex found in muliple resting states?")
          mycplx = Complex(sequence=domseq, structure=domstr, name=cplxname)
        self._testtube.add_complex(mycplx)
        self._enum_to_ttube[cx.name] = mycplx.name
    return 

def set_enum_args(enum, args):
  """Transfer options to self._enumerator object. 
  
  Do NOT set Nuskell-defaults values here. Defaults are set with the argparse
  object of peppercorn: nuskell/include/peppercorn/enumerator.py 
  
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

