# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@caltech.edu)
#
# Preprocessing and interface to peppercorn enumerator
#

from nuskell.objects import TestTube, Complex, Reaction

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
  """ A Wrapperclass to communicate between ``Nuskell'' and ``Peppercorn''.

  This class reads a ``Nuskell'' TestTube() object, translates it into a
  ``Peppercorn'' Enumerator() object, enumerates, and updates both the
  Enumerator() and the TestTube() object with the new reactions.

  Alternatively, it may simply read a ``Peppercorn'' Enumerator() object and
  return a ``Nuskell'' TestTube() object. For example, one might start with an
  enumerated DSD network, transfer it to a TestTube() object for sequence design
  and then load it back into the Enumerator().

  Both objects refer to the same complexes and domains.
  """

  condensed = True

  def __init__(self, testtube=None, enumerator=None, pargs=None, nargs=None,
      rename = False):

    if not (bool(testtube) != bool(enumerator)) : # NOT XOR
      raise ValueError("Need to specify either TestTube() or Enumerator(), but not both!")
    
    # Autoname resets Peppercorn() complex names and replace it with TestTube()
    # complex names. Two dictionaries store the mapping of names in TestTube()
    # and Enumerator() objects.
    self._autoname = rename
    self._enum_prefix = 'e'
    self._enum_append = '' if self._autoname else '_'
    self._enumN_to_ttubeO = dict()
    self._ttubeN_to_enumO = dict()

    # If True, the enumerator has been called and the TestTube is up-to-date.
    # only flips if used by the enumerate() function in *this* instance.
    self._processed = False 

    if testtube :
      self._testtube = testtube
      self._enumerator = self.testtube_to_enumerator(self._testtube) #NOTE Complexes only!!!
    elif enumerator :
      raise NotImplementedError
      self._enumerator = enumerator
      self._testtube = self.enumerator_to_testtube(self._enumerator)

    set_peppercorn_args(self._enumerator, pargs)

  @property
  def testtube(self):
    if not self._processed and self._enumerator._resting_states != None:
      raise ValueError('Enumerator has to be called within the IO object.')
    return self._testtube

  @property
  def enumerator(self):
    return self._enumerator

  def enumerate(self):
    # You may call self._enumerator.enumerate() directly, but then the testtube
    # object will not be updated afterwards. So this wrapper makes sure
    # everythings done properly, and signals that using self._processed=True
    self._enumerator.enumerate()
    self._testtube += self.enumerator_to_testtube(self._enumerator)
    self._processed = True

  def testtube_to_enumerator(self, testtube): # Does not add reactions!
    """Initialize the peppercorn enumerator object.

    Args:
      testtube <nuskell.objects.TestTube()>: complexes, strands and domains
      tN_to_eO <dict()>: Mapping from testtube names to enumerator Objects()

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
    for cplx in testtube.complexes:
      cplx_strands = []
      for s in cplx.lol_sequence :
        ns = dom_to_strand[tuple(map(str,s))]
        cplx_strands.append(strands[ns])
      complex_structure = peputils.parse_dot_paren(''.join(cplx.structure))
      complex = peputils.Complex(cplx.name, cplx_strands, complex_structure)
      complex.check_structure()
      complexes[cplx.name] = complex		

      self._enumN_to_ttubeO[cplx.name] = cplx
      self._ttubeN_to_enumO[cplx.name] = complex
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

  def _get_reaction_networkx(self):
    """ Merge the reaction network into self.testtube object. """
    pass

  def enumerator_to_testtube(self, enumerator):
    """Merge Enumerator() complexes into the TestTube() object. 

    Note that the enumerator cannot generate new domains or new strands, it
    only finds new complexes. This function adds all complexes that are
    present in resting states in the system. Transient states are ignored.
    """
    condensed = TestTubePeppercornIO.condensed
    ttcomplexes = dict()

    domains = enumerator.domains
    # The enumerator cannot find new domains, so if there is a testtube object
    # already, then we might as well save some effort and use the testtube
    # domains directly
    if self._testtube : 
      dNames = map(lambda x: x.name, domains)
      tDomains = self._testtube.domains
      assert all(map(lambda d: d in tDomains, dNames))
      domains = tDomains
    else :
      # translate domains to nuskell.Domain format, easy...
      raise NotImplementedError

    if condensed :
      enum_complexes = enumerator.resting_complexes
    else :
      # resting_complexes + transient_complexes
      enum_complexes = enumerator.complexes

    for cx in enum_complexes:
      if cx.name in self._enumN_to_ttubeO :
        ttcplx = self._enumN_to_ttubeO[cx.name] 
        ttcomplexes[cx.name] = (ttcplx, None, None)
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
        ttcplx = Complex(sequence=domseq, structure=domstr, 
            prefix=self._enum_prefix)
      else :
        cplxname = self._ecplx_rename(cx.name)
        if cplxname in ttcomplexes:
          raise ValueError("Complex found in muliple resting states?")
        ttcplx = Complex(sequence=domseq, structure=domstr, name=cplxname)

      ttcomplexes[ttcplx.name] = (ttcplx, None, None) 
      self._enumN_to_ttubeO[cx.name] = ttcplx
      self._ttubeN_to_enumO[ttcplx.name] = cx

    crn = dict()
    if condensed:
      condensed = condense_resting_states(enumerator, 
          compute_rates = True, k_fast=enumerator.k_fast)
      reactions = condensed['reactions']
      for r in reactions:
        # NOTE: takes only the *first complex* of resting states as representative
        react = []
        for rs in r.reactants:
          react.append(self._enumN_to_ttubeO[rs.complexes[0].name])
        prod = []
        for rs in r.products:
          prod.append(self._enumN_to_ttubeO[rs.complexes[0].name])
        rxn = Reaction(react, prod, rate=r.rate(), rtype = r.name)
        crn[rxn.name] = rxn
    else :
      for r in enumerator.reactions:
        react = []
        for re in r.reactants :
          react.append(self._enumN_to_ttubeO[re.name])
        prod = []
        for pr in r.products :
          prod.append(self._enumN_to_ttubeO[pr.name])
        rxn = Reaction(react, prod, rate=r.rate(), rtype = r.name)
        crn[rxn.name] = rxn

    return TestTube(complexes=ttcomplexes, reactions=crn)

def set_peppercorn_args(enum, args):
  """Transfer options to self._enumerator object. 

  Do NOT change default values here. These are supposed to be the defaults of
  peppercorn!  Defaults for nuskell or any other script using this library are
  set with the argparse object of your script, e.g. nuskell: scripts/nuskell.
  
  """

  if hasattr(args, 'verbose'):
    import logging
    logger = logging.getLogger()
    if args.verbose == 1:
      logger.setLevel(logging.INFO)
    elif args.verbose == 2:
      logger.setLevel(logging.DEBUG)
    elif args.verbose >= 3:
      logger.setLevel(logging.NOTSET)

  if hasattr(args, 'MAX_COMPLEX_SIZE'):
    enum.MAX_COMPLEX_SIZE = args.MAX_COMPLEX_SIZE
  else :
    enum.MAX_COMPLEX_SIZE = 100

  if hasattr(args, 'MAX_COMPLEX_COUNT'):
    enum.MAX_COMPLEX_COUNT = args.MAX_COMPLEX_COUNT
  else :
    enum.MAX_COMPLEX_COUNT = 1000

  if hasattr(args, 'MAX_REACTION_COUNT'):
    enum.MAX_REACTION_COUNT = args.MAX_REACTION_COUNT
  else :
    enum.MAX_REACTION_COUNT = 5000

  if hasattr(args, 'REJECT_REMOTE'):
    enum.REJECT_REMOTE = args.REJECT_REMOTE
  else :
    enum.REJECT_REMOTE = False

  if hasattr(args, 'ignore_branch_3way') and args.ignore_branch_3way:
    if reactions.branch_3way in enum.FAST_REACTIONS:
      enum.FAST_REACTIONS.remove(reactions.branch_3way)

  if hasattr(args, 'ignore_branch_4way') and args.ignore_branch_4way:
    if reactions.branch_4way in enum.FAST_REACTIONS:
      enum.FAST_REACTIONS.remove(reactions.branch_4way)

  if hasattr(args, 'RELEASE_CUTOFF_1_1'):
    enum.RELEASE_CUTOFF_1_1 = args.RELEASE_CUTOFF_1_1
  else :
    enum.RELEASE_CUTOFF_1_1 = 6

  if hasattr(args, 'RELEASE_CUTOFF_1_N'):
    enum.RELEASE_CUTOFF_1_N = args.RELEASE_CUTOFF_1_N
  else :
    enum.RELEASE_CUTOFF_1_N = 6

  if hasattr(args, 'RELEASE_CUTOFF'):
    if args.RELEASE_CUTOFF is not None :
      enum.RELEASE_CUTOFF_1_1 = args.RELEASE_CUTOFF 
      enum.RELEASE_CUTOFF_1_N = args.RELEASE_CUTOFF
      enum.RELEASE_CUTOFF = args.RELEASE_CUTOFF
  else :
    enum.RELEASE_CUTOFF = None

  if hasattr(args, 'UNZIP'):
    enum.UNZIP = args.UNZIP
  else :
    enum.UNZIP = True

  if hasattr(args, 'LEGACY_UNZIP'):
    enum.LEGACY_UNZIP = args.LEGACY_UNZIP
  else :
    enum.LEGACY_UNZIP = False

  if hasattr(args, 'k_slow'):
    enum.k_slow = args.k_slow
  else :
    enum.k_slow = 0.0

  if hasattr(args, 'k_fast'):
    enum.k_fast = args.k_fast
  else :
    enum.k_slow = 0.0

  return

