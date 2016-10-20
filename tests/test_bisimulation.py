
import unittest
import collections as c 

from nuskell.parser import parse_crn_string, split_reversible_reactions
import nuskell.verifier.crn_bisimulation_equivalence as bisimulation

class BisimulationTests(unittest.TestCase):
  """Bisimulation Testing Class:

  Compares *formal* CRNs with *enumerated* CRNs.

  Note: Translation and enumeration are not part of this testing class. The
        correct translation/enumeration has to be checked elsewhere! It is not
        necessary to refer to translation schemes here at all, but it does make
        sense for reproducability.
  """
  def setUp(self):
    # preprocessing for unittesting

    # e.g. load data into dictionary, maybe from a crn directory:
    # crn/roessler_formal.crn
    # crn/roessler_qian2011.crn
    # crn/roessler_qian20xx.crn 

    # self.crnequiv['formal']   = [[crn1], [crn2], [crn3]]
    # self.crnequiv['qian2011'] = [[crn1], [crn2], [crn3]]
    pass

  def tearDown(self):
    # clean up even if unittests failed
    pass

  def test_interface(self):
    """A sample test to aggree on a new interface for bisimulation.  

    Simply switch to old=False.
    """
    fcrn = "A->B"
    ecrn = "A<=>i19; i19<=>i39+X; i39->i71+i72"

    (fcrn, fs, cs) = parse_crn_string(fcrn) 
    fcrn = split_reversible_reactions(fcrn)

    (ecrn, w1, w2) = parse_crn_string(ecrn) 
    ecrn = split_reversible_reactions(ecrn)

    old = False
    if old :
      self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False))
      self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False, permcheck='pspace'))
      self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False, permcheck='whole'))
    else :
      partial = dict()
      partial['A'] = c.Counter(A=1)
      partial['B'] = c.Counter(B=1)

      fcrn = [[c.Counter(part) for part in rxn] for rxn in fcrn]
      ecrn = [[c.Counter(part) for part in rxn] for rxn in ecrn]

      # Since the interpretation maps species of eCRN to species of fCRN, I
      # would prefer to have that order in the arguments.
      out = bisimulation.test(fcrn, ecrn, fs,
            interpretation=partial, permissive='loop-search',
                              verbose=False)
      self.assertTrue(out[0])
      out = bisimulation.test(fcrn, ecrn, fs,
            interpretation=partial, permissive='whole-graph',
            verbose=False)
      self.assertTrue(out[0])
      out = bisimulation.test(fcrn, ecrn, fs,
            interpretation=partial, permissive='depth-first',
                              permissive_depth=8, verbose=False)
      self.assertTrue(out[0])

      # A function that does not say so, should not modify its arguments:
      argcheck = dict()
      argcheck['A'] = c.Counter(A=1)
      argcheck['B'] = c.Counter(B=1)
      self.assertDictEqual(partial, argcheck)

  def test_example(self):
    '''A simple example of finding a bisimulation from group meeting.'''

    fcrn = "A + B -> C + D ; A + C -> B + D"
    icrn = "x1 -> x2 ; x3 + x4 <=> x5 ; x2 -> x6 + x8 ; x5 -> x7 ; " + \
           "x3 <=> x6 ; x9 <=> x10 ; x10 + x4 <=> x1 ; x7 -> x9 + x8"

    (fcrn, fs, _) = parse_crn_string(fcrn)
    fcrn = split_reversible_reactions(fcrn)
    fcrn = [[c.Counter(part) for part in rxn] for rxn in fcrn]

    (icrn, _, _) = parse_crn_string(icrn)
    icrn = split_reversible_reactions(icrn)
    icrn = [[c.Counter(part) for part in rxn] for rxn in icrn]

    out = bisimulation.test(fcrn, icrn, fs)
    self.assertTrue(out[0])

    partial = dict()
    partial['x2'] = c.Counter(['B','D'])
    partial['x3'] = c.Counter(['C'])

    out = bisimulation.test(fcrn, icrn, fs, interpretation=partial,
                            permissive='loop-search', verbose=False)
    self.assertTrue(not out[0])
    self.assertTrue(out[1][1] >= 0)

    del partial['x3']

    out = bisimulation.test(fcrn, icrn, fs, permissive='loop-search',
                            interpretation=partial, verbose=False)
    self.assertTrue(out[0])
