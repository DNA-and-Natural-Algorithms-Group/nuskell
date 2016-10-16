
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

    old = True
    if old :
      self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False))
    else :
      partial = c.defaultdict(c.Counter)
      partial['A'] = c.Counter(A=1)
      partial['B'] = c.Counter(B=1)

      self.assertTrue(bisimulation.test(fcrn, ecrn, iterpretation=partial, 
              formal_species=None, method='pspace', verbose=False))

      # A function that does not say so, should not modify its arguments:
      argcheck = c.defaultdict(c.Counter)
      argcheck['A'] = c.Counter(A=1)
      argcheck['B'] = c.Counter(B=1)
      self.assertDictEqual(partial, argcheck)

