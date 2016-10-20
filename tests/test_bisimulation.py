import unittest
from collections import Counter

from nuskell.parser import parse_crn_string, split_reversible_reactions
from nuskell.parser import parse_crn_file
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

    # Load data from crn/ dictionary into self.crnequiv dictionary:
    # self.crnequiv['formal']   = [(crn1, fs), (crn2, fs), (crn3, fs)]
    # self.crnequiv['qian2011'] = [(crn1, fs), (crn2, fs), (crn3, fs)]

    self.equivalence = dict()
    (fcrn, formal, constant) = parse_crn_file('tests/crns/roessler_formal.crn')
    (icrn, waste1, waste2) = parse_crn_file('tests/crns/roessler_qian2011_gen.crn')
    self.equivalence['formal'] = [(fcrn, formal)]
    self.equivalence['qian2011_gen'] = [(icrn, formal)]

  def tearDown(self):
    # clean up even if unittests failed
    pass

  def dont_test_interface(self):
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
      # WHAT??
      self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False))
      self.assertFalse(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False, permcheck='whole-graph'))
    else :
      # Test with default optional arguments:
      # interpretation = {}
      # permissive = 'whole-graph'
      # permissive_limit = None
      # verbose = False
      self.assertTrue(bisimulation.test(fcrn, ecrn, fs))

      partial = dict()
      partial['A'] = Counter(A=1)
      partial['B'] = Counter(B=1)

      # NOTE: The new interface might have two return values:
      # correct, best = bisimulation.test(fcrn, ecrn, fs)

      # testing permissive condition with: while-graph, no bound on nr. of cycles
      self.assertTrue(bisimulation.test(fcrn, icrn, fs, iterpretation=dict(), 
        permissive='whole-graph', permissive_limit=None, verbose=False))

      # testing permissive condition with: loop-search, no bound on len. of paths
      self.assertTrue(bisimulation.test(fcrn, icrn, fs, iterpretation=dict(), 
        permissive='loop-search', permissive_limit=None, verbose=False))

      # testing permissive condition with: depth-first, bound to 6 trivial reactions
      self.assertTrue(bisimulation.test(fcrn, icrn, fs, iterpretation=dict(), 
        permissive='depth-first', permissive_limit=6, verbose=False))

      # A function that does not say so, should not modify its arguments:
      argcheck = dict()
      argcheck['A'] = Counter(A=1)
      argcheck['B'] = Counter(B=1)
      self.assertDictEqual(partial, argcheck)

  def dont_test_equivalence(self):
    #TODO: This takes too long, but it should verify
    (fcrn, fs) = self.equivalence['formal'][0]
    (ecrn, fs) = self.equivalence['qian2011_gen'][0]

    self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False))

  def dont_test_example_groupmeeting(self):
    # give a proper citation
    fcrn = "A + B -> C + D; A + C -> B + D"
    ecrn = "x1->x2; x3+x4<=>x5; x2->x6+x8;"+ \
        "x5->x7; x3<=>x6; x9<=>x10;" + \
        "x10+x4<=>x1; x7->x9+x8"

    (fcrn, fs, cs) = parse_crn_string(fcrn) 
    fcrn = split_reversible_reactions(fcrn)

    (ecrn, w1, w2) = parse_crn_string(ecrn) 
    ecrn = split_reversible_reactions(ecrn)

    #NOTE: this should actually verify correct, shouldn't it?
    #self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False))

    #TODO: call with an interpretation:
    partial = dict(Counter)
    partial['x1'] = Counter(B=1, D=1)
    partial['x2'] = Counter(B=1, D=1)
    partial['x3'] = Counter(B=1)
    partial['x4'] = Counter(A=1)
    partial['x5'] = Counter(A=1,B=1)
    partial['x6'] = Counter(B=1)
    partial['x7'] = Counter(A=1,B=1)
    partial['x8'] = Counter(D=1)
    partial['x9'] = Counter(C=1)
    partial['x10']= Counter(C=1)

