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

  def test_interface(self):
    """A sample test to aggree on a new interface for bisimulation.  

    Simply switch to old=False.
    """
    fcrn = "A->B"
    ecrn = "A<=>i19; i19<=>i39+X; i39->i71+i72"

    (fcrn, fs, _) = parse_crn_string(fcrn) 
    fcrn = split_reversible_reactions(fcrn)

    (ecrn, _, _) = parse_crn_string(ecrn) 
    ecrn = split_reversible_reactions(ecrn)

    partial = dict()
    partial['A'] = Counter(A=1)
    partial['B'] = Counter(B=1)

    fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
    ecrn = [[Counter(part) for part in rxn] for rxn in ecrn]

    # TODO
    out = bisimulation.test(fcrn, ecrn, fs, interpretation=partial, 
            permissive='loop-search', verbose=False)
    self.assertTrue(out[0])

    out = bisimulation.test(fcrn, ecrn, fs, interpretation=partial, 
            permissive='whole-graph', verbose=False)
    self.assertTrue(out[0])

    out = bisimulation.test(fcrn, ecrn, fs, interpretation=partial, 
            permissive='depth-first', permissive_depth=8, verbose=False)
    self.assertTrue(out[0])

    # A function that does not say so, should not modify its arguments:
    argcheck = dict()
    argcheck['A'] = Counter(A=1)
    argcheck['B'] = Counter(B=1)
    self.assertDictEqual(partial, argcheck)

  def test_example(self):
    """A simple example of finding a bisimulation from group meeting."""

    fcrn = "A + B -> C + D ; A + C -> B + D"
    icrn = "x1 -> x2 ; x3 + x4 <=> x5 ; x2 -> x6 + x8 ; x5 -> x7 ; " + \
           "x3 <=> x6 ; x9 <=> x10 ; x10 + x4 <=> x1 ; x7 -> x9 + x8"

    (fcrn, fs, _) = parse_crn_string(fcrn)
    fcrn = split_reversible_reactions(fcrn)
    fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]

    (icrn, _, _) = parse_crn_string(icrn)
    icrn = split_reversible_reactions(icrn)
    icrn = [[Counter(part) for part in rxn] for rxn in icrn]

    v, _ = bisimulation.test(fcrn, icrn, fs)
    self.assertTrue(v)

    # Test wrong partial interpretation
    partial = dict()
    partial['x2'] = Counter(['B','D'])
    partial['x3'] = Counter(['C'])

    out = bisimulation.test(fcrn, icrn, fs, interpretation=partial,
                            permissive='loop-search', verbose=False)
    self.assertFalse(out[0])
    self.assertTrue(out[1][1] >= 0)

    del partial['x3']

    out = bisimulation.test(fcrn, icrn, fs, permissive='loop-search',
                            interpretation=partial, verbose=False)
    self.assertTrue(out[0])

  def test_weird_behavior(self):
    """ a follow-up on testing the groupmeeting example """

    fcrn = "A + B -> C + D ; A + C -> B + D"
    icrn = "x1 -> x2 ; x3 + x4 <=> x5 ; x2 -> x6 + x8 ; x5 -> x7 ; " + \
           "x3 <=> x6 ; x9 <=> x10 ; x10 + x4 <=> x1 ; x7 -> x9 + x8"

    # First correct interpretation
    inter1 = {'x1': Counter({'A': 1, 'B': 1}), 
              'x2': Counter({'C': 1, 'D': 1}), 
              'x3': Counter({'C': 1}), 
              'x4': Counter({'A': 1}), 
              'x5': Counter({'A': 1, 'C': 1}),
              'x6': Counter({'C': 1}), 
              'x7': Counter({'B': 1, 'D': 1}), 
              'x8': Counter({'D': 1}), 
              'x9': Counter({'B': 1}), 
              'x10': Counter({'B': 1})} 
    pinter1= {'x7': Counter({'B': 1, 'D': 1})} 

    # Second correct interpretation
    inter2 = {'x1': Counter({'A': 1, 'C': 1}), 
              'x2': Counter({'B': 1, 'D': 1}), 
              'x3': Counter({'B': 1}), 
              'x4': Counter({'A': 1}), 
              'x5': Counter({'A': 1, 'B': 1}),
              'x6': Counter({'B': 1}), 
              'x7': Counter({'C': 1, 'D': 1}), 
              'x8': Counter({'D': 1}), 
              'x9': Counter({'C': 1}), 
              'x10': Counter({'C': 1})}
    pinter2= {'x7': Counter({'C': 1, 'D': 1})}

    # CRN preprocessing
    (fcrn, fs, _) = parse_crn_string(fcrn)
    fcrn = split_reversible_reactions(fcrn)
    fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]

    (icrn, _, _) = parse_crn_string(icrn)
    icrn = split_reversible_reactions(icrn)
    icrn = [[Counter(part) for part in rxn] for rxn in icrn]

    # NOTE: Correct behavior
    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter1,
                            permissive='whole-graph')
    self.assertTrue(v)
    self.assertDictEqual(inter1, i1)

    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter1,
                            permissive='loop-search')
    self.assertTrue(v)
    self.assertDictEqual(inter1, i1)

    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter1,
                            permissive='depth-first')
    self.assertTrue(v)
    self.assertDictEqual(inter1, i1)

    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                              permissive='loop-search')
    self.assertTrue(v)
    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                            permissive='whole-graph')
    self.assertTrue(v)
    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                            permissive='depth-first')
    self.assertTrue(v)

    del inter1['x3']
    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                            permissive='loop-search')
    self.assertTrue(v)

    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                            permissive='whole-graph')
    self.assertTrue(v)
    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                            permissive='depth-first')
    self.assertTrue(v)

    def dont_test_equivalence(self):
      (fcrn, fs) = self.equivalence['formal'][0]
      (ecrn, fs) = self.equivalence['qian2011_gen'][0]

