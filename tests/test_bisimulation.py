import os
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
    # Skip slow unittests unless you have a lot of time.
    self.skip_slow = True

  def tearDown(self):
    # clean up even if unittests failed
    pass

  def _parse_crn_file(self, filename):
    # Load data from filename into the proper bisimulation format
    if not os.path.isfile(filename):
      raise Exception("File for unittest is missing: {}".format(filename))
    (crn, formal, _) = parse_crn_file(filename)
    crn = split_reversible_reactions(crn)
    return ([[Counter(part) for part in rxn] for rxn in crn], formal)

  def _parse_crn_string(self, string):
    (crn, formal, _) = parse_crn_string(string)
    crn = split_reversible_reactions(crn)
    return ([[Counter(part) for part in rxn] for rxn in crn], formal)

  def test_equivalence(self):
    if self.skip_slow : return
    (fcrn, fs) = self._parse_crn_file('tests/crns/roessler_formal.crn')
    (icrn, _) = self._parse_crn_file('tests/crns/roessler_qian2011_gen.crn')

    v, i = bisimulation.test(fcrn, icrn, fs)
    self.assertTrue(v)

  def test_interface(self):
    #A sample test to aggree on a new interface for bisimulation.  
    fcrn = "A->B"
    ecrn = "A<=>i19; i19<=>i39+X; i39->i71+i72"

    (fcrn, fs) = self._parse_crn_string(fcrn) 
    (ecrn, _) = self._parse_crn_string(ecrn) 

    partial = dict()
    partial['A'] = Counter(A=1)
    partial['B'] = Counter(B=1)

    v, i = bisimulation.test(fcrn, ecrn, fs, interpretation=partial, 
            permissive='loop-search', verbose=False)
    self.assertTrue(v)

    v, i = bisimulation.test(fcrn, ecrn, fs, interpretation=partial, 
            permissive='whole-graph', verbose=False)
    self.assertTrue(v)

    v, i = bisimulation.test(fcrn, ecrn, fs, interpretation=partial, 
            permissive='depth-first', permissive_depth=8, verbose=False)
    self.assertTrue(v)

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

    v, i = bisimulation.test(fcrn, icrn, fs, interpretation=partial,
                            permissive='loop-search', verbose=False)
    self.assertFalse(v)

    del partial['x3']

    v, _ = bisimulation.test(fcrn, icrn, fs, permissive='loop-search',
                            interpretation=partial, verbose=False)
    self.assertTrue(v)

  def test_example_with_results(self):
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
    self.assertDictEqual(inter1, i1)

    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                            permissive='whole-graph')
    self.assertTrue(v)
    self.assertDictEqual(inter1, i1)

    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                            permissive='depth-first')
    self.assertTrue(v)
    self.assertDictEqual(inter1, i1)

    v, i2 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter2,
                            permissive='loop-search')
    self.assertTrue(v)
    self.assertDictEqual(inter2, i2)

    v, i2 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter2,
                            permissive='whole-graph')
    self.assertTrue(v)
    self.assertDictEqual(inter2, i2)

    v, i2 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter2,
                            permissive='depth-first')
    self.assertTrue(v)
    self.assertDictEqual(inter2, i2)

  def test_mini1(self):
    # Two valid interpretations
    fcrn = "B + B -> B"
    icrn = "B <=> x1; B + x1 -> x2 + x3; x2 -> B + x4"

    (fcrn, fs) = self._parse_crn_string(fcrn)
    (icrn, _) = self._parse_crn_string(icrn)

    inter = {'B': Counter(['B']),
             'x2': Counter(['B','B'])}

    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation = inter)
    self.assertTrue(v)

    inter = {'B': Counter(['B']),
             'x2': Counter(['B'])}

    v, i2 = bisimulation.test(fcrn, icrn, fs, interpretation = inter)
    self.assertTrue(v)

  def test_species_names(self):
    #TODO: naming species in certain ways breaks bisimulation
    fcrn = "A+C->A+B"

    #NOTE: replace x2 with e45 and it will not terminate!
    #icrn = "A <=> x1 + x2; C+x1 <=> x3 + x4; x3 -> A + B + x5"
    icrn = "A <=> x1 + e45; C + x1 <=> x3 + x4; x3 -> A + B + x5"

    (fcrn, fs) = self._parse_crn_string(fcrn)
    (icrn, _) = self._parse_crn_string(icrn)

    inter = {'A': Counter(['A']),
             'B': Counter(['B']),
             'C': Counter(['C'])}

    v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation = inter)
    self.assertTrue(v)

if __name__ == '__main__':
  unittest.main()

