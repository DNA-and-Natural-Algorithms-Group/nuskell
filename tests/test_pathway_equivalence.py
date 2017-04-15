import os
import unittest

from nuskell.parser import parse_crn_string, split_reversible_reactions
from nuskell.parser import parse_crn_file
import nuskell.verifier.crn_pathway_equivalence as pathway_equivalence

class PathwayEquivalenceTests(unittest.TestCase):
  """Pathway decomposition Testing Class:

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
    crn, formal, _, _ = parse_crn_file(filename)
    crn = split_reversible_reactions(crn)
    return ([[Counter(part) for part in rxn] for rxn in crn], formal)

  def _parse_crn_string(self, string):
    crn, formal, _, _ = parse_crn_string(string)
    crn = split_reversible_reactions(crn)
    return (crn, formal)

  def test_STW17_intro(self):
    crn1 = "A+B -> C+D; C+A -> C+C"
    crn2 = "A<=>i; i+B<=>j; i+j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
    crn3 = "A<=>i; i+B<=>j; j<=>C+k; k->D; C+A<=>m+n; m+n->C+C"
    crn4 = "A->i; i+B<=>j; j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
    crn5 = "A<=>i; i+B<=>j; j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
    crn6 = "A+g1<=>i+g2; i+B<=>j+g3; g4+j->C+k+w1; g5+k<=>D+w2; C+A<=>m+n; g6+m+n->C+C+w3"

    (crn1, fs) = self._parse_crn_string(crn1)
    (crn2, _) = self._parse_crn_string(crn2)
    (crn3, _) = self._parse_crn_string(crn3)
    (crn4, _) = self._parse_crn_string(crn4)
    (crn5, _) = self._parse_crn_string(crn5)
    (crn6, _) = self._parse_crn_string(crn6)

    inter = { 'A':'A', 'B':'B', 'C':'C', 'D':'D'}
 
    self.assertFalse(pathway_equivalence.test((crn1,fs),(crn2,fs), inter))
    self.assertFalse(pathway_equivalence.test((crn1,fs),(crn3,fs), inter))
    self.assertFalse(pathway_equivalence.test((crn1,fs),(crn4,fs), inter))
    self.assertTrue(pathway_equivalence.test((crn1,fs),(crn5,fs), inter))
    self.assertFalse(pathway_equivalence.test((crn1,fs),(crn6,fs), inter))

  def test_STW17_ex1(self):
    pass
 
if __name__ == '__main__':
  unittest.main()

