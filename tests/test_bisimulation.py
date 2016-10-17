
import unittest
import collections as c 

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

    self.equivalence = c.defaultdict(list)
    (fcrn, formal, constant) = parse_crn_file('tests/crns/roessler_formal.crn')
    (icrn, waste1, waste2) = parse_crn_file('tests/crns/roessler_qian2011_gen.crn')
    self.equivalence['formal'].append((fcrn, formal))
    self.equivalence['qian2011_gen'].append((icrn, formal))

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
      # Default: verbose=True
    else :
      # Test with default optional arguments:
      # verbose = False
      # method = ??
      # formal_species = None
      # interpretation = {}
      self.assertTrue(bisimulation.test(fcrn, ecrn))

      partial = c.defaultdict(c.Counter)
      partial['A'] = c.Counter(A=1)
      partial['B'] = c.Counter(B=1)

      self.assertTrue(bisimulation.test(fcrn, ecrn, formal_species=fs,
        iterpretation=partial, method='pspace', verbose=False))

      # A function that does not say so, should not modify its arguments:
      argcheck = c.defaultdict(c.Counter)
      argcheck['A'] = c.Counter(A=1)
      argcheck['B'] = c.Counter(B=1)
      self.assertDictEqual(partial, argcheck)

  def dont_test_equivalence(self):
    #TODO: This takes too long, but it should verify
    (fcrn, fs) = self.equivalence['formal'][0]
    (ecrn, fs) = self.equivalence['qian2011_gen'][0]

    self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), verbose=False))

  def test_example_groupmeeting(self):
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

    #TODO: call with the correct interpretation:
    partial = c.defaultdict(c.Counter)
    partial['x1'] = c.Counter(B=1, D=1)
    partial['x2'] = c.Counter(B=1, D=1)
    partial['x3'] = c.Counter(B=1)
    partial['x4'] = c.Counter(A=1)
    partial['x5'] = c.Counter(A=1,B=1)
    partial['x6'] = c.Counter(B=1)
    partial['x7'] = c.Counter(A=1,B=1)
    partial['x8'] = c.Counter(D=1)
    partial['x9'] = c.Counter(C=1)
    partial['x10']= c.Counter(C=1)

    inter = []
    inter.append([['x1'], [['B'], ['D']]])
    inter.append([['x2'], [['B'], ['D']]])
    inter.append([['x3'], [['B']]])
    inter.append([['x4'], [['A']]])
    inter.append([['x5'], [['A','B']]])
    inter.append([['x6'], [['B']]])
    inter.append([['x7'], [['A'],['B']]])
    inter.append([['x8'], [['D']]])
    inter.append([['x9'], [['C']]])
    inter.append([['x10'],[['C']]])
    #self.assertTrue(bisimulation.test((fcrn,fs), (ecrn,fs), inter=inter, verbose=False))

