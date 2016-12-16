import unittest
from collections import Counter

from nuskell import translate, enumerateTT, verify
from nuskell.verifier import preprocess
from nuskell.parser import parse_crn_string, split_reversible_reactions

@unittest.skip("skipping slow tests")
class CardelliSchemes(unittest.TestCase):

  def setUp(self):
    self.cFJ_original = 'schemes/original/cardelli2011_FJ.ts'
    self.cNM_original = 'schemes/original/cardelli2011_NM.ts'
    self.c2D_original = 'schemes/original/cardelli2013_2D.ts'
 
  def test_TT_ApB_XpY_irrev(self):
    input_crn = "A + B -> X + Y"

    (fcrn, fs, _) = parse_crn_string(input_crn)

    scheme = self.cFJ_original

    solution, _ = translate(input_crn, scheme)
    enum_solution, enum_crn = enumerateTT(solution)

    icrn, interpret = preprocess(fcrn, enum_crn, fs, solution, enum_solution)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='bisimulation')
    self.assertTrue(v)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='pathway')
    self.assertTrue(v)

  def test_TT_ApB_XpY_rev(self):
    input_crn = "A + B <=> X + Y"
    (fcrn, fs, _) = parse_crn_string(input_crn)
    fcrn = split_reversible_reactions(fcrn)
    scheme = self.cFJ_original

    solution, _ = translate(input_crn, scheme)
    enum_solution, enum_crn = enumerateTT(solution)

    icrn, interpret = preprocess(fcrn, enum_crn, fs, solution, enum_solution)
    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='bisimulation')
    self.assertTrue(v)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='pathway')
    self.assertTrue(v)

  def test_FF_A_B_irrev(self):
    input_crn = "A -> B"
    (fcrn, fs, _) = parse_crn_string(input_crn)
    fcrn = split_reversible_reactions(fcrn)
    scheme = self.cFJ_original

    solution, _ = translate(input_crn, scheme)
    enum_solution, enum_crn = enumerateTT(solution)

    icrn, interpret = preprocess(fcrn, enum_crn, fs, solution, enum_solution)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='bisimulation')
    self.assertFalse(v)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='pathway')
    self.assertFalse(v)

  def test_TF_A_B_rev(self):
    input_crn = "A <=> B"
    (fcrn, fs, _) = parse_crn_string(input_crn)
    fcrn = split_reversible_reactions(fcrn)
    scheme = self.cFJ_original

    solution, _ = translate(input_crn, scheme)
    enum_solution, enum_crn = enumerateTT(solution)

    icrn, interpret = preprocess(fcrn, enum_crn, fs, solution, enum_solution)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='bisimulation')
    self.assertTrue(v)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='pathway')
    self.assertFalse(v)

  def test_FF_A_f_rev(self):
    input_crn = "A <=> "
    (fcrn, fs, _) = parse_crn_string(input_crn)
    fcrn = split_reversible_reactions(fcrn)

    scheme = self.cFJ_original

    solution, _ = translate(input_crn, scheme)
    enum_solution, enum_crn = enumerateTT(solution)

    icrn, interpret = preprocess(fcrn, enum_crn, fs, solution, enum_solution)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='bisimulation')
    self.assertFalse(v)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='pathway')
    self.assertFalse(v)

    #TODO: see why bisimulation is not correct!!
    scheme = self.cNM_original

    solution, _ = translate(input_crn, scheme)
    enum_solution, enum_crn = enumerateTT(solution)

    icrn, interpret = preprocess(fcrn, enum_crn, fs, solution, enum_solution)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='bisimulation')
    self.assertFalse(v)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='pathway')
    self.assertTrue(v)


  def test_FF_ApB_f_rev(self):
    input_crn = "A +B <=> "
    (fcrn, fs, _) = parse_crn_string(input_crn)
    fcrn = split_reversible_reactions(fcrn)

    scheme = self.c2D_original

    solution, _ = translate(input_crn, scheme)
    enum_solution, enum_crn = enumerateTT(solution)

    icrn, interpret = preprocess(fcrn, enum_crn, fs, solution, enum_solution)

    # TODO: takes surprisingly long!!
    #v, i = verify(fcrn, icrn, fs, interpret=interpret, method='bisimulation')
    #self.assertFalse(v)

    v, i = verify(fcrn, icrn, fs, interpret=interpret, method='pathway')
    self.assertFalse(v)



if __name__ == '__main__':
  unittest.main()

