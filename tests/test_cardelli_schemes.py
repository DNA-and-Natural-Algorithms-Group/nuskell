import unittest
from collections import Counter

from nuskell import translate, verify
from nuskell.verifier import removeSpecies
from nuskell.parser import parse_crn_string, split_reversible_reactions

@unittest.skip('slow and TODO')
class CardelliSchemes(unittest.TestCase):
  def setUp(self):
    self.cFJ_original = 'schemes/original/cardelli2011_FJ.ts'
    self.cNM_original = 'schemes/original/cardelli2011_NM.ts'
    self.c2D_original = 'schemes/original/cardelli2013_2D.ts'

  def tearDown(self):
    pass
 
  def _get_verification_data(self, input_crn, scheme):
    (fcrn, fs, _) = parse_crn_string(input_crn)
    solution, _ = translate(input_crn, scheme)
    fuels = map(str, solution.present_species(exclude=fs))
    solution.enumerate_reactions()
    interpretation = solution.interpret_species(fs, prune=True)
    icrn = []
    for r in solution.reactions:
      rxn = [map(str,r.reactants),map(str,r.products)]
      icrn.append(rxn)
    vcrn = removeSpecies(icrn, fuels)
    fcrn = split_reversible_reactions(fcrn)
    return fcrn, vcrn, fs, interpretation
 
  def test_TT_ApB_XpY_irrev(self):
    input_crn = "A + B -> X + Y"

    fcrn, vcrn, fs, interpret = self._get_verification_data(input_crn, self.cFJ_original)
   
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation')
    self.assertTrue(v)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway')
    self.assertTrue(v)

  def test_TT_ApB_XpY_rev(self):
    input_crn = "A + B <=> X + Y"

    fcrn, vcrn, fs, interpret = self._get_verification_data(input_crn, self.cFJ_original)

    #TODO: does not terminate
    #v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation')
    #self.assertTrue(v)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway')
    self.assertTrue(v)

  def test_FF_A_B_irrev(self):
    input_crn = "A -> B"

    fcrn, vcrn, fs, interpret = self._get_verification_data(input_crn, self.cFJ_original)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation')
    self.assertFalse(v)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway')
    self.assertFalse(v)

  def test_TF_A_B_rev(self):
    input_crn = "A <=> B"

    fcrn, vcrn, fs, interpret = self._get_verification_data(input_crn, self.cFJ_original)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation')
    self.assertTrue(v)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway')
    self.assertFalse(v)

  def test_FF_A_f_rev(self):
    input_crn = "A <=> "
    fcrn, vcrn, fs, interpret = self._get_verification_data(input_crn, self.cFJ_original)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation')
    self.assertTrue(v)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway')
    self.assertFalse(v)

    fcrn, vcrn, fs, interpret = self._get_verification_data(input_crn, self.cNM_original)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation')
    self.assertTrue(v)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway')
    self.assertTrue(v)

  def test_FF_ApB_f_rev(self):
    input_crn = "A +B <=> "
    fcrn, vcrn, fs, interpret = self._get_verification_data(input_crn, self.c2D_original)

    #TODO: does not terminate
    #v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation')
    #self.assertTrue(v)

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway')
    self.assertFalse(v)

if __name__ == '__main__':
  unittest.main()

