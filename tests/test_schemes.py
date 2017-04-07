import unittest
import argparse

from nuskell import translate, verify, printCRN
from nuskell.verifier import removeSpecies, removeRates
from nuskell.parser import parse_crn_string, split_reversible_reactions

def get_verification_data(input_crn, scheme, pargs=None):
  (fcrn, fs, _) = parse_crn_string(input_crn)
  solution, _ = translate(input_crn, scheme)
  fuels = map(str, solution.present_complexes(exclude=fs))
  solution.enumerate_reactions(pargs)
  interpretation = solution.interpret_species(fs, prune=True)
  icrn = []
  for r in solution.reactions:
    rxn = [map(str,r.reactants),map(str,r.products)]
    icrn.append(rxn)
  vcrn = removeSpecies(icrn, fuels)
  fcrn = split_reversible_reactions(fcrn)
  fcrn = removeRates(fcrn)
  return fcrn, vcrn, fs, interpretation

#@unittest.skip("correct")
class SpontaneousReactions(unittest.TestCase):
  def setUp(self):
    # Enumerator-args
    parser = argparse.ArgumentParser()
    self.args = parser.parse_args([])

    # Schemes
    self.cFJ_original = 'schemes/original/cardelli2011_FJ.ts'
    self.cNM_original = 'schemes/original/cardelli2011_NM.ts'
    self.c2D_original = 'schemes/original/cardelli2013_2D.ts'
    self.cFJ_noGC     = 'schemes/original/cardelli2011_FJ_noGC.ts'
    self.cNM_noGC     = 'schemes/original/cardelli2011_NM_noGC.ts'

    self.qian2011_gen = 'schemes/generalized/qian2011_gen.ts'
    self.srin2017_phd = 'schemes/generalized/srinivas2017_phd.ts'
    self.solo2010_v1  = 'schemes/generalized/soloveichik2010_v1.ts'

    self.qian2011     = 'schemes/incorrect/qian2011.ts'

  def test_A_f_irrev(self):
    input_crn = "A -> "

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_f_A_irrev(self):
    input_crn = "-> A"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_A_f_rev(self):
    input_crn = "A <=> "

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_ApB_f_irrev(self):
    input_crn = "A + B -> "

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, None)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_f_ApB_irrev(self):
    input_crn = " -> A + B"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_ApB_f_rev(self):
    input_crn = "A + B <=> "

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertTrue(v in (True, None)) #TODO: Verification sometimes doesn't terminate
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, None)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, None)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

#@unittest.skip("correct")
class UnimolecularReactions(unittest.TestCase):
  def setUp(self):
    # Enumerator-args
    parser = argparse.ArgumentParser()
    self.args = parser.parse_args([])

    # Schemes
    self.cFJ_original = 'schemes/original/cardelli2011_FJ.ts'
    self.cNM_original = 'schemes/original/cardelli2011_NM.ts'
    self.c2D_original = 'schemes/original/cardelli2013_2D.ts'
    self.cFJ_noGC     = 'schemes/original/cardelli2011_FJ_noGC.ts'
    self.cNM_noGC     = 'schemes/original/cardelli2011_NM_noGC.ts'

    self.qian2011_gen = 'schemes/generalized/qian2011_gen.ts'
    self.srin2017_phd = 'schemes/generalized/srinivas2017_phd.ts'
    self.solo2010_v1  = 'schemes/generalized/soloveichik2010_v1.ts'

    self.qian2011     = 'schemes/incorrect/qian2011.ts'
 
  def test_A_B_irrev(self):
    input_crn = "A -> B"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_A_B_rev(self):
    input_crn = "A <=> B"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

#@unittest.skip("correct")
class BimolecularReactions(unittest.TestCase):
  def setUp(self):
    # Enumerator-args
    parser = argparse.ArgumentParser()
    self.args = parser.parse_args([])

    self.cFJ_original = 'schemes/original/cardelli2011_FJ.ts'
    self.cNM_original = 'schemes/original/cardelli2011_NM.ts'
    self.c2D_original = 'schemes/original/cardelli2013_2D.ts'
    self.cFJ_noGC     = 'schemes/original/cardelli2011_FJ_noGC.ts'
    self.cNM_noGC     = 'schemes/original/cardelli2011_NM_noGC.ts'

    self.qian2011_gen = 'schemes/generalized/qian2011_gen.ts'
    self.srin2017_phd = 'schemes/generalized/srinivas2017_phd.ts'
    self.solo2010_v1  = 'schemes/generalized/soloveichik2010_v1.ts'

    self.qian2011     = 'schemes/incorrect/qian2011.ts'

  def tearDown(self):
    self.args = None

  def test_ApB_XpY_irrev(self):
    input_crn = "A + B -> X + Y"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, None)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_ApB_XpY_rev(self):
    input_crn = "A + B <=> X + Y"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertTrue(v in (True, None)) #TODO: Verification sometimes doesn't terminate)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertTrue(v in (True, None)) #TODO: Verification sometimes doesn't terminate)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, None)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, None)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_ApB_ApA_irrev(self):
    input_crn = "A + B -> A + A"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, None)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_ApA_ApB_irrev(self):
    input_crn = "A + A -> A + B"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, None)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
  
    with self.assertRaises(RuntimeError):
      fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)

    self.args.REJECT_REMOTE = True
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, 
        self.solo2010_v1, pargs=self.args)
    self.args.REJECT_REMOTE = False

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, False)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

  def test_ApA_ApB_rev(self):
    input_crn = "A + A <=> A + B"

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, None)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, None)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, None)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.c2D_original)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, None)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cFJ_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, None)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.cNM_noGC)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, None)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, False)

    with self.assertRaises(RuntimeError):
      fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1)
    self.args.REJECT_REMOTE = True
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.solo2010_v1, pargs=self.args)
    self.args.REJECT_REMOTE = False

    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)

    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.qian2011_gen)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)
 
    fcrn, vcrn, fs, interpret = get_verification_data(input_crn, self.srin2017_phd)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='bisimulation', timeout=60)
    self.assertEqual(v, True)
    v, i = verify(fcrn, vcrn, fs, interpret=interpret, method='pathway', timeout=60)
    self.assertEqual(v, True)


if __name__ == '__main__':
  unittest.main()

