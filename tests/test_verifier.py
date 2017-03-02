
import unittest
import collections as c 

from nuskell.parser import parse_crn_string, split_reversible_reactions
from nuskell.parser import parse_crn_file
from nuskell.objects import Complex
import nuskell.verifier.verifier as verifier

class VerificationPreprocessingTests(unittest.TestCase):
  #Tests the processing of enumerated CRNs prior to verification.

  def setUp(self):
    # preprocessing for unittesting
    self.cplx1 = Complex(sequence=['?', 't0', 'd2'], structure=['.', '.', '.'])
    self.cplx2 = Complex(sequence=['d19', 't0', 'd2'], structure=['.', '.', '.'])
    self.cplx3 = Complex(sequence=['d19', 't1', 'd2'], structure=['.', '.', '.'])
    self.cplx4 = Complex(sequence=['?', 't0', 'd2'], structure=['.', '(', '('])
    self.cplx5 = Complex(sequence=['d19', '?', 'd2'], structure=['.', '.', '.'])

    self.cplx11 = Complex(sequence=['?', 't2', 't3', '+', 'd4', 't5', 'd11', 't0', 'd12', 't6', '+', 't5*', 'd4*', 't3*', 't2*', 'd1i*', 't0*'], 
        structure = ['.', '(', '(', '+', '(', '(', '.', '.', '.', '.', '+', ')', ')', ')', ')', '.', '.'])
    self.cplx22 = Complex(sequence=['d1', 't2', 't3', '+', 'd4', 't5', 'd11', 't0', 'd12', 't6', '+', 't5*', 'd4*', 't3*', 't2*', 'd1i*', 't0*'], 
        structure = ['.', '(', '(', '+', '(', '(', '.', '.', '.', '.', '+', ')', ')', ')', ')', '.', '.'])
    self.cplx33 = Complex(sequence=['d4', 't5', 'd11', 't0', 'd12', 't6', '+', 't5*', 'd4*', 't3*', 't2*', 'd1i*', 't0*', '+', 'd1', 't2', 't3'], 
        structure = ['(', '(', '.', '.', '.', '.', '+', ')', ')', '(', '(', '.', '.', '+', '.', ')', ')'])

    self.input_fs = None
    self.init_cplxs = None
    self.enum_cplxs = None

  def tearDown(self):
    # clean up even if unittests failed
    pass

  @unittest.skip("Depricated")
  def test_patternMatch(self):
    # single stranded complexes
    self.assertTrue(verifier.patternMatch(self.cplx1,self.cplx2))
    self.assertTrue(verifier.patternMatch(self.cplx1,self.cplx5))
    self.assertFalse(verifier.patternMatch(self.cplx1,self.cplx3))
    self.assertFalse(verifier.patternMatch(self.cplx4,self.cplx2))

    # single stranded complexes
    self.assertTrue(verifier.patternMatch(self.cplx11,self.cplx22))
    self.assertTrue(verifier.patternMatch(self.cplx11,self.cplx33))

  def test_removeDuplicates(self):
    crn1 = "A -> B; B <=> A; C -> D; B -> A + X; A -> A; Y + B -> Z + i + A; B + Y -> i + A + Z"
    crn2 = "A -> A; A -> B; B -> A; B -> A + X; B + Y -> A + Z + i; C -> D"

    (crn, fs, cs) = parse_crn_string(crn1)
    ir_crn1 = split_reversible_reactions(crn)

    (crn, fs, cs) = parse_crn_string(crn2)
    ir_crn2 = split_reversible_reactions(crn)

    self.assertEqual(verifier.removeDuplicates(ir_crn1), ir_crn2)


if __name__ == '__main__':
  unittest.main()

