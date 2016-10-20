
import unittest
import collections as c 

from nuskell.parser import parse_crn_string, split_reversible_reactions
from nuskell.parser import parse_crn_file
import nuskell.verifier.verifier as verifier

class VerificationPreprocessingTests(unittest.TestCase):
  #Tests the processing of enumerated CRNs prior to verification.

  def setUp(self):
    # preprocessing for unittesting
    pass 

  def tearDown(self):
    # clean up even if unittests failed
    pass

  def test_rotate(self):
    pass

  def dont_test_patternMatch(self):
    """Test the patterMatch function
    """
    x = [['?', ['t0'], ['d2']], ['.', '.', '.']]
    y = [[['d19'], ['t0'], ['d2']], ['.', '.', '.']]

    self.assertTrue(verifier.patternMatch(x,y))
    self.assertTrue(verifier.pM2(x,y))

    x = [['?', ['t0'], ['d2']], ['.', '.', '.']]
    y = [[['d19'], ['t1'], ['d2']], ['.', '.', '.']]
    self.assertFalse(verifier.patternMatch(x,y))
    self.assertFalse(verifier.pM2(x,y))

    x = [['?', ['t0'], ['d2']], ['.', '(', '(']]
    y = [[['d19'], ['t0'], ['d2']], ['.', '.', '.']]
    self.assertFalse(verifier.patternMatch(x,y))
    self.assertFalse(verifier.pM2(x,y))

    cx = [
        ['?', ['t2'], ['t3'], '+', ['d4'], ['t5'], ['d11'], ['t0'], ['d12'], ['t6'], '+', 
          ['t5', '*'], ['d4', '*'], ['t3', '*'], ['t2', '*'], ['d1', '*'], ['t0', '*']], 
        ['.', '(', '(', '+', '(', '(', '.', '.', '.', '.', '+', ')', ')', ')', ')', '.', '.']]
    cy = [
        [['d1'], ['t2'], ['t3'], '+', ['d4'], ['t5'], ['d11'], ['t0'], ['d12'], ['t6'], '+', 
          ['t5', '*'], ['d4', '*'], ['t3', '*'], ['t2', '*'], ['d1', '*'], ['t0', '*']], 
        ['.', '(', '(', '+', '(', '(', '.', '.', '.', '.', '+', ')', ')', ')', ')', '.', '.']]

    self.assertTrue(verifier.patternMatch(cx,cy))
    self.assertTrue(verifier.pM2(cx,cy))

