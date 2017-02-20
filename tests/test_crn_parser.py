import unittest
from nuskell.parser import parse_crn_string, split_reversible_reactions

class TestCRNparser(unittest.TestCase):
  def setUp(self):
    self.cr1="A + B -> X + Y"
    self.cr2="A <=> X"
    self.cr3=" -> X"
    self.cr4="B -> "
    self.fs="formal = {A, B, X} "
    self.cs="constant = {Y} "


  def tearDown(self):
    pass

  def test_split_reversible_reactions(self):
    self.crn1 = self.cr1 + "\n" + self.cr2
    (crn, fs, cs) = parse_crn_string(self.crn1)
    res = [[['A', 'B'], ['X', 'Y']], [['A'], ['X']], [['X'], ['A']]]
    self.assertEqual(split_reversible_reactions(crn), res, 'split irreversible reactions')

  def test_parse_reaction_string(self):
    """ Testing crn string parser """
    res1 = ([
             [['A', 'B'], ['X', 'Y'], 'irreversible']
            ], ['A', 'X', 'B', 'Y'], [])

    self.assertEqual(parse_crn_string(self.cr1), res1,
        'single chemical reaction 1')

    res2 = ([
             [['A'], ['X'], 'reversible']
            ], ['A', 'X'], [])

    self.assertEqual(parse_crn_string(self.cr2), res2,
        'single chemical reaction 2')

    res3 = ([
             [[], ['X'], 'irreversible']
            ], ['X'], [])

    self.assertEqual(parse_crn_string(self.cr3), res3,
        'single chemical reaction 3')

    res4 = ([
             [['B'], [], 'irreversible']
            ], ['B'], [])

    self.assertEqual(parse_crn_string(self.cr4), res4,
        'single chemical reaction 4')

  def test_parse_crn_string(self):
    crn1 = self.cr1 + "\n" + self.cr2
    crn1res = ([
                [['A', 'B'], ['X', 'Y'], 'irreversible'],
                [['A'], ['X'], 'reversible']
               ], ['A', 'X', 'B', 'Y'], [])
    self.assertEqual(parse_crn_string(crn1), crn1res,
        'chemical reaction network 1')

    crn2 = self.cr1 + "\n" + self.cr2  + "\n" + self.fs + "\n" + self.cs
    crn2res = ([
                [['A', 'B'], ['X', 'Y'], 'irreversible'],
                [['A'], ['X'], 'reversible']
               ], ['A', 'X', 'B'], ['Y'])

    self.assertEqual(parse_crn_string(crn2), crn2res,
        'chemical reaction network 2')

  def test_parse_crn_string_oneline(self):
   crn1 = "A+C->X+Y; Y<=>L; L -> C+A"
   out1 = ([
              [['A', 'C'], ['X', 'Y'], 'irreversible'], 
              [['Y'], ['L'], 'reversible'], 
              [['L'], ['C', 'A'], 'irreversible']
             ], ['A', 'X', 'C', 'Y', 'L'], [])
   self.assertEqual(parse_crn_string(crn1), out1,
        'oneline CRN format')

  #NOTE: this default behavior is not necessary, it might change
   crn2 = "A+C->X+Y; Y<=>L; L -> C+A; formal={A, x, y}; constant={L, C}"
   out2 = ([
            [['A', 'C'], ['X', 'Y'], 'irreversible'], 
            [['Y'], ['L'], 'reversible'], 
            [['L'], ['C', 'A'], 'irreversible']
           ], ['A', 'x', 'y'], ['C', 'L'])

   self.assertEqual(parse_crn_string(crn2), out2,
        'oneline CRN format')

if __name__ == '__main__':
  unittest.main()

