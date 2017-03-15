import unittest
from nuskell.parser import parse_crn_string, split_reversible_reactions

class TestCRNparser(unittest.TestCase):
  def setUp(self):
    pass

  def tearDown(self):
    pass

  def test_parse_reaction_string(self):
    """ Testing crn string parser """

    crn = "A+B->X+Y"
    res = ([[['A', 'B'], ['X', 'Y'], [None]]], ['A', 'B', 'X', 'Y'], [])
    self.assertEqual(parse_crn_string(crn), res, 'single chemical reaction 1')

    crn = "A<=>B"
    res = ([[['A'], ['B'], [None,None]]], ['A', 'B'], [])
    self.assertEqual(parse_crn_string(crn), res, 'single chemical reaction 2')

    crn = "->X"
    res = ([[[], ['X'], [None]]], ['X'], [])
    self.assertEqual(parse_crn_string(crn), res, 'single chemical reaction 3')

    crn = "B->"
    res = ([[['B'], [], [None]]], ['B'], [])
    self.assertEqual(parse_crn_string(crn), res, 'single chemical reaction 4')

  def test_split_reversible_reactions(self):
    crn = "A+B->X+Y\nA<=>X"
    (crn, fs, cs) = parse_crn_string(crn)
    res = [[['A', 'B'], ['X', 'Y'], [None]], [['A'], ['X'], [None]], [['X'], ['A'],[None]]]
    self.assertEqual(split_reversible_reactions(crn), res, 'split irreversible reactions')

  def test_parse_crn_string(self):
    crn = "A+B->X+Y\nA<=>X"
    res = ([[['A', 'B'], ['X', 'Y'], [None]], 
            [['A'], ['X'], [None, None]]], 
           ['A', 'B', 'X', 'Y'], [])

    self.assertEqual(parse_crn_string(crn), res,
        'chemical reaction network 1')

    crn =  "A+B->X+Y;A<=>X;formal={A,B,X};constant={Y}"
    res = ([
                [['A', 'B'], ['X', 'Y'], [None]],
                [['A'], ['X'], [None, None]]
               ], ['A', 'B', 'X'], ['Y'])
    self.assertEqual(parse_crn_string(crn), res,
        'chemical reaction network 2')

  def test_parse_crn_string_oneline(self):
   crn = "A+C->X+Y; Y<=>L; L -> C+A"
   out = ([[['A', 'C'], ['X', 'Y'], [None]], 
           [['Y'], ['L'], [None, None]], 
           [['L'], ['C', 'A'], [None]]
          ], ['A', 'C', 'L', 'X', 'Y'], [])
   self.assertEqual(parse_crn_string(crn), out,
        'oneline CRN format')

   crn = "A+C->X+Y; formal={A1, A2, B2}"
   out = ([[['A', 'C'], ['X', 'Y'], [None]]],
          ['A1', 'A2', 'B2'], ['A','C','X','Y'])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

   crn = "A+C->X+X; formal={A1, A2, B2}"
   out = ([[['A', 'C'], ['X', 'X'], [None]]],
          ['A1', 'A2', 'B2'], ['A','C','X'])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

   crn = "A+C->2X; formal={A1, A2, B2}"
   out = ([[['A', 'C'], ['X', 'X'], [None]]],
          ['A1', 'A2', 'B2'], ['A','C','X'])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

   crn = "A+C->X+Y; formal={A1, A2, B2, X, Y}; constant={}"
   out = ([[['A', 'C'], ['X', 'Y'], [None]]],
          ['A1', 'A2', 'B2', 'X', 'Y'], [])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

   crn = "A+C->X+Y; Y<=>L; L->C+A; formal={A, x, y}; constant={L, C}"
   out = ([
            [['A', 'C'], ['X', 'Y'], [None]], 
            [['Y'], ['L'], [None, None]], 
            [['L'], ['C', 'A'], [None]]
           ], ['A', 'x', 'y'], ['C', 'L'])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

if __name__ == '__main__':
  unittest.main()

