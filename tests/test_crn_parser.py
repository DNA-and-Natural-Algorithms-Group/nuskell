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
    res = ([[['A', 'B'], ['X', 'Y'], [None]]], ['A', 'B', 'X', 'Y'], ['A', 'B', 'X', 'Y'], [])
    self.assertEqual(parse_crn_string(crn), res, 'single chemical reaction 1')

    crn = "A<=>B"
    res = ([[['A'], ['B'], [None,None]]], ['A', 'B'], ['A', 'B'], [])
    self.assertEqual(parse_crn_string(crn), res, 'single chemical reaction 2')

    crn = "->X"
    res = ([[[], ['X'], [None]]], ['X'], ['X'], [])
    self.assertEqual(parse_crn_string(crn), res, 'single chemical reaction 3')

    crn = "B->"
    res = ([[['B'], [], [None]]], ['B'], ['B'], [])
    self.assertEqual(parse_crn_string(crn), res, 'single chemical reaction 4')

  def test_split_reversible_reactions(self):
    crn = "A+B->X+Y\nA<=>X"
    crn, fs, signals, fuels = parse_crn_string(crn)
    res = [[['A', 'B'], ['X', 'Y'], [None]], [['A'], ['X'], [None]], [['X'], ['A'],[None]]]
    self.assertEqual(split_reversible_reactions(crn), res, 'split irreversible reactions')

  def test_parse_crn_string(self):
    crn = "A+B->X+Y\nA<=>X"
    res = ([[['A', 'B'], ['X', 'Y'], [None]], 
            [['A'], ['X'], [None, None]]], 
           ['A', 'B', 'X', 'Y'],
           ['A', 'B', 'X', 'Y'], [])

    self.assertEqual(parse_crn_string(crn), res,
        'chemical reaction network 1')

    crn =  "A+B->X+Y;A<=>X;formals={A,B,X};signals={Y}"
    res = ([[['A', 'B'], ['X', 'Y'], [None]],
            [['A'], ['X'], [None, None]]], 
            ['A', 'B', 'X', 'Y'], ['Y'], [])
    self.assertEqual(parse_crn_string(crn), res,
        'chemical reaction network 2')

  def test_parse_crn_string_oneline(self):
   crn = "A+C->X+Y; Y<=>L; L -> C+A"
   out = ([[['A', 'C'], ['X', 'Y'], [None]], 
           [['Y'], ['L'], [None, None]], 
           [['L'], ['C', 'A'], [None]]
          ], ['A', 'C', 'L', 'X', 'Y'], 
          ['A', 'C', 'L', 'X', 'Y'], [])
   self.assertEqual(parse_crn_string(crn), out,
        'oneline CRN format')

   crn = "A+C->X+Y; formals={A1, A2, B2}"
   out = ([[['A', 'C'], ['X', 'Y'], [None]]],
          ['A', 'A1', 'A2', 'B2', 'C','X','Y'],
          ['A', 'A1', 'A2', 'B2', 'C','X','Y'], [])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

   crn = "A+C->X+X; formals={A1, A2, B2}"
   out = ([[['A', 'C'], ['X', 'X'], [None]]],
          ['A', 'A1', 'A2', 'B2', 'C','X'],
          ['A', 'A1', 'A2', 'B2', 'C','X'], [])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

   crn = "A+C->2X; formals={A1, A2, B2}"
   out = ([[['A', 'C'], ['X', 'X'], [None]]],
          ['A', 'A1', 'A2', 'B2', 'C','X'],
          ['A', 'A1', 'A2', 'B2', 'C','X'], [])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

   crn = "A+C->X+Y; formals={A1, A2, B2, X, Y}; fuels={}"
   out = ([[['A', 'C'], ['X', 'Y'], [None]]],
          ['A', 'A1', 'A2', 'B2', 'C', 'X', 'Y'],
          ['A', 'A1', 'A2', 'B2', 'C', 'X', 'Y'], [])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

   crn = "A+C->X+Y; Y<=>L; L->C+A; formals={A, x, y}; signals={L, C}"
   out = ([
            [['A', 'C'], ['X', 'Y'], [None]], 
            [['Y'], ['L'], [None, None]], 
            [['L'], ['C', 'A'], [None]]
           ], ['A', 'C', 'L', 'X', 'Y', 'x', 'y'], ['C', 'L'], [])
   self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

if __name__ == '__main__':
  unittest.main()

