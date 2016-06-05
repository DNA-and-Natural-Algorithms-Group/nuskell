import unittest
from nuskell.parser import parse_crn_string

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

  def test_parse_crn_string(self):
    """ Testing crn string parser """
    print self.cr1, parse_crn_string(self.cr1)
    res1 = (
        [['irreversible', ['A', 'B'], ['X', 'Y']]], 
        ['A', 'X', 'B', 'Y'], [])

    self.assertEqual(res1, parse_crn_string(self.cr1), 
        'single chemical reaction 1')

    res2 = ([['reversible', ['A'], ['X']]], ['A', 'X'], [])

    self.assertEqual(res2, parse_crn_string(self.cr2), 
        'single chemical reaction 2')

    res4 = ([['irreversible', [], ['X']]], ['X'], [])

    self.assertEqual(res4, parse_crn_string(self.cr3), 
        'single chemical reaction 3')

    print self.cr4, parse_crn_string(self.cr4)
    res4 = ([['irreversible', ['B'], []]], ['B'], [])

    self.assertEqual(res4, parse_crn_string(self.cr4), 
        'single chemical reaction 4')

    self.crn1 = self.cr1 + "\n" + self.cr2

    crn1res = ([['irreversible', ['A', 'B'], ['X', 'Y']],
        ['reversible', ['A'], ['X']]], ['A', 'X', 'B', 'Y'], [])
    self.assertEqual(crn1res, parse_crn_string(self.crn1), 
        'chemical reaction network 1')

    # TODO: not sure if this behavior is desired... species 
    # 'Y' remains unassigned
    self.crn2 = self.cr1 + "\n" + self.cr2  + "\n" + self.fs + "\n" + self.cs
    crn2res = ([['irreversible', ['A', 'B'], ['X', 'Y']],
        ['reversible', ['A'], ['X']]], ['A', 'X', 'B'], ['Y'])
    print 'p', parse_crn_string(self.crn2)
    self.assertEqual(crn2res, parse_crn_string(self.crn2), 
        'chemical reaction network 2')


