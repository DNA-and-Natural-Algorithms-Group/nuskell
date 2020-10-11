#!/usr/bin/env python
#
#  tests/dsdcompiler/test_crn_parser.py
#  NuskellCompilerProject
#
import unittest
from pyparsing import ParseException
from nuskell.dsdcompiler.crn_parser import parse_crn_string, Reaction

class TestCRNparser(unittest.TestCase):
    def test_parse_reaction_string(self):
        """ Testing crn string parser """
        crn = "A+B->X+Y"
        pcrn = [Reaction(['A', 'B'], ['X', 'Y'], 1, 0)]
        pfs = {'A': (None, None), 'X': (None, None), 'B': (None, None), 'Y': (None, None)}
        ecrn, efs = parse_crn_string(crn)

        self.assertEqual(ecrn, pcrn, 'single chemical reaction 1a')
        self.assertEqual(sorted(efs), sorted(pfs), 'single chemical reaction 1b')

        crn = "A<=>B"
        pcrn = [Reaction(['A'], ['B'], 1, 1)]
        pfs = {'A': (None, None), 'B': (None, None)}
        ecrn, efs = parse_crn_string(crn)
        self.assertEqual(ecrn, pcrn, 'single chemical reaction 2a')
        self.assertEqual(sorted(efs), sorted(pfs), 'single chemical reaction 2b')

        crn = "->X"
        pcrn = [Reaction([], ['X'], 1, 0)]
        pfs = {'X': (None, None)}
        ecrn, efs = parse_crn_string(crn)
        self.assertEqual(ecrn, pcrn, 'single chemical reaction 3a')
        self.assertEqual(sorted(efs), sorted(pfs), 'single chemical reaction 3b')

        crn = "B->"
        pcrn = [Reaction(['B'], [], 1, 0)]
        pfs = {'B': (None, None)}
        ecrn, efs = parse_crn_string(crn)
        self.assertEqual(ecrn, pcrn, 'single chemical reaction 4a')
        self.assertEqual(sorted(efs), sorted(pfs), 'single chemical reaction 4b')

    def test_parse_crn_string_oneline(self):
        crn = "A+C->X+Y; Y<=>L; L -> C+A"
        ecrn = [Reaction(['A', 'C'], ['X', 'Y'], 1, 0), 
                Reaction(['Y'], ['L'], 1, 1), 
                Reaction(['L'], ['C', 'A'], 1, 0)]
        efs = {'A': (None, None), 
               'C': (None, None), 
               'X': (None, None), 
               'Y': (None, None), 
               'L': (None, None)}
        pcrn, pfs = parse_crn_string(crn)
        self.assertEqual(pcrn, ecrn)
        self.assertDictEqual(pfs, efs)

    def test_deprecated_keywords(self):
        with self.assertRaises(ParseException):
            crn = "A+C->X+Y; formals={A1, A2, B2}"
            pcrn, pfs = parse_crn_string(crn)

if __name__ == '__main__':
    unittest.main()
