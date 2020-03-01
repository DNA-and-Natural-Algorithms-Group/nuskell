#!/usr/bin/env python
#
#  test_crn_parser.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

import unittest
from nuskell.crnutils import parse_crn_string, split_reversible_reactions, Reaction


class TestCRNparser(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

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

    def test_split_reversible_reactions(self):
        crn = "A+B->X+Y\nA<=>X"
        pcrn = [Reaction(['A', 'B'], ['X', 'Y'], 1, 0), Reaction(['A'], ['X'], 1, 1)]
        pfs = {'A': (None, None), 'B': (None, None), 'X': (None, None), 'Y': (None, None)}
        ecrn, efs = parse_crn_string(crn)
        self.assertEqual(ecrn, pcrn, 'single chemical reaction 5a')
        self.assertEqual(sorted(efs), sorted(pfs), 'single chemical reaction 5b')
        splitcrn = [Reaction(['A', 'B'], ['X', 'Y'], 1, 0), 
                    Reaction(['A'], ['X'], 1, 0),
                    Reaction(['X'], ['A'], 1, 0)]
        self.assertEqual(sorted(split_reversible_reactions(ecrn)), 
                         sorted(splitcrn), 'split irreversible reactions')

    def todo_test_parse_crn_string(self):
        # This does not work atm
        crn = "A+B->X+Y;A<=>X;formals={A,B,X};signals={Y}"
        res = ([[['A', 'B'], ['X', 'Y'], [1]], [['A'], ['X'], [1, 1]]], ['A', 'B', 'X', 'Y'])
        self.assertEqual(parse_crn_string(crn), res, 'chemical reaction network 2')

    def todo_test_parse_crn_string_oneline(self):
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
               ['A', 'A1', 'A2', 'B2', 'C', 'X', 'Y'],
               ['A', 'A1', 'A2', 'B2', 'C', 'X', 'Y'], [])
        self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

        crn = "A+C->X+X; formals={A1, A2, B2}"
        out = ([[['A', 'C'], ['X', 'X'], [None]]],
               ['A', 'A1', 'A2', 'B2', 'C', 'X'],
               ['A', 'A1', 'A2', 'B2', 'C', 'X'], [])
        self.assertEqual(parse_crn_string(crn), out, 'oneline CRN format')

        crn = "A+C->2X; formals={A1, A2, B2}"
        out = ([[['A', 'C'], ['X', 'X'], [None]]],
               ['A', 'A1', 'A2', 'B2', 'C', 'X'],
               ['A', 'A1', 'A2', 'B2', 'C', 'X'], [])
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
