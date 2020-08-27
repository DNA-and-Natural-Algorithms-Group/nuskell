#!/usr/bin/env python
#
#  test_crnutils.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

import unittest

from nuskell.crnutils import Reaction
from nuskell.crnutils import parse_crn_string, split_reversible_reactions
from nuskell.crnutils import crn_to_standard, assign_crn_species


class TestCRN_utils(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_crn_to_standard_format1(self):
        """ Testing CRN format 1"""
        crn = [[['A', 'A', 'B'], ['B', 'C', 'C', 'C']],
               [['X'], ['A', 'C']]]

        new = crn_to_standard(crn)

        assert new[0] == Reaction(['A', 'A', 'B'], ['B', 'C', 'C', 'C'], 1, 0)
        assert new[1] == Reaction(['X'], ['A', 'C'], 1, 0)

    def test_crn_to_standard_format2(self):
        """ Testing CRN format 1"""
        crn = [[['A', 'A', 'B'], ['B', 'C', 'C', 'C'], None],
               [['X'], ['A', 'C'], 3.5]]

        new = crn_to_standard(crn)

        assert new[0] == Reaction(['A', 'A', 'B'], ['B', 'C', 'C', 'C'], 1, 0)
        assert new[1] == Reaction(['X'], ['A', 'C'], 3.5, 0)

    def test_crn_to_standard_format3(self):
        """ Testing CRN format 1"""
        crn = [[['A', 'A', 'B'], ['B', 'C', 'C', 'C'], [None, None]],
               [['X'], ['A', 'C'], [None]]]

        new = crn_to_standard(crn)

        assert new[0] == Reaction(['A', 'A', 'B'], ['B', 'C', 'C', 'C'], 1, 1)
        assert new[1] == Reaction(['X'], ['A', 'C'], 1, 0)

    def test_crn_to_standard_format4(self):
        """ Testing CRN format 1"""
        crn = [[['A', 'A', 'B'], ['B', 'C', 'C', 'C'], [5, 6.5]],
               [['X'], ['A', 'C'], [7, 0]]]

        new = crn_to_standard(crn)

        assert new[0] == Reaction(['A', 'A', 'B'], ['B', 'C', 'C', 'C'], 5, 6.5)
        assert new[1] == Reaction(['X'], ['A', 'C'], 7, 0)

    def test_crn_to_standard_current(self):
        crn = [Reaction(['A', 'A', 'B'], ['B', 'C', 'C', 'C'], 5, 6.5),
               Reaction(['X'], ['A', 'C'], 7, 0)]
        new = crn_to_standard(crn)
        assert new[0] == Reaction(['A', 'A', 'B'], ['B', 'C', 'C', 'C'], 5, 6.5)
        assert new[1] == Reaction(['X'], ['A', 'C'], 7, 0)

    def test_crn_to_standard_invalid(self):
        crn = [[['A', 'A', 'B'], ['B', 'C', 'C', 'C'], [0, 0]]]
        with self.assertRaises(AssertionError) as e:
            new = crn_to_standard(crn)
        crn = [[['A', 'A', 'B'], ['B', 'C', 'C', 'C'], [0, None]]]
        with self.assertRaises(AssertionError) as e:
            new = crn_to_standard(crn)
        crn = [[['A', 'A', 'B'], ['B', 'C', 'C', 'C'], 0]]
        with self.assertRaises(AssertionError) as e:
            new = crn_to_standard(crn)


    def test_assign_species(self):
        crn = "A + B -> C + D"
        (crn, fs) = parse_crn_string(crn)
        fs = set(fs.keys())

        crn = "A + B -> B + w1; w1 -> w1 + w2"
        (crn, _) = parse_crn_string(crn)
        I, W, Wn = assign_crn_species(crn, fs)

        assert I == set()
        assert W == set(['w1', 'w2'])
        assert Wn == set(['w1'])

        crn = """ A + B -> i1 + C + w1; i1 -> D + w2; w1 + w2 -> w3 """
        (crn, _) = parse_crn_string(crn)
        I, W, Wn = assign_crn_species(crn, fs)

        assert I == set(['i1'])
        assert W == set(['w1', 'w2', 'w3'])
        assert Wn == set(['w1', 'w2'])
