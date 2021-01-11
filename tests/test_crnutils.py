#!/usr/bin/env python
#
#  test_crnutils.py
#  NuskellCompilerProject
#
import unittest

from nuskell.crnutils import (Reaction, 
                              split_reversible_rxns, 
                              combine_reversible_rxns,
                              remove_species,
                              remove_trivial_rxns,
                              remove_duplicate_rxns)

class TestCRN_utils(unittest.TestCase):
    def test_split_reversible_rxns(self):
        crn = [Reaction(['A', 'B'], ['B', 'B'], 1, 1)]
        crn2 = [Reaction(['A', 'B'], ['B', 'B'], 1, 0), 
                Reaction(['B', 'B'], ['A', 'B'], 1, 0)]
        assert split_reversible_rxns(crn) == crn2
        assert split_reversible_rxns(crn2) == crn2

    def test_combine_reversible_rxns(self):
        crn = [Reaction(['A', 'B'], ['B', 'B'], 1, 0), 
               Reaction(['B', 'B'], ['A', 'B'], 1, 0)]
        crn2 = [Reaction(['A', 'B'], ['B', 'B'], 1, 1)]
        assert combine_reversible_rxns(crn) == crn2

    def test_remove_species(self):
        crn = [Reaction(['A', 'B'], ['B', 'B'], 1, 1)]
        crn2 = [Reaction(['A'], [], 1, 1)]
        assert remove_species(crn, ['B']) == crn2

    def test_remove_trivial_rxns(self):
        crn = [Reaction(['A', 'B'], ['B', 'B'], 1, 1),
               Reaction(['B'], ['B'], 1, 0)]
        crn2 = [Reaction(['A', 'B'], ['B', 'B'], 1, 1)]
        assert remove_trivial_rxns(crn) == crn2

    def test_remove_duplicate_rxns(self):
        crn = [Reaction(['A', 'B'], ['B', 'B'], 1, 1),
               Reaction(['A', 'B'], ['B', 'B'], 1, 0)]
        crn2 = [Reaction(['A', 'B'], ['B', 'B'], 1, 1)]
        assert remove_duplicate_rxns(crn) == crn2
