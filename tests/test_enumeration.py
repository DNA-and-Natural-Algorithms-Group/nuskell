#!/usr/bin/env python
#
#  test_enumeration.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

import unittest
import argparse

from peppercornenumerator import Enumerator, PepperCondensation
from peppercornenumerator.objects import PepperDomain, PepperComplex, PepperReaction, clear_memory

class EnumerationTests(unittest.TestCase):
    """ Test the results of the peppercorn enumerator.

    As the peppercorn enumerator is developed on its own, this testing class
    ensures that known examples hold and the input/ouput formats remained the
    same after every update of peppercorn.
    """

    def setUp(self):
        parser = argparse.ArgumentParser()
        self.args = parser.parse_args([])

        self.args.verbose = 0

        # Use default behavior
        self.args.max_complex_size = 50
        self.args.max_complex_count = 1000
        self.args.max_reaction_count = 5000

        # Reduce enumerator output
        self.args.reject_remote = False
        self.args.ignore_branch_3way = False
        self.args.ignore_branch_4way = False

        # Use default behavior
        self.args.release_cutoff_1_1 = 6
        self.args.release_cutoff_1_N = 6
        self.args.release_cutoff = None

        self.args.k_slow = 0.0
        self.args.k_fast = 0.0

        self.args.max_helix = True

    def tearDown(self):
        clear_memory()

    def test_peppercorn_interface(self):
        """ Make sure "regular" peppercornenumerator interface did not change """

        # corresponding pil-kernel notation
        # length t0 = 5
        # length d1 = 15
        #
        # c1 = t0 d1
        # c2 = d1( + ) t0*
        # c3 = t0( d1 + d1( + ) )
        # c4 = t0( d1( + ) )
        # c5 = d1

        t0 = PepperDomain('t0', length=5)
        t0_ = PepperDomain('t0*', length=5)
        d1 = PepperDomain('d1', length=15)
        d1_ = PepperDomain('d1*', length=15)
        domains = [t0, t0_, d1, d1_]

        c1 = PepperComplex([t0, d1], ['.','.'], name='c1')
        c2 = PepperComplex([d1, '+', d1_, t0_], ['(', '+', ')', '.'], name='c2')
        c3 = PepperComplex([t0, d1, '+', d1, '+', d1_, t0_], ['(', '.', '+', '(', '+', ')', ')'], name='c3')
        c4 = PepperComplex([t0, d1, '+', d1_, t0_], ['(', '(', '+', ')', ')'], name='c4')
        c5 = PepperComplex([d1], ['.'], name='c5')
        complexes = [c1, c2]

        enum = Enumerator(complexes)

        enum.enumerate()

        # Get full output CRN
        reactions = enum.reactions

        self.assertEqual(len(reactions), 3)

        r1 = PepperReaction([c1, c2], [c3], 'bind21', memorycheck=False)
        r2 = PepperReaction([c3], [c1, c2], 'open', memorycheck=False)
        r3 = PepperReaction([c3], [c4, c5], 'branch-3way', memorycheck=False)

        self.assertEqual(sorted(enum.reactions), sorted([r1, r2, r3]))

        ###########################
        # Get condensed output CRN
        enum.condense()
        reactions = enum.condensed_reactions

        self.assertEqual(len(reactions), 1)
        rC = PepperReaction([c1, c2], [c4, c5], 'condensed', memorycheck=False)
        self.assertEqual(reactions, [rC])


if __name__ == '__main__':
    unittest.main()
