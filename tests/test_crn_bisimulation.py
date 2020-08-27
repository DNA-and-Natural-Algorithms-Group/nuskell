#!/usr/bin/env python
#
#  test_bisimulation.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

import os
import unittest
from collections import Counter

from nuskell.crnutils import parse_crn_file, parse_crn_string, split_reversible_reactions
import nuskell.verifier.crn_bisimulation_equivalence as bisimulation

SKIP_SLOW = True

class BisimulationTests(unittest.TestCase):
    """Bisimulation Testing Class:
    Compares *formal* CRNs with *enumerated* CRNs.

    Note: Translation and enumeration are not part of this testing class. The
    correct translation/enumeration has to be checked elsewhere! It is not
    necessary to refer to translation schemes here at all, but it does make
    sense for reproducability.
    """

    def setUp(self):
        pass

    def tearDown(self):
        # clean up even if unittests failed
        pass

    def _parse_crn_file(self, filename):
        # Load data from filename into the proper bisimulation format
        if not os.path.isfile(filename):
            raise Exception(
                "File for unittest is missing: {}".format(filename))
        crn, formal = parse_crn_file(filename)
        crn = split_reversible_reactions(crn)
        return ([[Counter(part) for part in rxn[:2]] for rxn in crn], formal)

    def _parse_modular_crn_file(self, filename, reversible=False):
        if not os.path.isfile(filename):
            raise Exception(
                "File for unittest is missing: {}".format(filename))
        crn, formal = parse_crn_file(filename)
        if reversible:
            crn = [split_reversible_reactions([rxn]) for rxn in crn]
        else:
            crn = [[rxn] for rxn in split_reversible_reactions(crn)]
        return ([[[Counter(part) for part in rxn[:2]] for rxn in module]
                 for module in crn], formal)

    def _parse_crn_string(self, string):
        crn, formal = parse_crn_string(string)
        crn = split_reversible_reactions(crn)
        return ([[Counter(part) for part in rxn[:2]] for rxn in crn], formal)

    @unittest.skipIf(SKIP_SLOW, "skipping slow tests")
    def test_roessler_qian_equivalence(self):
        (fcrn, fs) = self._parse_crn_file('tests/crns/roessler_01.crn')
        (icrn, _) = self._parse_crn_file('tests/crns/icrns/roessler_qian2011_gen.crn')

        v, i = bisimulation.test(fcrn, icrn, fs)
        self.assertTrue(v)

    def test_roessler_modular_equivalence(self):
        (fcrns, fs) = self._parse_modular_crn_file(
            'tests/crns/roessler_01.crn', reversible=True)
        icrns = [0 for rxn in fcrns]
        for i in range(len(icrns)):
            (icrns[i], _) = self._parse_crn_file(
                'tests/crns/icrns/roessler_qian2011_module' + str(i + 1) + '.crn')

        partial = {sp: Counter({sp: 1}) for sp in fs}
        backup = {sp: Counter({sp: 1}) for sp in fs}
        v, i = bisimulation.testModules(fcrns, icrns, fs, partial,
                                        ispCommon=set(fs))
        self.assertTrue(v)
        # TODO: A function that does not say so, should not modify its arguments.
        #self.assertDictEqual(partial, backup)

    def todo_test_cardelliNM_modularity(self):
        # echo "->A; B->" | nuskell --ts cardelli2011_NM.ts --verify-timeout 5
        # --verify pathway bisimulation modular-bisimulation bisimulation
        pass

    def test_roessler_qian_equivalence_with_modinterpretation(self):
        (fcrn, fs) = self._parse_crn_file('tests/crns/roessler_01.crn')
        (icrn, _) = self._parse_crn_file('tests/crns/icrns/roessler_qian2011_gen.crn')
        (fcrns, _) = self._parse_modular_crn_file(
            'tests/crns/roessler_01.crn', reversible=True)
        icrns = [0 for rxn in fcrns]
        for i in range(len(icrns)):
            (icrns[i], _) = self._parse_crn_file(
                'tests/crns/icrns/roessler_qian2011_module' + str(i + 1) + '.crn')

        partial = {sp: Counter({sp: 1}) for sp in fs}
        v, i = bisimulation.testModules(fcrns, icrns, fs, partial,
                                        ispCommon=set(fs))
        self.assertTrue(v)

        v, i = bisimulation.test(fcrn, icrn, fs, interpretation=i)
        self.assertTrue(v)

    @unittest.skipIf(SKIP_SLOW, "skipping slow tests")
    def test_qingdongthesis_solo(self):
        # An example where the choice of the permchecker matters ...
        (fcrn, fs) = self._parse_crn_file('tests/crns/crn6.crn')
        (icrn, _) = self._parse_crn_file('tests/crns/icrns/crn6_qingdong_thesis.crn')

        inter_01 = {'i778': Counter(['Y']),
                    'i575': Counter(['X']),
                    'i599': Counter(['C']),
                    'i2232': Counter(['A']),
                    'i73': Counter(['B'])}

        inter_02 = {'i842': Counter(['Y', 'X', 'A']),
                    'i394': Counter(['X', 'Y', 'X']),
                    'i119': Counter(['X', 'B', 'A']),
                    'i2300': Counter(['A', 'C']),
                    'i778': Counter(['Y']),
                    'i575': Counter(['X']),
                    'i599': Counter(['C']),
                    'i2232': Counter(['A']),
                    'i73': Counter(['B'])}

        if True: # NOTE: These tests complete in less than 10 minutes
            v, _ = bisimulation.test(fcrn, icrn, fs, permissive='whole-graph')
            self.assertTrue(v)

            v, _ = bisimulation.test(fcrn, icrn, fs, permissive='depth-first')
            self.assertTrue(v)

        if False: # NOTE: These might not even finish in an overnight run ...
            v, _ = bisimulation.test(fcrn, icrn, fs, permissive='loop-search')
            self.assertTrue(v)

            v, _ = bisimulation.test(
                fcrn, icrn, fs, interpretation=inter_01, permissive='loop-search')
            self.assertTrue(v)

        if True: # These tests pass quite fast!
            v, _ = bisimulation.test(fcrn, icrn, fs, interpretation=inter_01)
            self.assertTrue(v)

            v, _ = bisimulation.test(
                fcrn, icrn, fs, interpretation=inter_01, permissive='depth-first')
            self.assertTrue(v)

            v, _ = bisimulation.test(
                fcrn, icrn, fs, interpretation=inter_02, permissive='whole-graph')
            self.assertTrue(v)

            v, _ = bisimulation.test(
                fcrn, icrn, fs, interpretation=inter_02, permissive='depth-first')
            self.assertTrue(v)

            v, _ = bisimulation.test(
                fcrn, icrn, fs, interpretation=inter_02, permissive='loop-search')
            self.assertTrue(v)

    def test_example_01(self):
        # A sample test to aggree on a new interface for bisimulation.
        fcrn = "A->B"
        ecrn = "A<=>i19; i19<=>i39+X; i39->i71+i72"

        (fcrn, fs) = self._parse_crn_string(fcrn)
        (ecrn, _) = self._parse_crn_string(ecrn)

        partial = dict()
        partial['A'] = Counter(A=1)
        partial['B'] = Counter(B=1)

        v, i = bisimulation.test(fcrn, ecrn, fs, interpretation=partial,
                                 permissive='loop-search', verbose=False)
        self.assertTrue(v)

        v, i = bisimulation.test(fcrn, ecrn, fs, interpretation=partial,
                                 permissive='whole-graph', verbose=False)
        self.assertTrue(v)

        v, i = bisimulation.test(fcrn, ecrn, fs, interpretation=partial,
                                 permissive='depth-first', permissive_depth=8, verbose=False)
        self.assertTrue(v)

        # A function that does not say so, should not modify its arguments:
        argcheck = dict()
        argcheck['A'] = Counter(A=1)
        argcheck['B'] = Counter(B=1)
        self.assertDictEqual(partial, argcheck)

    @unittest.skipIf(SKIP_SLOW, "skipping slow tests")
    def test_example_02(self):
        # """A simple example of finding a bisimulation from group meeting."""

        fcrn = "A + B -> C + D ; A + C -> B + D"
        icrn = "x1 -> x2 ; x3 + x4 <=> x5 ; x2 -> x6 + x8 ; x5 -> x7 ; " + \
               "x3 <=> x6 ; x9 <=> x10 ; x10 + x4 <=> x1 ; x7 -> x9 + x8"

        (fcrn, fs) = self._parse_crn_string(fcrn)

        (icrn, _) = self._parse_crn_string(icrn)

        v, _ = bisimulation.test(fcrn, icrn, fs)
        self.assertTrue(v)

        # Test wrong partial interpretation
        partial = dict()
        partial['x2'] = Counter(['B', 'D'])
        partial['x3'] = Counter(['C'])

        v, i = bisimulation.test(fcrn, icrn, fs, interpretation=partial,
                                 permissive='loop-search', verbose=False)
        self.assertFalse(v)

        del partial['x3']

        v, _ = bisimulation.test(fcrn, icrn, fs, permissive='loop-search',
                                 interpretation=partial, verbose=False)
        self.assertTrue(v)

    @unittest.skipIf(SKIP_SLOW, "skipping slow tests")
    def test_example_03(self):
        #""" a follow-up on testing the groupmeeting example """

        fcrn = "A + B -> C + D ; A + C -> B + D"
        icrn = "x1 -> x2 ; x3 + x4 <=> x5 ; x2 -> x6 + x8 ; x5 -> x7 ; " + \
               "x3 <=> x6 ; x9 <=> x10 ; x10 + x4 <=> x1 ; x7 -> x9 + x8"

        # First correct interpretation
        inter1 = {'x1': Counter({'A': 1, 'B': 1}),
                  'x2': Counter({'C': 1, 'D': 1}),
                  'x3': Counter({'C': 1}),
                  'x4': Counter({'A': 1}),
                  'x5': Counter({'A': 1, 'C': 1}),
                  'x6': Counter({'C': 1}),
                  'x7': Counter({'B': 1, 'D': 1}),
                  'x8': Counter({'D': 1}),
                  'x9': Counter({'B': 1}),
                  'x10': Counter({'B': 1})}
        pinter1 = {'x7': Counter({'B': 1, 'D': 1})}

        # Second correct interpretation
        inter2 = {'x1': Counter({'A': 1, 'C': 1}),
                  'x2': Counter({'B': 1, 'D': 1}),
                  'x3': Counter({'B': 1}),
                  'x4': Counter({'A': 1}),
                  'x5': Counter({'A': 1, 'B': 1}),
                  'x6': Counter({'B': 1}),
                  'x7': Counter({'C': 1, 'D': 1}),
                  'x8': Counter({'D': 1}),
                  'x9': Counter({'C': 1}),
                  'x10': Counter({'C': 1})}
        pinter2 = {'x7': Counter({'C': 1, 'D': 1})}

        # CRN preprocessing
        (fcrn, fs) = self._parse_crn_string(fcrn)
        (icrn, _) = self._parse_crn_string(icrn)

        # NOTE: Correct behavior
        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter1,
                                  permissive='whole-graph')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter1,
                                  permissive='loop-search')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter1,
                                  permissive='depth-first')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                                  permissive='loop-search')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                                  permissive='whole-graph')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter1,
                                  permissive='depth-first')
        self.assertTrue(v)
        self.assertDictEqual(inter1, i1)

        v, i2 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter2,
                                  permissive='loop-search')
        self.assertTrue(v)
        self.assertDictEqual(inter2, i2)

        v, i2 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter2,
                                  permissive='whole-graph')
        self.assertTrue(v)
        self.assertDictEqual(inter2, i2)

        v, i2 = bisimulation.test(fcrn, icrn, fs, interpretation=pinter2,
                                  permissive='depth-first')
        self.assertTrue(v)
        self.assertDictEqual(inter2, i2)

    def test_example_04(self):
        # Two valid interpretations
        fcrn = "B + B -> B"
        icrn = "B <=> x1; B + x1 -> x2 + x3; x2 -> B + x4"

        (fcrn, fs) = self._parse_crn_string(fcrn)
        (icrn, _) = self._parse_crn_string(icrn)

        ifull1 = {'B': Counter(['B']),
                  'x1': Counter(['B']),
                  'x2': Counter(['B', 'B']),
                  'x3': Counter(),
                  'x4': Counter()}

        ipart1 = {'B': Counter(['B']),
                  'x2': Counter(['B', 'B'])}

        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=ipart1)
        self.assertTrue(v)
        self.assertDictEqual(i1, ifull1)

        ifull2 = {'B': Counter(['B']),
                  'x1': Counter(['B']),
                  'x2': Counter(['B']),
                  'x3': Counter(),
                  'x4': Counter()}

        ipart2 = {'B': Counter(['B']),
                  'x2': Counter(['B'])}

        v, i2 = bisimulation.test(fcrn, icrn, fs, interpretation=ipart2)
        self.assertTrue(v)
        self.assertDictEqual(i2, ifull2)

    def test_example_05(self):
        # Issue fixed: Naming species in certain ways broke bisimulation
        fcrn = "A+C->A+B"

        #icrn = "A <=> x1 + x2; C+x1 <=> x3 + x4; x3 -> A + B + x5"
        icrn = "A <=> x1 + e45; C + x1 <=> x3 + x4; x3 -> A + B + x5"

        (fcrn, fs) = self._parse_crn_string(fcrn)
        (icrn, _) = self._parse_crn_string(icrn)

        inter = {'A': Counter(['A']),
                 'B': Counter(['B']),
                 'C': Counter(['C'])}

        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter)
        self.assertTrue(v)

    def test_example_GC(self):
        # GC schemes do not produce a correct crn bisimulation ...
        fcrn = "A + B <=> X + Y"
        icrn = """ 
                A <=> i22 [kf = 503452, kr = 503452]
                i59 <=> i139 [kf = 503452, kr = 503452]
                i45 -> i351 + i352 [k = 757794]
                i22 + B <=> i45 + i44 [kf = 503452, kr = 503452]
                i44 <=> i60 + i59 [kf = 503452, kr = 503452]
                i60 -> i104 + i105 [k = 757794]
                i139 <=> i227 + X [kf = 503452, kr = 503452]
                i227 <=> i269 + Y [kf = 503452, kr = 503452]
                i269 -> i338 + i339 [k = 757794]
               """

        (fcrn, fs) = self._parse_crn_string(fcrn)
        (icrn, _) = self._parse_crn_string(icrn)

        inter = {'A': Counter(['A']),
                 'B': Counter(['B']),
                 'X': Counter(['X']),
                 'Y': Counter(['Y']),
                 'i22': Counter(['A']),
                 'i44': Counter(['A', 'B']),
                 'i59': Counter(['A', 'B']),
                 'i139': Counter(['A', 'B']),
                 'i227': Counter(['Y']),
                 'i269': Counter(),
                 'i60': Counter(), 
                 'i104': Counter(),
                 'i105': Counter(),
                 'i45': Counter(), 
                 'i351': Counter(),
                 'i352': Counter()}

        v, i1 = bisimulation.test(fcrn, icrn, fs, interpretation=inter)

        assert v is False

if __name__ == '__main__':
    unittest.main()
