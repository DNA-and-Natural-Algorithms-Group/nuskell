#!/usr/bin/env python
#
#  test_pathway_equivalence.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

import os
import unittest
from collections import Counter
from nuskell.crnutils import genCRN, assign_crn_species

from nuskell.verifier.basis_finder import my_parse_crn
from nuskell.verifier.basis_finder import BasisFinderError, NoFormalBasisError

import nuskell.verifier.crn_pathway_equivalence as pathway_equivalence
import nuskell.verifier.basis_finder as basis_finder

import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

SKIP = False

def remove_const(crn, const):
    for rxn in crn:
        for x in const:
            while x in rxn[0]:
                rxn[0].remove(x)
            while x in rxn[1]:
                rxn[1].remove(x)
    return crn

@unittest.skipIf(SKIP, "skipping tests")
class BisimulationTests(unittest.TestCase):
    def test_crn_STW2019_text01(self):
        crn1 = "A + B -> C"
        crn2 = """
            A <=> i # note this is different from paper, which says A -> i
            i + B1 <=> j1
            i + B2 <=> j2
            j1 -> C
            j2 -> C
            """
        fcrn, fs = my_parse_crn(crn1)
        icrn, _ = my_parse_crn(crn2)
        i, wastes, _ = assign_crn_species(icrn, fs)

        print(fs)
        print(wastes)

        ifs  = set(['A', 'B1', 'B2', 'C'])#, 'i': ['A']}

        inter = {'A' : ['A'], 'B1': ['B'], 'B2': ['B'], 'C': ['C']}#, 'i': ['A']}

        for r in genCRN(icrn, rates = False, interpretation = inter):
            print(r)

        basis_raw, basis_int = basis_finder.find_basis(icrn, ifs, interpretation = inter)

    def test_crn_JDW2019_F5(self):
        crn1 = "A + B -> C + D"
        crn2 = "A + B <=> C + D"
        crn3 = """
        A + f_ABCD <=> i_A_BCD + f_A
        B + i_A_BCD <=> i_AB_CD + f_B
        i_AB_CD + f_C <=> i_ABC_D + C
        i_ABC_D + f_D <=> i_ABCD_ + D
        i_ABCD_ + fi -> w_ABCD
        """

        fcrn, fs = my_parse_crn(crn2)
        icrn, _ = my_parse_crn(crn3)
        icrn = remove_const(icrn, ['f_ABCD', 'f_A', 'f_B', 'f_C', 'f_D'])

        i, wastes, _ = assign_crn_species(icrn, fs)

        #basis_raw, basis_int = basis_finder.find_basis(icrn, fs | wastes)

        print(fcrn)
        print(icrn)
        print(wastes)

@unittest.skipIf(SKIP, "skipping tests")
class PathwayEquivalenceTests(unittest.TestCase):
    """Pathway decomposition Testing Class:

    Compares *formal* CRNs with *enumerated* CRNs.

    Note: Translation and enumeration are not part of this testing class. The
          correct translation/enumeration has to be checked elsewhere! It is not
          necessary to refer to translation schemes here at all, but it does make
          sense for reproducability.
    """

    def test_STW17_intro(self):
        crn1 = "A+B -> C+D; C+A -> C+C"
        crn2 = "A<=>i; i+B<=>j; i+j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
        crn3 = "A<=>i; i+B<=>j; j<=>C+k; k->D; C+A<=>m+n; m+n->C+C"
        crn4 = "A->i; i+B<=>j; j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
        crn5 = "A<=>i; i+B<=>j; j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
        crn6 = "A+g1<=>i+g2; i+B<=>j+g3; g4+j->C+k+w1; g5+k<=>D+w2; C+A<=>m+n; g6+m+n->C+C+w3"
        crn6f = "A<=>i; i+B<=>j; j->C+k+w1; k<=>D+w2; C+A<=>m+n; m+n->C+C+w3"

        (crn1, fs) = my_parse_crn(crn1)
        (crn2, _) = my_parse_crn(crn2)
        (crn3, _) = my_parse_crn(crn3)
        (crn4, _) = my_parse_crn(crn4)
        (crn5, _) = my_parse_crn(crn5)
        (crn6, _) = my_parse_crn(crn6)
        (crn6f, _) = my_parse_crn(crn6f)

        inter = {'A': 'A', 'B': 'B', 'C': 'C', 'D': 'D'}

        self.assertFalse(pathway_equivalence.test((crn1, fs), (crn2, fs), inter))
        self.assertFalse(pathway_equivalence.test((crn1, fs), (crn3, fs), inter))
        self.assertFalse(pathway_equivalence.test((crn1, fs), (crn4, fs), inter))
        self.assertTrue( pathway_equivalence.test((crn1, fs), (crn5, fs), inter))
        self.assertFalse(pathway_equivalence.test((crn1, fs), (crn6, fs), inter))
        self.assertFalse(pathway_equivalence.test((crn1, fs), (crn6f, fs), inter))

if __name__ == '__main__':
    unittest.main()

