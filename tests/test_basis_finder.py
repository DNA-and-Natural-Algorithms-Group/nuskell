#!/usr/bin/env python
#
#  test_basis_finder.py
#  NuskellCompilerProject
#
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

import unittest
from collections import Counter
from nuskell.crnutils import parse_crn_string
from nuskell.crnutils import split_reversible_reactions

from nuskell.verifier.basis_finder import Path, find_basis, my_parse_crn
from nuskell.verifier.basis_finder import BasisFinderError, NoFormalBasisError

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestPathProperties(unittest.TestCase):
    def test_wmax_examples(self):
        fs = set(['A', 'B', 'C', 'X', 'Y', 'Z'])

        w, b, iR, fR = 0, 2, 0, 1
        wmax = w * iR + b
        assert wmax == 2
        p1 = 'C -> A + i'
        path0 = Path([], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path1
        assert pathC.is_prime
        assert pathC.width == 2 <= wmax

        w, b, iR, fR = 2, 4, 3, 1
        assert iR + fR == b
        wmax = w * iR + b
        assert wmax == 10 
        p0 = 'A -> A + i'
        p1 = 'A + 3i -> 3i + B'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path0 + path0 + path1
        assert pathC.is_prime
        assert pathC.width == 4 <= wmax
        pathC = path0 + path0 + path0 + path0 + path1
        assert not pathC.is_prime

        w, b, iR, fR = 2, 4, 3, 1
        assert iR + fR == b
        wmax = w * iR + b
        assert wmax == 10 
        p0 = 'C -> A + i'
        p1 = 'A + 3i -> 3i + B'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path0 + path0 + path1
        assert pathC.is_prime
        assert pathC.width == 6 <= wmax
        pathC = path0 + path0 + path0 + path0 + path1
        assert not pathC.is_prime

        w, b, iR, fR = 1, 2, 1, 1
        wmax = w * iR + b
        assert wmax == 3
        p0 = 'A -> i'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.width == w
        p1 = 'X + i -> l'
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path1
        assert pathC.is_prime 
        assert pathC.width == 2 <= wmax
        rxn1 = 'i -> l + k'
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path1
        assert pathC.is_prime 
        assert pathC.width == 2 <= wmax

        w, b, iR, fR = 3, 4, 4, None
        wmax = w * iR + b
        assert wmax == 16
        p0 = 'A -> i; B -> j; B -> j; 2j + i -> k'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.width == w
        p1 = '4k -> C'
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path0 + path0 + path0 + path1
        assert pathC.is_prime 
        assert pathC.width == 12 <= wmax 
        p1 = '3X + k -> C'
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path0 + path0 + path0 + path1
        assert not pathC.is_prime 

    def test_imax_examples(self):
        fs = set(['A', 'B', 'C', 'X', 'Y', 'Z'])

        i, b, iR, fR = 1, 4, 3, 1 
        assert iR + fR == b
        imax = i * iR + fR
        assert imax == 4
        p0 = 'A -> A + i'
        p1 = 'A + 3i -> 3i + B'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path0 + path0 + path1
        assert pathC.is_prime
        assert len(pathC.S0) == 1 <= imax
        pathC = path0 + path0 + path0 + path0 + path1
        assert not pathC.is_prime

        i, b, iR, fR = 1, 4, 3, 1 
        assert iR + fR == b
        imax = i * iR + fR
        assert imax == 4
        p0 = 'C -> A + i'
        p1 = 'A + 3i -> 3i + B'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path0 + path0 + path1
        assert pathC.is_prime
        assert len(pathC.S0) == 3 <= imax
        pathC = path0 + path0 + path0 + path0 + path1
        assert not pathC.is_prime

        i, b, iR, fR = 1, 4, 3, 1 
        assert iR + fR == b
        imax = i * iR + fR
        assert imax == 4
        p0 = 'C -> A + i'
        p1 = 'X + 3i -> 3i + Y'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path0 + path0 + path1
        assert pathC.is_prime
        assert len(pathC.S0) == 4 <= imax
        pathC = path0 + path0 + path0 + path0 + path1
        assert not pathC.is_prime

        i, b, iR, fR = 2, 5, 3, 2 
        assert iR + fR == b
        imax = i * iR + fR
        assert imax == 8
        p0 = 'B + C -> A + i'
        p1 = 'X + Y + 3i -> 3i + Z'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        pathC = path0 + path0 + path0 + path1
        assert pathC.is_prime
        assert len(pathC.S0) == 8 <= imax
        pathC = path0 + path0 + path0 + path0 + path1
        assert not pathC.is_prime

    def test_utilities(self):
        fs = ['A', 'B', 'X']
        lpath = [[['A'],['i']],
                 [['i'],['j']],
                 [['j'],['i']],
                 [['i'],['k', 'l']],
                 [['k'],['B']],
                 [['l'],['A']]]
        path = Path(lpath, fs)
        hash(path)
        assert path.is_formal
        assert path.is_semiformal
 
    def test_errors(self):
        fs = ['A', 'B', 'C', 'X', 'Y']
        lpath = [[['A'],['i']],
                 [['i'],['j']],
                 [['j'],['i']],
                 [['i'],['k', 'l']],
                 [['k'],['B']],
                 [['l'],['A']]]
        path = Path(lpath)
        with self.assertRaises(BasisFinderError) as e:
            path.formal_closure

        lpath = [[['X'], ['i', 'k', 'l']],
                 [['l'], ['p', 'q']],
                 [['i', 'k', 'l'], ['A']],
                 [['p', 'q'], ['k']],
                 [['k'], ['l']]]
        path = Path(lpath, fs)
        assert path.minimal_initial_state == Counter({'X': 1, 'l': 1}) 
        assert path.final_state == Counter({'A': 1, 'l': 1})
        assert path.width == 5
        assert path.formal_closure == Counter({'X': 1, 'A': 1})
        assert path.rfs == set({('A',)})
        with self.assertRaises(BasisFinderError) as e:
            path.dfs

        lpath = [[['A'],['i']],
                 [['i', 'B', 'B'],['j']],
                 [['X'], ['i', 'k', 'l']],
                 [['l'], ['p', 'q']],
                 [['i', 'B'],['k', 'l', 'X']],
                 [['i', 'k', 'l'], ['A']],
                 [['p', 'q'], ['k']],
                 [['k'], ['l']],
                 [['j'],['i', 'B']],
                 [['l'],['A']]]
        path = Path(lpath, fs)
        RFS = set([('A', 'A', 'B', 'X'),])
        assert RFS == path.rfs
        with self.assertRaises(BasisFinderError) as e:
            path.dfs

        lpath = [[['A'], ['i']],
                [['i'], ['i', 'j']],
                [['i'], ['i', 'j']],
                [['i'], ['i', 'j']],
                [['j'], ['k']],
                [['i'], ['l']],
                [['k', 'l'], ['B']],
                [['i', 'i'], ['C']]]
        path = Path(lpath, fs)
        RFS = set([('B', 'C')])
        assert RFS == path.rfs
        with self.assertRaises(BasisFinderError) as e:
            path.dfs

    def test_is_linear(self):
        fs = set('ABCXY')
        def nonwastes(crn, fs):
            species = set().union(*[set().union(*rxn) for rxn in crn])
            intermediates = species - fs
            reactants = set().union(*[set(r) for [r, p] in crn])
            wastes = intermediates - reactants 
            nonwastes = intermediates - wastes
            return nonwastes
 
        p0 = 'A -> i; i-> j; j->i'
        p0, _ = my_parse_crn(p0)
        path = Path(p0)
        nonw = nonwastes(p0, fs)
        assert set(['i', 'j']) == nonw
        assert path.is_linear(nonw) is True

        p0 = 'A -> i; i-> k + l; k -> B; l -> A'
        p0, _ = my_parse_crn(p0)
        path = Path(p0)
        nonw = nonwastes(p0, fs)
        assert set(['i', 'l', 'k']) == nonw
        assert path.is_linear(nonw) is False

        p0 = 'A -> i; i-> k + l; k -> B'
        p0, _ = my_parse_crn(p0)
        path = Path(p0)
        nonw = nonwastes(p0, fs)
        assert set(['i', 'k']) == nonw
        assert path.is_linear(nonw) is True

        p0 = 'A -> k; i -> k; k + k -> B'
        p0, _ = my_parse_crn(p0)
        path = Path(p0)
        nonw = nonwastes(p0, fs)
        assert set(['i', 'k']) == nonw
        assert path.is_linear(nonw) is False


    def test_equality(self):
        fs = set('ABCXY')

        # same signature
        p0 = 'A -> i'
        p1 = 'A -> i; i-> j; j->i'
        p0, _ = my_parse_crn(p0)
        p1, _ = my_parse_crn(p1)
        assert Path(p0, fs) == Path(p1, fs)

        p0 = 'A -> i29; i29 -> i46 + i47'
        p1 = 'A -> i29; i29 -> i46 + i47; i47 -> i70'
        p2 = 'A -> i29; i29 -> i46 + i47; i47 -> i70; i70 -> i47'
        p0, _ = my_parse_crn(p0)
        p1, _ = my_parse_crn(p1)
        p2, _ = my_parse_crn(p2)
        path0 = Path(p0, fs)
        path1 = Path(p1, fs)
        path2 = Path(p2, fs)
        assert path0 != path1
        assert path0 == path2

        # same signature
        p0 = 'A -> i; B -> j'
        p1 = 'B -> j; A -> i'
        p0, _ = my_parse_crn(p0)
        p1, _ = my_parse_crn(p1)
        assert Path(p0, fs) == Path(p1, fs)

        # same signature
        p0 = 'A -> i; i -> B'
        p1 = 'A -> i; i -> j; j -> B'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        assert path0 == path1

        # different width, different path
        p0 = 'A -> i; i -> i + k; k + i -> B'
        p1 = 'A -> p; p -> k + l; k + l -> B'
        p2 = 'A -> p; p -> k; k -> B'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        path2 = Path(my_parse_crn(p2)[0], fs)
        assert path0.minimal_initial_state == path1.minimal_initial_state \
                == path2.minimal_initial_state
        assert path0.final_state == path1.final_state == path2.final_state
        assert path0.width == path1.width != path2.width
        assert path0.formal_closure == path1.formal_closure == path2.formal_closure
        assert path0.dfs == path1.dfs == path2.dfs
        assert path0.rfs == path1.rfs == path2.rfs
        assert path0 == path1 != path2

        # TODO different formal closure

        # Regular vs. non regular final states.
        p0 = 'A -> i; B -> C + j'#; C + i + j -> k'
        p1 = 'B -> i; A -> C + j'#; C + i + j -> k'
        p2 = 'A -> C + j; B -> i'#; C + i + j -> k'
        path0 = Path(my_parse_crn(p0)[0], fs)
        path1 = Path(my_parse_crn(p1)[0], fs)
        path2 = Path(my_parse_crn(p2)[0], fs)
        assert path0.minimal_initial_state == path1.minimal_initial_state == path2.minimal_initial_state
        assert path0.final_state == path1.final_state == path2.final_state
        assert path0.width == path1.width == path2.width
        assert path0.formal_closure == path1.formal_closure == path2.formal_closure
        assert path0.dfs == path1.dfs == path2.dfs
        assert path0.rfs == path1.rfs != path2.rfs
        assert path0 == path1 != path2

    def test_formal_closure(self):
        fs = set(['A', 'B', 'C', 'X'])
        path = [[['A'],['i']],
                [['i'],['j']],
                [['j'],['i']],
                [['i'],['k', 'l']],
                [['k'],['B']],
                [['l'],['A']]]
        fclosure = Counter({'A': 1, 'B': 1})
        assert fclosure == Path(path, fs).formal_closure

        path = [[['A'],['i']],
                [['i'],['j']],
                [['j'],['i']],
                [['i'],['k', 'l']],
                [['l'],['A']]]
        fclosure = Counter(['A'])
        assert fclosure == Path(path, fs).formal_closure

        path = [[['A'],['i']],
                [['i', 'B', 'B'],['j']],
                [['j'],['i', 'B']],
                [['i', 'B'],['k', 'l', 'X']],
                [['l'],['A']]]
        fclosure = Counter({'A': 1, 'B': 2, 'X': 1})
        assert fclosure == Path(path, fs).formal_closure

        path = [[['A'],['i']],
                [['i'],['A', 'j']],
                [['A', 'j'], ['B']]]
        fclosure = Counter({'A': 1, 'B': 1})
        assert fclosure == Path(path, fs).formal_closure

    def test_regular_final_states(self):
        """ STW2019 - RFS examples. """
        fs = set('ABCXY')
        pX = 'A -> i; B + i -> j; j -> B + i; i + B -> C'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({('C', ), ('B', 'C')})

        pX = 'A -> i; B + i -> j; j -> X + k'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({('X',)})

        pX = 'A -> i; i -> B + j; B + j -> k'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({('B',)})

        pX = 'A -> m; m -> i; i -> A + j; j -> p; p -> j; A + j -> k; k -> l; l -> B'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({('A', 'B'), ('B',)})

        pX = 'A -> i; i -> A + j; j -> p; p -> j; A + j -> B'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({('A', 'B'), ('B',)})

        pX = 'A -> i; B + i -> j; j -> B + i; B + i -> j; j -> C + k; C + k -> j; j -> B + i; i -> A'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({('A', 'B', 'C')})

        # Temporarily non-regular (1/2)
        p0 = 'A -> i; B -> j; i + j -> k; k -> X + l; X -> m'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.rfs == frozenset({('X',)})
        # Temporarily non-regular (2/2)
        p0 = 'A -> i; B -> j; i + j -> k; k -> X + l; X -> m; m + l -> X + Y'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.rfs == frozenset({('X','Y')})

        # Two turning points, same result ...
        p0 = 'B -> i; A + i -> j + A; j -> k; A + k -> A + C'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.rfs == frozenset({('A','C')})

        # No turning point ...
        p0 = 'X -> i; B -> B + j; A -> A + p; i + j + p + B -> B + Y'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.rfs == frozenset()

    def dont_test_extended_path_properties(self):
        # NOTE: Potential extensions for the path object. 
        # Most of these properties do not exist, but 
        # they may become useful ...
        fs = set('ABCXYZ')

        p0 = 'i -> j + B; j -> C'
        path = Path(my_parse_crn(p0)[0], fs)
        assert path.is_prime
        assert not path.is_formal
        assert path.is_intermediate
        assert path.is_closing

        p0 = 'i -> j + B; j + B -> C'
        path = Path(my_parse_crn(p0)[0], fs)
        assert path.is_prime
        assert not path.is_formal
        assert path.is_intermediate
        assert not path.is_closing

        p0 = 'i -> k + B; j -> l; k + l -> C'
        path = Path(my_parse_crn(p0)[0], fs)
        assert path.is_prime
        assert not path.is_formal
        assert path.is_intermediate
        assert path.is_closing

        p0 = 'i -> k + B; j -> l'
        path = Path(my_parse_crn(p0)[0], fs)
        assert not path.is_prime
        assert not path.is_formal
        assert path.is_intermediate
        assert not path.is_closing

        p0 = 'i -> j + B; j -> l'
        path = Path(my_parse_crn(p0)[0], fs)
        assert not path.is_prime
        assert not path.is_formal
        assert path.is_intermediate
        assert not path.is_closing

    def test_decomposable_final_states(self):
        fs = set('AVWXY')

        p1 = 'A -> i29'
        # ((i29),()) => ()
        dfs1 = frozenset()
        p2 = 'A -> i29; A -> i29'
        # [[['i29'], ['i29']]]
        # ((i29, i29), ()), ((i29,), (i29,))
        dfs2 = frozenset({('i29',)})
        p3 = 'A -> i29; A -> i29; i29 -> A'
        # [[['A'], ['i29']], [['i29'], ['A']]]
        # ((A, i29), ()), ((A,), (i29,))
        dfs3 = frozenset({('i29',), ('A',)})
        p4 = 'A -> i29; A -> i29; i29 -> A; A -> i29'
        # [[['A'], ['i29', 'i29']], [['A', 'i29'], ['i29']], [['i29'], ['i29']],  ... ]
        # ((i29, i29), ()), ((A, i29), (i29)), ((i29,), (i29,)), ((A,), (i29, i29))
        dfs4 = frozenset({('i29',), ('i29', 'i29'), ('A',), ('A', 'i29')})
        p5 = 'A -> i29; A -> i29; i29 -> A; A -> i29; i29 -> A'
        # [[['A'], ['A', 'i29']], [['A'], ['i29']], [['A', 'A'], ['i29']], ... ]
        # ((A, i29), ()), ((A, A), (i29,)), ((A, i29), (A,)), ((A,), (i29,)), ((A,), (A, i29))
        dfs5 = frozenset({('i29',), ('A', 'A'), ('A',), ('A', 'i29')})
        path1 = Path(my_parse_crn(p1)[0], fs)
        path2 = Path(my_parse_crn(p2)[0], fs)
        path3 = Path(my_parse_crn(p3)[0], fs)
        path4 = Path(my_parse_crn(p4)[0], fs)
        path5 = Path(my_parse_crn(p5)[0], fs)
        assert path1.dfs == dfs1
        assert path2.dfs == dfs2
        assert path3.dfs == dfs3
        assert path4.dfs == dfs4
        assert path5.dfs == dfs5

        pX = 'A -> i29; i29 -> i46 + i47; i47 -> i70; i46 -> X + Y; i70 -> A + i89; A -> i29; i29 -> i46 + i47; i47 -> i70; i46 -> X + Y; i70 -> i47; A + i89 -> i70'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.signature == (('A', 'A'), 
                                   ('X', 'X', 'Y', 'Y', 'i47', 'i70'), 
                                   7, 
                                   ('A', 'A', 'X', 'X', 'Y', 'Y'), 
                                   frozenset({
                                       ('X', 'Y', 'i47'), 
                                       ('X', 'Y', 'i70')}), 
                                   frozenset())

        pX = '-> i3;  -> i3; A -> i7; i3 -> ; i7 -> i16 + i8; i16 -> i17'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.dfs == frozenset({('i17', 'i3', 'i8'), ('i17', 'i8'), ('i3',), ()})

        pX = 'X -> i3; X -> i3; A -> i7; i3 -> X'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.dfs == frozenset({('i7',), ('X',), ('X', 'i3'), ('i3',), ('X', 'i7'), ('i3', 'i7')})

        pX = '-> i3; -> i3; A -> i7; i3 -> '
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.dfs == frozenset({('i3', 'i7'), ('i7',), (), ('i3',)})

    def test_more_pathway_properties(self):
        fs = ['A', 'B', 'C', 'X', 'Y']
        lpath = [[['A'],['i']],
                 [['i'],['j']],
                 [['j'],['i']],
                 [['i'],['k', 'l']],
                 [['k'],['B']],
                 [['l'],['A']]]
        path = Path(lpath, fs)
        assert path.minimal_initial_state == Counter({'A': 1})
        assert path.final_state == Counter({'A': 1, 'B': 1})
        assert path.width == 2
        assert path.formal_closure == Counter({'A': 1, 'B': 1})
        assert path.dfs == set()
        assert path.rfs == frozenset({('A', 'B')})
 
        lpath = [[['A'],['i']],
                 [['i'],['j']],
                 [['j'],['i']],
                 [['i'],['k', 'l']],
                 [['l'],['A']]]
        path = Path(lpath, fs)
        assert path.minimal_initial_state == Counter({'A': 1})
        assert path.final_state == Counter({'A': 1, 'k': 1})
        assert path.width == 2
        assert path.formal_closure == Counter({'A': 1})
        assert path.dfs == set()
        assert path.rfs == frozenset({('A',)})

        lpath = [[['A'],['i']],
                 [['i', 'B', 'B'],['j']],
                 [['j'],['i', 'B']],
                 [['i', 'B'],['k', 'l', 'X']],
                 [['l'],['A']]]
        path = Path(lpath, fs)
        assert path.minimal_initial_state == Counter({'A': 1, 'B': 2})
        assert path.final_state == Counter({'A': 1, 'k': 1, 'X': 1})
        assert path.width == 3
        assert path.formal_closure == Counter({'A': 1, 'B': 2, 'X': 1})
        assert path.dfs == set()
        assert path.rfs == frozenset({('A', 'B', 'X'), ('A', 'X')})

        lpath = [[['A'],['i']],
                 [['i'],['A', 'j']],
                 [['A', 'j'], ['B']]]
        path = Path(lpath, fs)
        assert path.minimal_initial_state == Counter({'A': 1})
        assert path.final_state == Counter({'B': 1})
        assert path.width == 2
        assert path.formal_closure == Counter({'A': 1, 'B': 1})
        assert path.dfs == set()
        assert path.rfs == frozenset({('B',), ('A', 'B')})

        cpath = [[Counter({'B': 1}), Counter({'i11': 1})], 
                 [Counter({'B': 1}), Counter({'i11': 1})]]
        path = Path(cpath, fs)
        assert path.minimal_initial_state == Counter({'B': 2}) 
        assert path.final_state == Counter({'i11': 2})
        assert path.width == 2
        assert path.formal_closure == Counter({'B': 2})
        assert path.dfs == frozenset({('i11',)})
        assert path.rfs == frozenset({()})

        lpath = [[['A'], ['i']],
                 [['i'], ['j', 'j']],
                 [['B'], ['j']],
                 [['j'], ['B', 'C']],
                 [['j', 'j'], ['i']],
                 [['i'], ['C']]]
        path = Path(lpath, fs)
        RFS = set([('B', 'C', 'C')])
        assert RFS == path.rfs
        DFS = set([('C',), ('B', 'C')])
        assert path.dfs == DFS

        lpath = [[['A'], ['m']],
                 [['A'], ['m']],
                 [['A'], ['m']],
                 [['A'], ['m']],
                 [['m'], ['A']]]
        path = Path(lpath, fs)
        DFS = {('m',), ('A', 'm'), ('A', 'm', 'm'), ('m', 'm'), ('m', 'm', 'm'), ('A',)}
        assert DFS == path.dfs
        assert path.rfs == set([('A',)])

        cpath = [[Counter({'A': 1}), Counter({'i29': 1})], 
                 [Counter({'A': 1}), Counter({'i29': 1})],
                 [Counter({'i29': 1}), Counter({'i46': 1, 'i47': 1})], 
                 [Counter({'i47': 1}), Counter({'i70': 1})], 
                 [Counter({'i29': 1}), Counter({'i46': 1, 'i47': 1})], 
                 [Counter({'i47': 1}), Counter({'i70': 1})], 
                 [Counter({'i46': 1}), Counter({'i126': 1, 'i127': 1})], 
                 [Counter({'i70': 1}), Counter({'i47': 1})], 
                 [Counter({'i70': 1}), Counter({'i47': 1})], 
                 [Counter({'i47': 1}), Counter({'i70': 1})], 
                 [Counter({'i46': 1, 'i47': 1}), Counter({'i29': 1})], 
                 [Counter({'i70': 1}), Counter({'A': 1, 'i89': 1})], 
                 [Counter({'i29': 1}), Counter({'i46': 1, 'i47': 1})], 
                 [Counter({'i47': 1}), Counter({'i70': 1})], 
                 [Counter({'A': 1, 'i89': 1}), Counter({'i70': 1})], 
                 [Counter({'i70': 1}), Counter({'i47': 1})], 
                 [Counter({'i70': 1}), Counter({'i47': 1})], 
                 [Counter({'i47': 1}), Counter({'i70': 1})], 
                 [Counter({'i46': 1, 'i47': 1}), Counter({'i29': 1})], 
                 [Counter({'i70': 1}), Counter({'A': 1, 'i89': 1})], 
                 [Counter({'i29': 1}), Counter({'i46': 1, 'i47': 1})], 
                 [Counter({'i47': 1}), Counter({'i70': 1})], 
                 [Counter({'A': 1, 'i89': 1}), Counter({'i70': 1})], 
                 [Counter({'i70': 1}), Counter({'i47': 1})], 
                 [Counter({'i70': 1}), Counter({'i47': 1})]]
        path = Path(cpath, fs)
        DFS = {('i126', 'i127', 'i47'), ('i46', 'i47')}
        assert DFS == path.dfs
        assert path.rfs == frozenset({(), ('A',)})

        cpath = [[Counter({'B': 1}), Counter({'i11': 1})], 
                 [Counter({'i11': 1}), Counter({'B': 1})], 
                 [Counter({'B': 1}), Counter({'i11': 1})]] 
        path = Path(cpath, fs)
        RFS = set([(), ('B',)])
        assert RFS == path.rfs
        DFS = set([('B',), ('i11',)])
        assert DFS == path.dfs

        lpath = [[['A'], ['i29']], 
                 [['i29'], ['i46', 'i47']], 
                 [['i47'], ['i70']], 
                 [['i70'], ['i89', 'A']], 
                 [['A'], ['i29']], 
                 [['i29'], ['i46', 'i47']], 
                 [['i47'], ['i70']], 
                 [['i89', 'A'], ['i70']], 
                 [['i70'], ['i47']], 
                 [['i70'], ['i47']], 
                 [['i46', 'i47'], ['i29']]]
        path = Path(lpath, fs)
        DFS = {('i46', 'i47'), ('i29',)}
        assert DFS == path.dfs
        assert path.rfs == frozenset({()})

        lpath =   [[['A'], ['i9']], 
                   [['A'], ['i9']], 
                   [['i9'], ['i6', 'i7']], 
                   [['i7'], ['i0']], 
                   [['i9'], ['i6', 'i7']], 
                   [['i7'], ['i0']], 
                   [['i6'], ['X', 'Y']], 
                   [['i0'], ['i7']]]
        path = Path(lpath, fs)
        DFS = {('i0', 'i6'), ('X', 'Y', 'i0'), ('X', 'Y', 'i7'), ('i6', 'i7')}
        assert DFS == path.dfs
        assert path.rfs == frozenset({('X', 'Y')})

        lpath = [[['A'], ['i']],
                 [['B', 'i'], ['j']],
                 [['j'], ['X', 'k']]]
        path = Path(lpath, fs)
        RFS = set([('X',)])
        assert RFS == path.rfs
        assert path.dfs == set()

        p0 = 'A -> i; i -> j + B; B + j -> k'
        path = Path(my_parse_crn(p0)[0], fs)
        RFS = set([('B',)])
        assert RFS == path.rfs
        assert path.dfs == set()

        lpath = [[['A'],['i']],
                 [['i'],['A', 'j']],
                 [['A', 'j'], ['B']]]
        path = Path(lpath, fs)
        RFS = set([('A', 'B'), ('B',)])
        assert RFS == path.rfs
        assert path.dfs == set()

@unittest.skipIf(SKIP, "skipping tests")
class TestHelperFunctions(unittest.TestCase):
    # from nuskell.verifier.basis_finder import crn_properties
    # from nuskell.verifier.basis_finder import get_formal_basis
    def test_my_parse_crn(self): 
        crn = """
            A <=> i
            i + B1 <=> j1
            i + B2 <=> j2
            j1 -> C
            j2 -> C
            """
        icrn, isp = my_parse_crn(crn)
        ecrn = [[['A'], ['i']], 
                [['i'], ['A']], 
                [['i', 'B1'], ['j1']], 
                [['j1'], ['i', 'B1']], 
                [['i', 'B2'], ['j2']],
                [['j2'], ['i', 'B2']],
                [['j1'], ['C']],
                [['j2'], ['C']]]
        esp = set(['B2', 'j1', 'C', 'A', 'B1', 'j2', 'i'])
        assert icrn == ecrn
        assert isp == esp

    def test_clean_crn(self):
        from nuskell.verifier.basis_finder import clean_crn
        crn = [[['A', 'C'], ['C', 'A']], 
                [['A', 'C'], ['C', 'C']], 
                [['A'], ['A']], 
                [['A', 'B'], ['A', 'B']], 
                [['C', 'A'], ['C', 'C']], 
                [['A', 'B'], ['C', 'D']], 
                [['D'], ['D']]]

        new1 = [[['A', 'C'], ['A', 'C']], 
                [['A', 'C'], ['C', 'C']], 
                [['A'], ['A']], 
                [['A', 'B'], ['A', 'B']], 
                [['A', 'B'], ['C', 'D']], 
                [['D'], ['D']]]
        assert new1 == clean_crn(crn, duplicates = True, trivial = False)

        new1 = [[['A', 'C'], ['C', 'C']], 
                [['A', 'C'], ['C', 'C']], 
                [['A', 'B'], ['C', 'D']]]
        assert new1 == clean_crn(crn, duplicates = False, trivial = True)

        new1 = [[['A', 'C'], ['C', 'C']], 
                [['A', 'B'], ['C', 'D']]]
        assert new1 == clean_crn(crn)

    def test_formal_intermediate(self):
        from nuskell.verifier.basis_finder import formal, intermediate, is_formal_state
        fs = set(['A', 'B', 'X'])
        s1 = list(Counter({'A': 2, 'i': 1, 'X': 1, 'k': 34}).elements())
        s2 = list(Counter({'A': 2, 'X': 1}).elements())
        s3 = list(Counter({'i': 1, 'k': 34}).elements())
        assert sorted(s2) == formal(s1, fs)
        assert sorted(s3) == intermediate(s1, fs)
        assert is_formal_state(s2, fs)
        assert not is_formal_state(s1, fs)

    def test_interpret(self):
        from nuskell.verifier.basis_finder import interpret
        fs = set(['A', 'B', 'X'])
        inter = {'i': ['A'], 
                 'j': [], 
                 'k': ['B', 'X']}
        s1 = list(Counter({'A': 2, 'i': 1, 'X': 1, 'j': 15, 'k': 4}).elements())
        s2 = list(Counter({'A': 2, 'X': 1}).elements())
        s3 = list(Counter({'i': 1, 'k': 4}).elements())

        i1 = list(Counter({'A': 3, 'X': 5, 'B': 4}).elements())
        i2 = list(Counter({'A': 2, 'X': 1}).elements())
        i3 = list(Counter({'A': 1, 'B': 4, 'X': 4}).elements())

        assert sorted(i1) == sorted(interpret(s1, inter))
        assert sorted(i2) == sorted(interpret(s2, inter))
        assert sorted(i3) == sorted(interpret(s3, inter))

    def test_is_subset(self):
        from nuskell.verifier.basis_finder import is_subset
        fs = Counter(['A', 'B', 'X'])
        s1 = Counter({'A': 2, 'i': 1, 'X': 1, 'B': 15, 'k': 4})
        s2 = Counter({'A': 2, 'X': 1})
        s3 = Counter({'i': 1, 'k': 4})

        assert is_subset(fs, s1)
        assert not is_subset(fs, s2)
        assert not is_subset(fs, s3)

    def test_is_tidy(self):
        from nuskell.verifier.basis_finder import tidy
        fs = set(['A', 'X'])
        crn = [[['X'], ['i', 'k', 'l']],
               [['l'], ['p', 'q']],
               [['i', 'k', 'l'], ['A']],
               [['p', 'q'], ['k']],
               [['k'], ['l']]]
        T = ['i', 'k', 'l']
        assert tidy([T], crn, fs)

        crn = [[['X'], ['i', 'k', 'l']],
               [['l'], ['p', 'q']],
               [['i', 'k', 'l'], ['A']],
               [['k'], ['l']],
               [['p', 'q'], ['k']]]
        T = ['i', 'k', 'l']
        assert tidy([T], crn, fs)

        crn = [[['X'], ['i', 'k', 'l']],
               [['l'], []],
               [['l'], ['l', 'l', 'l']],
               [['i', 'k', 'l', 'l', 'l', 'l'], ['A']],
               [['k'], ['j']]]
        T = ['i', 'k', 'l']
        assert tidy([T], crn, fs)

        crn = """
        i -> i + q
        i -> i + i
        q -> q + k
        10k + i + q -> X
        10i -> i
        2q -> q
        """
        crn, _ = my_parse_crn(crn)
        T = ['i']
        assert tidy([T], crn, fs)

        crn = """
        m -> 3i
        i -> 3k 
        k -> 2k + l
        3l -> X
        2i ->
        2k ->
        """
        crn, _ = my_parse_crn(crn)
        T = ['m']
        assert tidy([T], crn, fs)

        crn = """
        a -> a + 3i
        7i -> b + 3l
        7l -> a
        2a + 3b + 2l -> X
        """
        crn, _ = my_parse_crn(crn)
        T = ['a']
        assert tidy([T], crn, fs)

        # This tests a modification to the original codebase:
        crn = """
        -> i
        2i -> A
        """
        crn, _ = my_parse_crn(crn)
        T = ['i']
        assert tidy([T], crn, fs)

    def test_get_crn_modules(self):
        from nuskell.verifier.basis_finder import get_crn_modules
        fs = set(['A', 'B', 'C'])
        crn = """
        A <=> i503
        A + i383 -> i420
        A + i407 -> i420
        A + i472 <=> i420
        A <=> i683
        A + i383 -> i420
        A + i407 -> i420
        B <=> i165
        B + i43 <=> i17
        B + i503 -> i17
        B + i683 -> i17
        B <=> i157
        B + i503 -> i17
        B + i683 -> i17
        C <=> i383
        C + i157 -> i178
        C + i165 -> i178
        C + i232 <=> i178
        C <=> i407
        C + i157 -> i178
        C + i165 -> i178
        i232 -> C
        i43 -> B
        i472 -> A
        """
        crn, _ = my_parse_crn(crn)
        species = set().union(*[set().union(*rxn) for rxn in crn])
        intermediates = species - fs

        modules = get_crn_modules(crn, intermediates)
        assert len(modules) == 3

        m1 = """
        A <=> i503
        B + i503 -> i17
        B + i503 -> i17
        B + i43 <=> i17
        B + i683 -> i17
        B + i683 -> i17
        i43 -> B
        A <=> i683
        """
        m2 = """
        A + i383 -> i420
        A + i383 -> i420
        C <=> i383
        A + i407 -> i420
        A + i407 -> i420
        A + i472 <=> i420
        C <=> i407
        i472 -> A
        """
        m3 = """
        B <=> i165
        B <=> i157
        C + i157 -> i178
        C + i165 -> i178
        C + i232 <=> i178
        C + i157 -> i178
        C + i165 -> i178
        i232 -> C
        """
        m1, _ = my_parse_crn(m1)
        m2, _ = my_parse_crn(m2)
        m3, _ = my_parse_crn(m3)

        m1 = sorted([sorted(r), sorted(p)] for [r, p] in m1)
        m2 = sorted([sorted(r), sorted(p)] for [r, p] in m2)
        m3 = sorted([sorted(r), sorted(p)] for [r, p] in m3)

        sm = []
        for m in modules:
            sm.append(sorted([sorted(r), sorted(p)] for [r, p] in m))

        assert m1 in sm
        assert m2 in sm
        assert m3 in sm

@unittest.skipIf(SKIP, "skipping tests")
class TestTidyAccumulation(unittest.TestCase):
    def test_tidy_binary_counter(self):
        from nuskell.verifier.basis_finder import tidy
        # An example with width > 2^K
        # With every bit increase (C_K)
        # the species A accumulates by 2^K
        fs = set(['N', 'M'])
        crn = """
        N -> i
        i -> D + x00 + x10 + x20 + x30
        F + x00 + x10 + x20 + x30 -> M

        D -> C0 + A
        C4 -> F
        F + A -> F
        C0 + x00 -> D  + x01
        C0 + x01 -> C1 + x00
        C1 + x10 -> D  + x11
        C1 + x11 -> C2 + x10
        C2 + x20 -> D  + x21
        C2 + x21 -> C3 + x20
        C3 + x30 -> D  + x31
        C3 + x31 -> C4 + x30
        """
        crn, _ = my_parse_crn(crn)
        S = ['i', 'i', 'i']
        T = ['M']
        assert tidy([S], crn, fs)

    def test_nontidy_binc(self):
        from nuskell.verifier.basis_finder import tidy
        fs = set(['N', 'M'])
        crn = """
        N -> i
        i -> D + x00 + x10 + x20 + x30
        F + x00 + x10 + x20 + x30 -> M

        D -> C0 + A
        C4 -> F
        F + A -> F
        C0 + x00 -> D  + x01
        C0 + x01 -> C1 + x00
        C1 + x10 -> D  + x11
        C1 + x11 -> C2 + x10
        C2 + x20 -> D  + x21
        C2 + x21 -> C3 + x20
        C3 + x30 -> D  + x31
        C3 + x31 -> C4 + x30
        """
        crn, _ = my_parse_crn(crn)
        #p = "N -> i; N -> i; N -> i; N -> i; N -> i; i -> D + x00 + x10 + x20 + x30; i -> D + x00 + x10 + x20 + x30; D -> A + C0; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x21 -> C3 + x20; C2 + x21 -> C3 + x20; C3 + x30 -> D + x31; C3 + x31 -> C4 + x30; C4 -> F; A + F -> F; D -> A + C0; A + F -> F; A + F -> F; A + F -> F; A + F -> F; i -> D + x00 + x10 + x20 + x30; A + F -> F; D -> A + C0; A + F -> F; A + F -> F; A + F -> F; A + F -> F; i -> D + x00 + x10 + x20 + x30; A + F -> F; D -> A + C0; A + F -> F; A + F -> F; A + F -> F; A + F -> F; i -> D + x00 + x10 + x20 + x30; A + F -> F; D -> A + C0; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; C0 + x00 -> D + x01; D -> A + C0; A + F -> F; C0 + x00 -> D + x01; D -> A + C0; A + F -> F; F + x00 + x10 + x20 + x30 -> M; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C1 + x11 -> C2 + x10; C1 + x11 -> C2 + x10; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; C2 + x20 -> D + x21; C2 + x21 -> C3 + x20; C2 + x21 -> C3 + x20; C3 + x30 -> D + x31; C3 + x31 -> C4 + x30; C4 -> F; A + F -> F; D -> A + C0; A + F -> F; D -> A + C0; A + F -> F; D -> A + C0; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; F + x00 + x10 + x20 + x30 -> M; C0 + x00 -> D + x01; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C1 + x10 -> D + x11; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; C2 + x21 -> C3 + x20; C2 + x21 -> C3 + x20; C3 + x30 -> D + x31; C3 + x31 -> C4 + x30; C4 -> F; A + F -> F; D -> A + C0; A + F -> F; D -> A + C0; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; F + x00 + x10 + x20 + x30 -> M; C0 + x00 -> D + x01; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x21 -> C3 + x20; C2 + x21 -> C3 + x20; C3 + x30 -> D + x31; C3 + x31 -> C4 + x30; C4 -> F; A + F -> F; D -> A + C0; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; A + F -> F; F + x00 + x10 + x20 + x30 -> M; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x21 -> C3 + x20; C3 + x30 -> D + x31; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x20 -> D + x21; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x10 -> D + x11; D -> A + C0; C0 + x00 -> D + x01; D -> A + C0; C0 + x01 -> C1 + x00; C1 + x11 -> C2 + x10; C2 + x21 -> C3 + x20; C3 + x31 -> C4 + x30; C4 -> F; F + x00 + x10 + x20 + x30 -> M"
        #path = Path(my_parse_crn(p)[0], fs)
        #T = path.Sn
        T = ['A', 'A']
        assert not tidy([[x for x in T if x not in fs]], crn, fs)

    def test_tidy_accumulation_01(self):
        from nuskell.verifier.basis_finder import tidy
        fs = set(['A', 'B', 'C', 'X'])
        crn = """
        X -> i + k + l
        l -> 
        l -> 3l
        i -> i + i
        B + i + k + 4l -> A
        """
        crn, _ = my_parse_crn(crn)
        T = ['i', 'k', 'l']
        assert tidy([T], crn, fs, bound = 50) is not True

    def test_tidy_accumulation_02(self):
        from nuskell.verifier.basis_finder import tidy
        fs = set(['A', 'B', 'C', 'X'])
        crn = """
        i -> 8i + q
        q -> 8q + k
        i + q + k -> C
        2i -> i
        2q -> q
        """
        crn, _ = my_parse_crn(crn)
        T = ['i']
        assert tidy([T], crn, fs, bound = 9) is not True
        assert tidy([T], crn, fs, bound = 10)

    def test_tidy_accumulation_03(self):
        # This test is somewhat misleading. Due to optimization,
        # -> k; k -> A and -> i would be in different modules ...
        from nuskell.verifier.basis_finder import tidy
        fs = set(['A'])
        crn = """
        -> k
        -> i
        k -> A
        """
        crn, _ = my_parse_crn(crn)
        T = ['i']
        assert tidy([T], crn, fs, bound = 10) is not True

    def test_a_large_example(self):
        # Comes from cat crns/basic/basic_05.crn | nuskell --ts lakin2012_3D.ts --verify pathway --max-complex-count 10000 --max-reaction-count 10000 --verify-timeout 0 -vv
        # formal CRN: A <=> 
        from nuskell.verifier.basis_finder import tidy
        fs = {'i219', 'i3173', 'i4464', 'i335', 'i3586', 
              'A_1_', 'i5223', 'i3400', 'i508', 'i2743', 
              'i108', 'i2525', 'i4989', 'i5467', 'i4598', 
              'i5115', 'i3786', 'i3717', 'i5325', 'i4696', 'i3296'}
        crn = """
         -> i39
        A_1_ -> i757
        A_1_ + i1945 -> i1681
        A_1_ + i1959 -> i1649
        i104 -> i165
        i104 + i1092 -> i165 + i801
        i104 + i1100 -> i355 + i801
        i104 + i1124 -> i289 + i801
        i104 + i1132 -> i173 + i801
        i104 + i1200 -> i1124 + i801
        i104 + i1226 -> i1124
        i104 + i1234 -> i1124 + i77
        i104 + i1649 -> i1603 + i77
        i104 + i165 -> i289 + i77
        i104 + i1681 -> i1603
        i104 + i2942 -> i3507 + i77
        i104 + i2950 -> i3507
        i104 + i2969 -> i3507 + i801
        i104 + i355 -> i289
        i104 + i39 -> i173 + i77
        i104 + i412 -> i355
        i104 + i4170 -> i4877 + i77
        i104 + i4181 -> i4877
        i104 + i4197 -> i4877 + i801
        i104 + i594 -> i165 + i77
        i104 + i610 -> i355 + i77
        i104 + i757 -> i1603 + i801
        i104 + i97 -> i173
        i1092 -> i39 + i801
        i1092 + i2401 -> i4170 + i801
        i1092 + i2406 -> i2942 + i801
        i1100 -> i801 + i97
        i1100 + i2401 -> i4181 + i801
        i1100 + i2406 -> i2950 + i801
        i1124 -> i104 + i1226
        i1124 -> i297 + i801
        i1124 + i2401 -> i4250 + i801
        i1124 + i2406 -> i2985 + i801
        i1124 + i77 -> i104 + i1234
        i1124 + i801 -> i104 + i1200
        i1132 -> i1226
        i1132 -> i76 + i801
        i1132 + i2401 -> i4219 + i801
        i1132 + i2406 -> i2993 + i801
        i1132 + i77 -> i1234
        i1132 + i801 -> i1200
        i1200 -> i1132 + i801
        i1200 + i2401 -> i4197 + i801
        i1200 + i2406 -> i2969 + i801
        i1226 -> i1132
        i1226 + i2401 -> i4197
        i1226 + i2406 -> i2969
        i1234 -> i1132 + i77
        i1234 + i2401 -> i4197 + i77
        i1234 + i2406 -> i2969 + i77
        i1603 -> i104 + i1681
        i1603 -> i2401 + i2406 + i2525
        i1603 + i77 -> i104 + i1649
        i1603 + i801 -> i104 + i757
        i1617 -> i1681
        i1617 -> i2401 + i2406 + i2743
        i1617 + i77 -> i1649
        i1617 + i801 -> i757
        i1649 -> A_1_ + i1959
        i1649 -> i1617 + i77
        i1649 + i2401 -> i4299 + i77
        i1649 + i2406 -> i3009 + i77
        i165 -> i104
        i165 -> i297 + i77
        i165 + i2401 -> i4250 + i77
        i165 + i2406 -> i2985 + i77
        i165 + i77 -> i104 + i594
        i165 + i801 -> i104 + i1092
        i1681 -> A_1_ + i1945
        i1681 -> i1617
        i1681 + i2401 -> i4299
        i1681 + i2406 -> i3009
        i173 -> A_1_ + i104 + i219
        i173 -> i104 + i97
        i173 + i77 -> i104 + i39
        i173 + i801 -> i104 + i1132
        i2401 -> i4170
        i2401 + i2942 -> i4189 + i77
        i2401 + i2950 -> i4189
        i2401 + i2969 -> i4189 + i801
        i2401 + i355 -> i4250
        i2401 + i39 -> i4219 + i77
        i2401 + i412 -> i4181
        i2401 + i4170 -> i4869 + i77
        i2401 + i4181 -> i4869
        i2401 + i4197 -> i4869 + i801
        i2401 + i594 -> i4170 + i77
        i2401 + i610 -> i4181 + i77
        i2401 + i757 -> i4299 + i801
        i2401 + i97 -> i4219
        i2406 -> i2942
        i2406 + i2942 -> i3478 + i77
        i2406 + i2950 -> i3478
        i2406 + i2969 -> i3478 + i801
        i2406 + i355 -> i2985
        i2406 + i39 -> i2993 + i77
        i2406 + i412 -> i2950
        i2406 + i4170 -> i4903 + i77
        i2406 + i4181 -> i4903
        i2406 + i4197 -> i4903 + i801
        i2406 + i594 -> i2942 + i77
        i2406 + i610 -> i2950 + i77
        i2406 + i757 -> i3009 + i801
        i2406 + i97 -> i2993
        i289 -> A_1_ + i104 + i508
        i289 -> i104 + i355
        i289 + i77 -> i104 + i165
        i289 + i801 -> i104 + i1124
        i2942 -> i2406
        i2942 -> i3455 + i77
        i2942 + i77 -> i2406 + i594
        i2942 + i801 -> i1092 + i2406
        i2950 -> i2406 + i412
        i2950 -> i3455
        i2950 + i77 -> i2406 + i610
        i2950 + i801 -> i1100 + i2406
        i2969 -> i1226 + i2406
        i2969 -> i3455 + i801
        i2969 + i77 -> i1234 + i2406
        i2969 + i801 -> i1200 + i2406
        i297 -> A_1_ + i104 + i335
        i297 -> i355
        i297 + i77 -> i165
        i297 + i801 -> i1124
        i2985 -> A_1_ + i104 + i3400
        i2985 -> i2406 + i355
        i2985 + i77 -> i165 + i2406
        i2985 + i801 -> i1124 + i2406
        i2993 -> A_1_ + i104 + i3296
        i2993 -> i2406 + i97
        i2993 + i77 -> i2406 + i39
        i2993 + i801 -> i1132 + i2406
        i3009 -> i1681 + i2406
        i3009 -> i2401 + i2406 + i3173
        i3009 + i77 -> i1649 + i2406
        i3009 + i801 -> i2406 + i757
        i3455 -> A_1_ + i104 + i3786
        i3455 -> i2950
        i3455 + i77 -> i2942
        i3455 + i801 -> i2969
        i3478 -> A_1_ + i104 + i3717
        i3478 -> i2406 + i2950
        i3478 + i77 -> i2406 + i2942
        i3478 + i801 -> i2406 + i2969
        i3507 -> A_1_ + i104 + i3586
        i3507 -> i104 + i2950
        i3507 + i77 -> i104 + i2942
        i3507 + i801 -> i104 + i2969
        i355 -> i104 + i412
        i355 -> i297
        i355 + i77 -> i104 + i610
        i355 + i801 -> i104 + i1100
        i39 ->
        i39 -> i76 + i77
        i39 + i77 -> i594
        i39 + i801 -> i1092
        i412 -> i97
        i4170 -> i2401
        i4170 -> i4914 + i77
        i4170 + i77 -> i2401 + i594
        i4170 + i801 -> i1092 + i2401
        i4181 -> i2401 + i412
        i4181 -> i4914
        i4181 + i77 -> i2401 + i610
        i4181 + i801 -> i1100 + i2401
        i4189 -> A_1_ + i104 + i5467
        i4189 -> i2401 + i2950
        i4189 + i77 -> i2401 + i2942
        i4189 + i801 -> i2401 + i2969
        i4197 -> i1226 + i2401
        i4197 -> i4914 + i801
        i4197 + i77 -> i1234 + i2401
        i4197 + i801 -> i1200 + i2401
        i4219 -> A_1_ + i104 + i4696
        i4219 -> i2401 + i97
        i4219 + i77 -> i2401 + i39
        i4219 + i801 -> i1132 + i2401
        i4250 -> A_1_ + i104 + i4598
        i4250 -> i2401 + i355
        i4250 + i77 -> i165 + i2401
        i4250 + i801 -> i1124 + i2401
        i4299 -> i1681 + i2401
        i4299 -> i2401 + i2406 + i4464
        i4299 + i77 -> i1649 + i2401
        i4299 + i801 -> i2401 + i757
        i4869 -> A_1_ + i104 + i5325
        i4869 -> i2401 + i4181
        i4869 + i77 -> i2401 + i4170
        i4869 + i801 -> i2401 + i4197
        i4877 -> A_1_ + i104 + i5223
        i4877 -> i104 + i4181
        i4877 + i77 -> i104 + i4170
        i4877 + i801 -> i104 + i4197
        i4903 -> A_1_ + i104 + i5115
        i4903 -> i2406 + i4181
        i4903 + i77 -> i2406 + i4170
        i4903 + i801 -> i2406 + i4197
        i4914 -> A_1_ + i104 + i4989
        i4914 -> i4181
        i4914 + i77 -> i4170
        i4914 + i801 -> i4197
        i594 -> i39 + i77
        i610 -> i77 + i97
        i757 -> A_1_
        i757 -> i1617 + i801
        i76 -> A_1_ + i104 + i108
        i76 -> i97
        i76 + i77 -> i39
        i76 + i801 -> i1132
        i77 + i97 -> i610
        i801 + i97 -> i1100
        i97 -> i412
        i97 -> i76
        """
        crn, _ = my_parse_crn(crn)
        T = ['A_1_', 'i1945', 'i801']
        assert tidy([[t for t in T if t not in fs]], crn, fs, bound = 5) is not True


@unittest.skipIf(SKIP, "skipping tests")
class TestBasisFinder(unittest.TestCase):
    def test_no_basis_found(self):
        pf = 'A -> B; B -> C'
        pi = 'A <=> j; j -> i + C; C + i -> B; B -> C; i->j'
        p0 = 'A -> j; j -> i + C; C + i -> B'
        fcrn, fs = my_parse_crn(pf)
        icrn, _ = my_parse_crn(pi)
        with self.assertRaises(NoFormalBasisError) as e:
            basis_raw, basis_int = find_basis(icrn, fs)

        cf = 'A + B <=> C + D'
        ci = 'A <=> i; i + B <=> j; j <=> k + C; k <=> D'
        fcrn, fs = my_parse_crn(cf)
        icrn, _ = my_parse_crn(ci)
        with self.assertRaises(NoFormalBasisError) as e:
            basis_raw, basis_int = find_basis(icrn, fs)

    def test_STW2019_F3(self):
        crn1 = "A + B -> D + E; D -> F"
        crn2 = """
            A <=> i
            B <=> j
            i + j  -> k
            k -> D + l
            l -> E
            D <=> m
            m -> F
            """
        fcrn, fs = my_parse_crn(crn1)
        icrn, _ = my_parse_crn(crn2)

        basis_raw = [[['A', 'B'], ['D', 'E']], [['D'], ['F']]]
        basis_raw = sorted(sorted(rxn) for rxn in basis_raw)
        basis_int = []

        basis_raw1, basis_int1 = find_basis(icrn, fs, modular = True)
        basis_raw1 = sorted(sorted(rxn) for rxn in basis_raw1)

        basis_raw2, basis_int2 = find_basis(icrn, fs, modular = False)
        basis_raw2 = sorted(sorted(rxn) for rxn in basis_raw1)

        assert basis_raw1 == basis_raw2 == basis_raw
        assert basis_int1 == basis_int2 == basis_int

    def test_STW2019_F4(self):
        crn1 = "A -> B"
        crn2 = """
            A1 -> i
            i -> B1 + W
            A2 -> j
            j -> B2
            W + j -> B1
            """
        fcrn, fs = my_parse_crn(crn1)
        icrn, _ = my_parse_crn(crn2)
        
        inter = {'A1': ['A'], 
                 'A2': ['A'], 
                 'B1': ['B'],
                 'B2': ['B'],
                 'W': []} 

        braw1 = braw2 = [[['A1'], ['B1', 'W']], [['A2'], ['B2']], [['A2', 'W'], ['B1']]]
        bint1 = []
        bint2 = [[['A'], ['B']]]

        fs2 = set(inter.keys()) 
        basis_raw1, basis_int1 = find_basis(icrn, fs2, interpretation = None)

        assert sorted(braw1) == sorted(basis_raw1)
        assert sorted(bint1) == sorted(basis_int1)

        basis_raw2, basis_int2 = find_basis(icrn, fs2, interpretation = inter)
        assert sorted(braw2) == sorted(basis_raw2)
        assert sorted(bint2) == sorted(basis_int2)

    def test_STW2019_intro(self):
        crn1 = "A+B -> C+D; C+A -> C+C"
        crn2 = "A<=>i; i+B<=>j; i+j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
        crn3 = "A<=>i; i+B<=>j; j<=>C+k; k->D; C+A<=>m+n; m+n->C+C"
        crn4 = "A->i; i+B<=>j; j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
        crn5 = "A<=>i; i+B<=>j; j->C+k; k<=>D; C+A<=>m+n; m+n->C+C"
        crn6 = "A+g1<=>i+g2; i+B<=>j+g3; g4+j->C+k+w1; g5+k<=>D+w2; C+A<=>m+n; g6+m+n->C+C+w3"

        crn1, fs = my_parse_crn(crn1)
        crn2, _ = my_parse_crn(crn2)
        crn3, _ = my_parse_crn(crn3)
        crn4, _ = my_parse_crn(crn4)
        crn5, _ = my_parse_crn(crn5)
        crn6, _ = my_parse_crn(crn6)

        basis1, _ = find_basis(crn1, fs)
        basis1 = sorted(sorted(rxn) for rxn in basis1)

        with self.assertRaises(NoFormalBasisError) as e:
            find_basis(crn2, fs)

        with self.assertRaises(NoFormalBasisError) as e:
            find_basis(crn3, fs)

        with self.assertRaises(NoFormalBasisError) as e:
            find_basis(crn4, fs)

        basis5, _ = find_basis(crn5, fs)
        basis5 = sorted(sorted(rxn) for rxn in basis5)
        assert basis5 == basis1

        basis6, _ = find_basis(crn6, fs)
        basis6 = sorted(sorted(rxn) for rxn in basis6)
        assert basis6 != basis1

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

        inter1 = {'A': ['A'], 
                 'B1': ['B'],
                 'B2': ['B'],
                 'C': ['C']}
        fs1 = set(inter1.keys()) 

        inter2 = {'A': ['A'], 
                 'B1': ['B'],
                 'B2': ['B'],
                 'i':  ['A'],
                 'j1': ['A', 'B'],
                 'j2': ['A', 'B'],
                 'C': ['C']}
        fs2 = set(inter2.keys()) 
        icrn = sorted([sorted(r), sorted(p)] for [r, p] in icrn)

        basis_raw, basis_int = find_basis(icrn, fs)
        assert basis_raw == []
        assert basis_int == []

        with self.assertRaises(NoFormalBasisError) as e:
            basis_raw, basis_int = find_basis(icrn, fs1, interpretation = None)

        basis_raw, basis_int = find_basis(icrn, fs2, interpretation = None)
        assert sorted(basis_raw) == icrn
        assert basis_int == []

        basis_raw, basis_int = find_basis(icrn, fs1, interpretation = inter1)
        br = [[['A', 'B1'], ['C']], 
               [['A', 'B1', 'B2'], ['B2', 'C']], 
               [['A', 'B1', 'B2'], ['B1', 'C']], 
               [['A', 'B2'], ['C']]]
        assert sorted(basis_raw) == sorted(br)
        assert sorted(basis_int) == fcrn

        basis_raw, basis_int = find_basis(icrn, fs2, interpretation = inter2)
        assert sorted(basis_raw) == icrn
        assert sorted(basis_int) == fcrn

    def test_crn_JDW2019_F5(self):
        crn1 = "A + B -> C + D"
        crn2 = "A + B <=> C + D"
        crn3 = """
        A <=> i_A_BCD
        B + i_A_BCD <=> i_AB_CD
        i_AB_CD <=> i_ABC_D + C
        i_ABC_D <=> i_ABCD_ + D
        i_ABCD_ + fi -> w_ABCD
        """
        icrn, _ = my_parse_crn(crn3)
        inter = {'A': ['A'], 
                 'B': ['B'], 
                 'C': ['C'], 
                 'D': ['D']}
        fs = set(inter.keys())
        with self.assertRaises(NoFormalBasisError) as e:
            basis_raw, basis_int = find_basis(icrn, fs)
        with self.assertRaises(NoFormalBasisError) as e:
            basis_raw, basis_int = find_basis(icrn, fs, interpretation = inter)

    @unittest.skip("skipping slow test")
    def test_lakin2016_2D_3I_accumulation(self):
        #TODO: does this terminate over night?
        fcrn = "B + C -> A"
        icrn = """
            B -> i47 + i48
            C + i47 -> i97 + i98
            i135 -> i560 + i561
            i135 + i154 -> i561 + i572
            i135 + i216 -> i560 + i582
            i136 -> i154 + i155
            i136 + i216 -> i155
            i136 + i582 -> i155 + i561
            i154 + i155 -> i136
            i154 + i216 ->
            i154 + i560 -> i572
            i154 + i582 -> i561
            i154 + i802 -> i136 + i572
            i155 -> i136 + i216
            i155 + i560 -> i216 + i802
            i216 + i802 -> i155 + i560
            i216 + i97 -> i132 + i135 + i155
            i47 + i48 -> B
            i560 + i582 -> i135 + i216
            i560 + i97 -> i132 + i135 + i802
            i582 -> i216 + i561
            i582 + i802 -> i135 + i155
            i802 -> i136 + i560
            i802 -> i155 + i572
            i97 -> i132 + i135 + i136
            i97 + i98 -> C + i47
            """
        fcrn, fs = my_parse_crn(fcrn)
        icrn, _ = my_parse_crn(icrn)
        with self.assertRaises(NoFormalBasisError) as e:
            basis_raw, basis_int = find_basis(icrn, fs)


class TestIntegratedHybrid(unittest.TestCase):
    pass

