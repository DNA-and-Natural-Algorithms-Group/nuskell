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
    def test_max_width_examples(self):
        #
        # TODO: find an example where (w + 1) * b = w_max
        #
        fs = set(['A', 'B', 'C', 'X', 'Y', 'Z'])

        w, b, wmax = 1, 2, 4
        p0 = 'A -> i; i -> j'
        path0 = Path(my_parse_crn(p0)[0], fs)
        # Any reaction with b == 2 that yields an undecomposable pathway can at
        # most increase the width to 2
        rxn1 = 'X + j -> l'
        #rxn1 = 'j -> l + k'
        path1 = Path(my_parse_crn(rxn1)[0], fs)
        path01 = path0 + path1
        assert path0.width == w
        assert path1.width == b
        assert path01.width == 2
        assert path01.width <= wmax

        w, b, wmax = 3, 4, 16
        p1 = 'A -> i; B -> j; B -> j; 2j + i -> k'
        path1 = Path(my_parse_crn(p1)[0], fs)
        rxn = '4k -> C' # yields w = 12
        rxn = Path(my_parse_crn(rxn)[0], fs)
        path = path1 + path1 + path1 + path1 + rxn
        # rxn = '3X + k-> C'
        # rxn = Path(my_parse_crn(rxn)[0], fs)
        # path = path1 + rxn # yields w = 6
        assert path1.width == w
        assert rxn.width == b
        assert path.width <= wmax # 12!

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
        assert path.rfs == frozenset({frozenset({('A', 1)})})
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
        RFS = frozenset([frozenset(Counter(x).items()) for x in RFS])
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
        RFS = frozenset([frozenset(Counter(x).items()) for x in RFS])
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
        assert pathX.rfs == frozenset({frozenset({('C', 1)}), frozenset({('B', 1), ('C', 1)})})

        pX = 'A -> i; B + i -> j; j -> X + k'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({frozenset({('X', 1)})})

        pX = 'A -> i; i -> B + j; B + j -> k'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({frozenset({('B', 1)})})

        pX = 'A -> m; m -> i; i -> A + j; j -> p; p -> j; A + j -> k; k -> l; l -> B'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({frozenset({('A', 1), ('B', 1)}), frozenset({('B', 1)})})

        pX = 'A -> i; i -> A + j; j -> p; p -> j; A + j -> B'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({frozenset({('A', 1), ('B', 1)}), frozenset({('B', 1)})})

        pX = 'A -> i; B + i -> j; j -> B + i; B + i -> j; j -> C + k; C + k -> j; j -> B + i; i -> A'
        pathX = Path(my_parse_crn(pX)[0], fs)
        assert pathX.rfs == frozenset({frozenset({('A', 1), ('C', 1), ('B', 1)})})

        # Temporarily non-regular (1/2)
        p0 = 'A -> i; B -> j; i + j -> k; k -> X + l; X -> m'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.rfs == frozenset({frozenset({('X', 1)})})
        # Temporarily non-regular (2/2)
        p0 = 'A -> i; B -> j; i + j -> k; k -> X + l; X -> m; m + l -> X + Y'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.rfs == frozenset({frozenset({('X', 1), ('Y', 1)})})

        # Two turning points, same result ...
        p0 = 'B -> i; A + i -> j + A; j -> k; A + k -> A + C'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.rfs == frozenset({frozenset({('C', 1), ('A', 1)})})

        # No turning point ...
        p0 = 'X -> i; B -> B + j; A -> A + p; i + j + p + B -> B + Y'
        path0 = Path(my_parse_crn(p0)[0], fs)
        assert path0.rfs == frozenset()

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
        assert pathX.signature == (frozenset({('A', 2)}), 
                                   frozenset({('Y', 2), ('i47', 1), ('i70', 1), ('X', 2)}), 
                                   7, 
                                   frozenset({('Y', 2), ('A', 2), ('X', 2)}), 
                                   frozenset({('X', 'Y', 'i47'), ('X', 'Y', 'i70')}), 
                                   frozenset())

    def test_deprecated_functions(self):
        fs = set('AVWXY')
        p0 = 'A -> i29; i29 -> i46 + i47; i47 -> i70; i46 -> X + Y; i70 -> A + i89; A -> i29; i29 -> i46 + i47; i47 -> i70; i70 -> A + i89; A -> i29'
        p1 = 'A -> i29; i29 -> i46 + i47; i47 -> i70; i46 -> X + Y; i70 -> A + i89; A -> i29; i29 -> i46 + i47; i47 -> i70; i70 -> A + i89; A -> i29; A + i89 -> i70'
        p2 = 'A -> i29; i29 -> i46 + i47; i47 -> i70; i46 -> X + Y; i70 -> A + i89; A -> i29; i29 -> i46 + i47; i47 -> i70; i70 -> A + i89; A -> i29; A + i89 -> i70; i70 -> i47'

        path0 = Path(my_parse_crn(p0)[0], fs)
        #print(path0)
        #print(old_signature(path0.lpath, fs))

        path1 = Path(my_parse_crn(p1)[0], fs)
        #print(path1)
        #print(old_signature(path1.lpath, fs))

        path2 = Path(my_parse_crn(p2)[0], fs)
        #print(path2)
        #print(old_signature(path2.lpath, fs))
        
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
        RFS = frozenset({frozenset({('A', 1), ('B', 1)})})
        assert path.rfs == RFS
 
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
        RFS = frozenset({frozenset({('A', 1)})})
        assert path.rfs == RFS

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
        RFS = frozenset({frozenset({('A', 1), ('X', 1)}), 
                         frozenset({('A', 1), ('B', 1), ('X', 1)})})
        assert RFS == path.rfs

        lpath = [[['A'],['i']],
                 [['i'],['A', 'j']],
                 [['A', 'j'], ['B']]]
        path = Path(lpath, fs)
        assert path.minimal_initial_state == Counter({'A': 1})
        assert path.final_state == Counter({'B': 1})
        assert path.width == 2
        assert path.formal_closure == Counter({'A': 1, 'B': 1})
        assert path.dfs == set()
        RFS = frozenset({frozenset({('B', 1)}), frozenset({('A', 1), ('B', 1)})})
        assert RFS == path.rfs

        cpath = [[Counter({'B': 1}), Counter({'i11': 1})], 
                 [Counter({'B': 1}), Counter({'i11': 1})]]
        path = Path(cpath, fs)
        assert path.minimal_initial_state == Counter({'B': 2}) 
        assert path.final_state == Counter({'i11': 2})
        assert path.width == 2
        assert path.formal_closure == Counter({'B': 2})
        assert path.dfs == frozenset({('i11',)})
        RFS = frozenset({frozenset()})
        assert RFS == path.rfs

        lpath = [[['A'], ['i']],
                 [['i'], ['j', 'j']],
                 [['B'], ['j']],
                 [['j'], ['B', 'C']],
                 [['j', 'j'], ['i']],
                 [['i'], ['C']]]
        path = Path(lpath, fs)
        RFS = set([('B', 'C', 'C')])
        RFS = frozenset([frozenset(Counter(x).items()) for x in RFS])
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
        RFS = frozenset({frozenset({('A', 1)})})
        assert RFS == path.rfs

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
        RFS = frozenset({frozenset(), frozenset({('A', 1)})})
        assert RFS == path.rfs

        cpath = [[Counter({'B': 1}), Counter({'i11': 1})], 
                 [Counter({'i11': 1}), Counter({'B': 1})], 
                 [Counter({'B': 1}), Counter({'i11': 1})]] 
        path = Path(cpath, fs)
        RFS = set([(), ('B',)])
        RFS = frozenset([frozenset(Counter(x).items()) for x in RFS])
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
        RFS = frozenset({frozenset()})
        assert RFS == path.rfs

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
        RFS = frozenset({frozenset({('X', 1), ('Y', 1)})})
        assert RFS == path.rfs

        lpath = [[['A'], ['i']],
                 [['B', 'i'], ['j']],
                 [['j'], ['X', 'k']]]
        path = Path(lpath, fs)
        RFS = set([('X',)])
        RFS = frozenset([frozenset(Counter(x).items()) for x in RFS])
        assert RFS == path.rfs
        assert path.dfs == set()

        p0 = 'A -> i; i -> j + B; B + j -> k'
        path = Path(my_parse_crn(p0)[0], fs)
        RFS = set([('B',)])
        RFS = frozenset([frozenset(Counter(x).items()) for x in RFS])
        assert RFS == path.rfs
        assert path.dfs == set()

        lpath = [[['A'],['i']],
                 [['i'],['A', 'j']],
                 [['A', 'j'], ['B']]]
        path = Path(lpath, fs)
        RFS = set([('A', 'B'), ('B',)])
        RFS = frozenset([frozenset(Counter(x).items()) for x in RFS])
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
        from nuskell.verifier.basis_finder import is_tidy
        fs = set(['A', 'X'])
        crn = [[['X'], ['i', 'k', 'l']],
               [['l'], ['p', 'q']],
               [['i', 'k', 'l'], ['A']],
               [['p', 'q'], ['k']],
               [['k'], ['l']]]
        T = ['i', 'k', 'l']
        assert is_tidy(T, crn, fs)

        crn = [[['X'], ['i', 'k', 'l']],
               [['l'], ['p', 'q']],
               [['i', 'k', 'l'], ['A']],
               [['k'], ['l']],
               [['p', 'q'], ['k']]]
        T = ['i', 'k', 'l']
        assert is_tidy(T, crn, fs)

        crn = [[['X'], ['i', 'k', 'l']],
               [['l'], []],
               [['l'], ['l', 'l', 'l']],
               [['i', 'k', 'l', 'l', 'l', 'l'], ['A']],
               [['k'], ['j']]]
        T = ['i', 'k', 'l']
        assert is_tidy(T, crn, fs)

        crn = [[['X'], ['i', 'k', 'l']],
               [['l'], []],
               [['l'], ['l', 'l', 'l']],
               [['B', 'i', 'k', 'l', 'l', 'l', 'l'], ['A']],
               [['k'], ['j']]]
        T = ['i', 'k', 'l']
        assert not is_tidy(T, crn, fs)

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

    def test_STW17_intro(self):
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

class TestIntegratedHybrid(unittest.TestCase):
    pass

