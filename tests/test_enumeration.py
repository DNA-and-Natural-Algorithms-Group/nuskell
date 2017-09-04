
import unittest
import argparse
from collections import Counter

import peppercornenumerator as pep
import peppercornenumerator.utils as peputils
import peppercornenumerator.reactions as reactions
from peppercornenumerator.enumerator import Enumerator
from peppercornenumerator.condense import condense_resting_states

import nuskell.enumeration as ne
import nuskell.parser as np
import nuskell.objects as objects

#from nuskell.parser import parse_crn_string, split_reversible_reactions
#import nuskell.verifier.crn_bisimulation_equivalence as bisimulation
#import nuskell.include.peppercorn.reactions as reactions


@unittest.skip('need to repair')
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
        self.args.MAX_COMPLEX_SIZE = 100
        self.args.MAX_COMPLEX_COUNT = 1000
        self.args.MAX_REACTION_COUNT = 1000

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
        pass

    def test_peppercorn_interface(self):
        """ Make sure peppercorn.utils did not change """

        t0 = peputils.Domain('t0', 5, sequence='H' * 5)
        t0_ = peputils.Domain('t0', 5, sequence='D' * 5, is_complement=True)
        d1 = peputils.Domain('d1', 15, sequence='H' * 15)
        d1_ = peputils.Domain('d1', 15, sequence='D' * 15, is_complement=True)
        domains = [t0, t0_, d1, d1_]

        s0 = peputils.Strand('s0', [t0, d1])
        s1 = peputils.Strand('s1', [d1])
        s2 = peputils.Strand('s2', [d1_, t0_])
        strands = [s0, s1, s2]

        c1s = peputils.parse_dot_paren('..')
        c1 = peputils.Complex('c1', [s0], c1s)
        c1.check_structure()

        c2s = peputils.parse_dot_paren('(+).')
        c2 = peputils.Complex('c2', [s1, s2], c2s)
        c2.check_structure()
        complexes = [c1, c2]

        enum = Enumerator(complexes)  # domains, strands, complexes)

        enum.enumerate()

        ###########################
        # Get full output CRN
        reactions = enum.reactions

        self.assertEqual(len(reactions), 3)

        r1_kernel = 't0 d1  +  d1( + ) t0*  ->  t0( d1 + d1( + ) )'
        r2_kernel = 't0( d1 + d1( + ) )  ->  t0( d1( + ) )  +  d1'
        r3_kernel = 't0( d1 + d1( + ) )  ->  t0 d1  +  d1( + ) t0*'

        p_cntr = Counter()
        r_cntr = Counter()
        for r in sorted(reactions):
            # print r.kernel_string()
            for cx in r.reactants:
                p_cntr += Counter(map(str, cx.strands))
            for cx in r.products:
                r_cntr += Counter(map(str, cx.strands))

            self.assertEqual(p_cntr, r_cntr)
            self.assertTrue(
                r.kernel_string() in [
                    r1_kernel,
                    r2_kernel,
                    r3_kernel])

        ###########################
        # Get condensed output CRN
        condensed = condense_resting_states(
            enum, compute_rates=True, k_fast=0.)
        reactions = condensed['reactions']

        self.assertEqual(len(reactions), 1)

        result_kernel = 't0 d1  +  d1( + ) t0*  ->  t0( d1( + ) )  +  d1'
        self.assertEqual(reactions[0].kernel_string(), result_kernel)

    def _TestTube_from_DOM(self, domstring):
        solution = objects.TestTube()

        [doms, cplxs] = np.parse_dom_string(domstring)
        domains = dict()
        for n, l in doms:
            domains[n] = objects.Domain(list('N' * int(l)), name=n)
            domains[n +
                    '*'] = domains[n].get_ComplementDomain(list('N' * int(l)))

        for n, se, ss in cplxs:
            for e, x in enumerate(se):
                x = ''.join(x)
                if x == '+':
                    continue
                se[e] = domains[x]

            cx = objects.Complex(se, ss, name=n)
            solution.add_complex(cx, (float('inf'), True))

        return solution

    @unittest.skip("needs to be fixed -> kernel string sometimes wrong")
    def test_under_the_hood1(self):
        # A test function that checks enumerated reaction networks by directly
        # calling the Enumerator object, not the output of the TestTubePeppercornIO
        # wrapper.
        domstring = """
    sequence tf : 5
    sequence bm : 15
    sequence tb : 5

    A :
    tf bm
    . .

    B :
    bm tb + tb* bm* tf*
    ( ( + ) ) .
    """

        solution = self._TestTube_from_DOM(domstring)
        peppercorn = ne.TestTubePeppercornIO(
            testtube=solution, pargs=self.args)
        self.assertIsInstance(pep.enumerator, Enumerator)

        pep.enumerate()

        ###########################
        # Get full output CRN
        reactions = pep.enumerator.reactions
        r1 = 'tf bm  +  bm( tb( + ) ) tf*  ->  tf( bm + bm( tb( + ) ) )'
        r2 = 'tf( bm( + tb* ) )  +  bm tb  ->  tf( bm( + bm tb( + ) ) )'
        r3 = 'tf( bm + bm( tb( + ) ) )  ->  tf( bm( + bm tb( + ) ) )'
        r4 = 'tf( bm( + bm tb( + ) ) )  ->  tf( bm + bm( tb( + ) ) )'
        r5 = 'tf( bm + bm( tb( + ) ) )  ->  tf bm  +  bm( tb( + ) ) tf*'
        r6 = 'tf( bm( + bm tb( + ) ) )  ->  tf( bm( + tb* ) )  +  bm tb'

        for r in sorted(reactions):
            print r.kernel_string()
            self.assertTrue(r.kernel_string() in [r1, r2, r3, r4, r5, r6])

        ###########################
        # Get condensed output CRN
        condensed = condense_resting_states(pep.enumerator)
        reactions = condensed['reactions']

        rc1 = 'tf bm  +  bm( tb( + ) ) tf*  ->  tf( bm( + tb* ) )  +  bm tb'
        rc2 = 'tf( bm( + tb* ) )  +  bm tb  ->  tf bm  +  bm( tb( + ) ) tf*'

        for r in sorted(reactions):
            self.assertTrue(r.kernel_string() in [rc1, rc2])

    @unittest.skip("needs to be fixed -> kernel string sometimes wrong")
    def test_under_the_hood2(self):
        domstring = """
    sequence tf : 5
    sequence bm : 15

    A :
    tf bm
    . .

    B :
    bm tf + tf* bm* tf*
    ( ( + ) ) .
    """

        solution = self._TestTube_from_DOM(domstring)
        peppercorn = ne.TestTubePeppercornIO(
            testtube=solution, pargs=self.args)
        self.assertIsInstance(pep.enumerator, Enumerator)

        pep.enumerate()

        ###########################
        # Get full output CRN
        reactions = pep.enumerator.reactions
        r1 = 'tf bm  +  tf( bm( + tf* ) )  ->  tf( bm( + tf( bm + ) ) )'
        r2 = 'tf bm  +  bm( tf( + ) ) tf*  ->  tf( bm + bm( tf( + ) ) )'
        r3 = 'tf( bm( + tf* ) )  +  bm tf  ->  tf( bm( + bm tf( + ) ) )'
        r4 = 'bm tf  +  bm( tf( + ) ) tf*  ->  bm tf( + bm( tf( + ) ) )'
        r5 = 'tf( bm + bm( tf( + ) ) )  ->  tf( bm( + bm tf( + ) ) )'
        r6 = 'tf( bm( + bm tf( + ) ) )  ->  tf( bm + bm( tf( + ) ) )'
        r7 = 'tf( bm( + tf( bm + ) ) )  ->  tf bm  +  tf( bm( + tf* ) )'
        r8 = 'tf( bm + bm( tf( + ) ) )  ->  tf bm  +  bm( tf( + ) ) tf*'
        r9 = 'tf( bm( + bm tf( + ) ) )  ->  tf( bm( + tf* ) )  +  bm tf'
        r10 = 'bm tf( + bm( tf( + ) ) )  ->  bm tf  +  bm( tf( + ) ) tf*'

        self.assertEqual(len(reactions), 10)
        for r in sorted(reactions):
            self.assertTrue(r.kernel_string() in [
                r1, r2, r3, r4, r5, r6, r7, r8, r9, r10])

        ###########################
        # Get condensed output CRN
        condensed = condense_resting_states(pep.enumerator)
        reactions = condensed['reactions']

        rc1 = 'tf bm  +  bm( tf( + ) ) tf*  ->  tf( bm( + tf* ) )  +  bm tf'
        rc2 = 'tf( bm( + tf* ) )  +  bm tf  ->  tf bm  +  bm( tf( + ) ) tf*'

        self.assertEqual(len(reactions), 2)
        for r in sorted(reactions):
            self.assertTrue(r.kernel_string() in [rc1, rc2])

    @unittest.skip("debugging stuff")
    def test_under_the_hood3(self):
        # Here we have a test example that illustrates the condese function of the
        # enumerator. The example uses complexes from the Soloveichik scheme,
        # reaction B+B->B. Note that some of these complexes are commented out!
        # It turns out that the enumerator does the correct thing, and the test is
        # kind of pointless. But it is nice to understand what is going on under the
        # hood of the condese function, and what problems you might find when
        # looking at the return statements.
        domstring = """
    sequence a : 6
    sequence x : 15
    sequence b : 6
    sequence y : 15

    B :
    y a x b
    . . . .

    #C_1_ :
    #x b a
    #. . .

    #C_2_ :
    #y a x b + a* y* b*
    #( ( . . + ) ) .

    C_3_ :
    x b a + x b y a + b* x* a* b* x* a*
    ( ( ( + ( ( . . + ) ) ) ) ) .
    """

        solution = self._TestTube_from_DOM(domstring)
        self.args.MAX_COMPLEX_SIZE = 4
        self.args.MAX_COMPLEX_COUNT = 9
        self.args.MAX_REACTION_COUNT = 10
        self.args.reject_remote = False
        self.args.ignore_branch_4way = True

        peppercorn = ne.TestTubePeppercornIO(
            testtube=solution, pargs=self.args)
        self.assertIsInstance(pep.enumerator, Enumerator)

        pep.enumerate()
        ###########################
        # Get full output CRN
        reactions = pep.enumerator.reactions

        count = 0
        ks_t_cp = dict()
        for r in sorted(reactions):
            print map(str, r.reactants), '->', map(str, r.products)
            print r.kernel_string()

        print '--'
        ###########################
        # Get condensed output CRN
        condensed = condense_resting_states(pep.enumerator)
        reactions = condensed['reactions']

        for r in sorted(reactions):
            print map(str, map(lambda x: x.complexes, r.reactants)), '->',
            print map(str, map(lambda x: x.complexes, r.products))
            print r.kernel_string()
            for rs in r.reactants:
                if len(rs.complexes) == 1:
                    print rs.complexes[0],
                else:
                    print 'rs' + str(rs),
            print '->',
            for rs in r.products:
                if len(rs.complexes) == 1:
                    print rs.complexes[0],
                else:
                    print 'rs' + str(rs),
            print

        # INCORRECT (misses resting state)
        for r in condensed['resting_states']:
            print 'r', r

        ## CORRECT (frozensets)
        for m in condensed['resting_state_map']:
            print 'm', m

        # CORRECT (all species)
        for t in condensed['resting_state_targets']:
            print 't', t

        # CORRECT (all resting states)
        for s in condensed['complexes_to_resting_states']:
            print 's', s


if __name__ == '__main__':
    unittest.main()
