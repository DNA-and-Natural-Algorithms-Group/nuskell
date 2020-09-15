#
# Unittests for nuskell.dsdcompiler.compiler
#
# Written by Stefan Badelt (bad-ants-fleet@posteo.eu).
#

import unittest
from nuskell.dsdcompiler.compiler import translate

class Test_Workflow(unittest.TestCase):
    def test_compile(self):
        crn = 'A + B -> C + D; A + A <=> C + A; C @i 5'
        ts = 'schemes/literature/soloveichik2010.ts'

        solution, modules = translate(crn, ts, modular = True)

        f = [x for x in solution if x[0] == 'f']
        s = [x for x in solution if x[0] != 'f']

        assert len(s) == 4
        assert len(f) == 9
        assert len(modules) == 2

        assert solution['A'].concentration is None
        assert solution['C'].concentration.mode == 'initial'
        assert solution['C'].concentration.value == 5
        assert solution['C'].concentration.unit == 'nM'
        assert solution['f1'].concentration.mode == 'constant'
        assert solution['f1'].concentration.value == 100
        assert solution['f1'].concentration.unit == 'nM'


if __name__ == '__main__':
    unittest.main()
