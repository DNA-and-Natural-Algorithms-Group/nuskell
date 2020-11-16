#
# Unittests for nuskell.dsdcompiler.compiler
#
# Written by Stefan Badelt (bad-ants-fleet@posteo.eu).
#

import unittest
from nuskell.dsdcompiler.objects import clear_memory
from nuskell.dsdcompiler.compiler import translate

class Test_Workflow(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_compile(self):
        crn = 'A + B -> C + D; A + A <=> C + A; C @i 5'
        ts = 'soloveichik2010.ts'

        solution, modules = translate(crn, ts, modular = True)

        f = [x for x in solution if x[0] == 'f']
        s = [x for x in solution if x[0] != 'f']

        assert len(s) == 4
        assert len(f) == 9
        assert len(modules) == 2
        assert solution['A'].concentration is None
        assert solution['C'].concentration[0] == 'initial'
        assert solution['C'].concentration[1] == 5
        assert solution['C'].concentration[2] == 'nM'
        assert solution['f1'].concentration[0] == 'constant'
        assert solution['f1'].concentration[1] == 100
        assert solution['f1'].concentration[2] == 'nM'


if __name__ == '__main__':
    unittest.main()
