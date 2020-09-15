#
# Unittests for nuskell.dsdcompiler.compiler
#
# Written by Stefan Badelt (bad-ants-fleet@posteo.eu).
#

import unittest
from nuskell.dsdcompiler.compiler import translate

class Test_Workflow(unittest.TestCase):

    def test_compile(self):
        crn = 'A + B -> C + D'
        ts = 'soloveichik2010.ts'

        a, b = translate(crn, ts)

        print(a)
        print(b)

if __name__ == '__main__':
    unittest.main()
