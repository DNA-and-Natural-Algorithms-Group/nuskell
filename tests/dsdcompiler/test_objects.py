#!/usr/bin/env python
#
#  test_objects.py
#  NuskellCompilerProject
#

import os
import sys
import unittest
from nuskell.dsdcompiler.objects import NuskellDomain, clear_memory
from nuskell.dsdcompiler.objects import NuskellComplex


SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class NuskellDomainObjectTest(unittest.TestCase):
    def setUp(self):
        pass
    def tearDown(self):
        clear_memory()

    def test_ComplementDomainInit(self):
        foo = NuskellDomain('foo', dtype='long')
        bar = ~foo
        moo = ~foo

        self.assertTrue(bar == moo, "bar is moo")
        self.assertTrue(bar is moo, "bar is moo")
        self.assertTrue(foo == ~bar, "foo is complement of bar")
        self.assertFalse(moo == ~bar, "moo is not complement of bar")
        self.assertFalse(foo == bar)
        self.assertTrue(bar.is_complement, "bar is complement")
        self.assertFalse(foo.is_complement, "foo is complement")

    def test_domains_of_domains(self):
        d1aa = NuskellDomain('d1aa', dtype='short')
        d1ab = NuskellDomain('d1ab', dtype='short')

        with self.assertRaises(NotImplementedError):
            d1a = d1aa + d1ab


@unittest.skipIf(SKIP, "skipping tests")
class NuskellComplexObjectTest(unittest.TestCase):
    def setUp(self):
        self.d1 = NuskellDomain('d1', dtype='short')
        self.d2 = NuskellDomain('d2', dtype='short')
        self.d3 = NuskellDomain('d3', dtype='short')
        self.d1c = ~self.d1
        self.d2c = ~self.d2
        self.d3c = ~self.d3

    def tearDown(self):
        clear_memory()

    def test_ComplexInit(self):
        with self.assertRaises(TypeError):
            foo = NuskellComplex()

    def test_ComplexDomains(self):
        foo = NuskellComplex(
                sequence = [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                            self.d1c, self.d3c, self.d1c, self.d2], 
                structure = list('..(+(+))..'))

        self.assertEqual(foo.sequence, 
                   [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                    self.d1c, self.d3c, self.d1c, self.d2])
        self.assertEqual(foo.lol_sequence, 
                  [[self.d1, self.d2, self.d3], [self.d1], 
                   [self.d1c, self.d3c, self.d1c, self.d2]])

        bar = NuskellComplex(
                sequence = [self.d1c, self.d3c, self.d1c, self.d2, '+', 
                            self.d1, self.d2, self.d3, '+', self.d1], 
                structure = list('((..+..)+)'), memorycheck=False)

        self.assertEqual(foo, bar)
        self.assertTrue(foo == bar)
        self.assertFalse(foo is bar)

    def test_names(self):
        foo = NuskellComplex(
                sequence = [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                            self.d1c, self.d3c, self.d1c, self.d2], 
                structure = list('..(+(+))..'))
        self.assertEqual(foo.name, 'cplx0')

        foo.name = 'bar'
        foo.name = 'foo'
        self.assertEqual(foo.name, 'foo')
        self.assertSetEqual(set(NuskellComplex.NAMES.keys()), set(['foo']))

        foo.name = 'bar'
        self.assertEqual(foo.name, 'bar')
        self.assertSetEqual(set(NuskellComplex.NAMES.keys()), set(['bar']))

    def test_rotations(self):
        foo = NuskellComplex(
                sequence = [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                            self.d1c, self.d3c, self.d1c, self.d2], 
                structure = list('..(+(+))..'))
        self.assertEqual(foo.rotate_once(), foo)
        self.assertEqual(foo.sequence, 
                   [self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2, '+', 
                    self.d1, self.d2, self.d3])
        self.assertEqual(foo.structure, list('(+)(..+..)'))

        for r in foo.rotate():
            self.assertEqual(r, foo)
        self.assertEqual(foo.sequence, 
                   [self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2, '+', 
                    self.d1, self.d2, self.d3])
        self.assertEqual(foo.structure, list('(+)(..+..)'))

if __name__ == '__main__':
    unittest.main()
