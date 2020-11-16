#!/usr/bin/env python
#
#  test_objects.py
#  NuskellCompilerProject
#

import os
import sys
import unittest
from nuskell.dsdcompiler.objects import (SingletonError, 
                                         NuskellDomain, 
                                         NuskellComplex,
                                         clear_memory)

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestNuskellDomain(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_ComplementDomainInit(self):
        foo = NuskellDomain('foo', dtype = 'long')
        bar = ~foo
        moo = ~foo

        self.assertTrue(bar == moo, "bar is moo")
        self.assertTrue(bar is moo, "bar is moo")
        self.assertTrue(foo == ~bar, "foo is complement of bar")
        self.assertFalse(moo == ~bar, "moo is not complement of bar")
        self.assertFalse(foo == bar)
        self.assertTrue(bar.is_complement, "bar is complement")
        self.assertFalse(foo.is_complement, "foo is complement")

@unittest.skipIf(SKIP, "skipping tests")
class TestNuskellComplex(unittest.TestCase):
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

        self.assertEqual(list(foo.sequence), 
                   [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                    self.d1c, self.d3c, self.d1c, self.d2])
        self.assertEqual(list(foo.strand_table), 
                  [[self.d1, self.d2, self.d3], [self.d1], 
                   [self.d1c, self.d3c, self.d1c, self.d2]])

        try:
            bar = NuskellComplex([self.d1c, self.d3c, self.d1c, self.d2, '+', self.d1, self.d2, self.d3, '+', self.d1], list('((..+..)+)'))
        except SingletonError as err:
            bar = err.existing

        self.assertTrue(foo == bar)
        self.assertTrue(foo is bar)

    def test_names(self):
        NuskellComplex.PREFIX = 'cplx'
        NuskellComplex.ID = 0
        foo = NuskellComplex([self.d1, self.d2, self.d3, '+', self.d1, '+', 
                                self.d1c, self.d3c, self.d1c, self.d2], 
                             list('..(+(+))..'))
        assert foo.name == 'cplx0'
        assert foo.name != 'foo'

        with self.assertRaises(SingletonError):
            foo.name = 'foo'

        seq = list(foo.sequence)
        sst = list(foo.structure)
        del foo

        foo = NuskellComplex([self.d1, self.d2, self.d3, '+', self.d1, '+', 
                              self.d1c, self.d3c, self.d1c, self.d2], list('..(+(+))..'), name = 'foo')
        assert foo.name == 'foo'

    def test_rotations(self):
        d1, d2, d3 = self.d1, self.d2, self.d3
        d1c, d2c, d3c = self.d1c, self.d2c, self.d3c
        foo = NuskellComplex([d1, d2, d3, '+', d1, '+', d1c, d3c, d1c, d2], 
                structure = list('..(+(+))..'))

        seq1 = [d1,  d2,  d3,  '+', d1,  '+', d1c, d3c, d1c, d2] 
        sst1 = ['.', '.', '(', '+', '(', '+', ')', ')', '.', '.']
        seq2 = [d1,  '+', d1c, d3c, d1c, d2,  '+',d1,  d2,  d3,] 
        sst2 = ['(', '+', ')', '(', '.', '.', '+','.', '.', ')']
        seq3 = [d1c, d3c, d1c, d2,  '+',d1,  d2,  d3, '+', d1] 
        sst3 = ['(', '(', '.', '.', '+','.', '.', ')', '+', ')']

        for e, r in enumerate(foo.rotate()):
            if e == 0:
                assert r[0] == seq1
                assert r[1] == sst1
            elif e == 1:
                assert r[0] == seq2
                assert r[1] == sst2
            else:
                assert r[0] == seq3
                assert r[1] == sst3

if __name__ == '__main__':
    unittest.main()
