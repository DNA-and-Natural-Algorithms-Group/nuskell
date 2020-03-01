#!/usr/bin/env python
#
#  test_objects.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function

import os
import sys
import unittest
import nuskell.objects as objects
from collections import Counter

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class NuskellDomainObjectTest(unittest.TestCase):
    # Test the Nuskell Version of DSDObjects
    def setUp(self):
        pass

    def tearDown(self):
        objects.clear_memory()

    def test_DomainInit(self):
        pass

    def test_ComplementDomainInit(self):
        foo = objects.NuskellDomain('foo', dtype='long')
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
        d1aa = objects.NuskellDomain('d1aa', dtype='short')
        d1ab = objects.NuskellDomain('d1ab', dtype='short')

        with self.assertRaises(NotImplementedError):
            d1a = d1aa + d1ab


@unittest.skipIf(SKIP, "skipping tests")
class NuskellComplexObjectTest(unittest.TestCase):

    def setUp(self):
        self.d1 = objects.NuskellDomain('d1', dtype='short')
        self.d2 = objects.NuskellDomain('d2', dtype='short')
        self.d3 = objects.NuskellDomain('d3', dtype='short')
        self.d1c = ~self.d1
        self.d2c = ~self.d2
        self.d3c = ~self.d3

    def tearDown(self):
        objects.clear_memory()

    def test_ComplexInit(self):
        # NOTE: There is no particular reason for this Error, so it might
        # change!
        with self.assertRaises(TypeError):
            foo = objects.NuskellComplex()

    def test_ComplexDomains(self):
        foo = objects.NuskellComplex(
                sequence = [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                            self.d1c, self.d3c, self.d1c, self.d2], 
                structure = list('..(+(+))..'))

        self.assertEqual(foo.sequence, 
                   [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                    self.d1c, self.d3c, self.d1c, self.d2])
        self.assertEqual(foo.lol_sequence, 
                  [[self.d1, self.d2, self.d3], [self.d1], 
                   [self.d1c, self.d3c, self.d1c, self.d2]])

        bar = objects.NuskellComplex(
                sequence = [self.d1c, self.d3c, self.d1c, self.d2, '+', 
                            self.d1, self.d2, self.d3, '+', self.d1], 
                structure = list('((..+..)+)'), memorycheck=False)

        self.assertEqual(foo, bar)
        self.assertTrue(foo == bar)
        self.assertFalse(foo is bar)

    def test_names(self):
        foo = objects.NuskellComplex(
                sequence = [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                            self.d1c, self.d3c, self.d1c, self.d2], 
                structure = list('..(+(+))..'))
        self.assertEqual(foo.name, 'cplx0')

        foo.name = 'bar'
        foo.name = 'foo'
        self.assertEqual(foo.name, 'foo')
        self.assertSetEqual(set(objects.NuskellComplex.NAMES.keys()), set(['foo']))

        foo.name = 'bar'
        self.assertEqual(foo.name, 'bar')
        self.assertSetEqual(set(objects.NuskellComplex.NAMES.keys()), set(['bar']))

    def test_rotations(self):
        foo = objects.NuskellComplex(
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


@unittest.skipIf(SKIP, "skipping tests")
class TestTubeTests(unittest.TestCase):
    def setUp(self):
        self.h1 = objects.NuskellDomain('h1', dtype='long')
        self.d1 = objects.NuskellDomain('d1', dtype='long')
        self.d2 = objects.NuskellDomain('d2', dtype='short')
        self.d3 = objects.NuskellDomain('d3', dtype='short')
        self.d1c = ~self.d1
        self.d2c = ~self.d2
        self.d3c = ~self.d3

        self.cplx1 = objects.NuskellComplex(
                [self.d1, self.d2, self.d3, '+', self.d1, '+', 
                 self.d1c, self.d3c, self.d1c, self.d2], 
                list('..(+(+))..'), name = 'cplx1')

        self.cplx2 = objects.NuskellComplex(
                [self.d1, self.d2, self.d3, '+', self.d3c, self.d2c, self.d2], 
                list('.((+)).'), name = 'cplx2')

        self.cplx3 = objects.NuskellComplex(
                [self.h1, self.d2, self.d3, '+', self.d3c, self.d2c, self.d2], 
                list('.((+)).'), name = 'cplx3')

    def tearDown(self):
        objects.clear_memory()

    def test_TestTubeInit(self):
        with self.assertRaises(AssertionError):
            complexes = { 
                    'name1' : [self.cplx1], 
                    'name2' : [self.cplx2] }
            foo = objects.TestTube(complexes)

        complexes = { 
                'name1' : [self.cplx1, None, None, None], 
                'name2' : [self.cplx2, None, None, None] }
        foo = objects.TestTube(complexes)

        complexes = { 
                'name1' : [self.cplx1, 0.1, True, 'fuel'], 
                'name2' : [self.cplx2, 1e-9, False, 'signal'] }
        foo = objects.TestTube(complexes)

        complexes = [self.cplx1, self.cplx2]
        foo = objects.TestTube(complexes)

        foo = objects.TestTube()
        foo.add_complex(self.cplx1, (None, None))
        self.assertTrue(foo.has_complex(self.cplx1))

    def test_TestTubeInterpret(self):
        complexes = [self.cplx1, self.cplx2, self.cplx3]
        foo = objects.TestTube(complexes)
        out = foo.interpret_species([self.cplx3.name], prune = False)
        exp = {'cplx3_1_' : Counter({'cplx3':1})}
        self.assertDictEqual(exp, out)
        self.assertEqual(sorted(foo.complexes), sorted([self.cplx1, self.cplx2]))

        objects.clear_memory()
        complexes = [self.cplx1, self.cplx2, self.cplx3]
        foo = objects.TestTube(complexes)
        out = foo.interpret_species([self.cplx3.name], prune = True)
        exp = {'cplx3_1_' : Counter({'cplx3':1})}
        self.assertDictEqual(exp, out)
        self.assertEqual(sorted(foo.complexes), sorted([self.cplx2]))

    def test_TestTubeSum(self):
        foo = objects.TestTube()
        bar = objects.TestTube()
        foo.add_complex(self.cplx1, (None, True))
        bar.add_complex(self.cplx2, (float('inf'), None))

        # Both versions should work
        foobar = foo + bar
        foobar = sum([foo, bar])

        # Check if original TestTubes remain unchanged
        self.assertTrue(foo.has_complex(self.cplx1))
        self.assertFalse(foo.has_complex(self.cplx2))

        self.assertTrue(bar.has_complex(self.cplx2))
        self.assertFalse(bar.has_complex(self.cplx1))

        # Check if new TestTube has everything
        self.assertTrue(foobar.has_complex(self.cplx1))
        self.assertTrue(foobar.has_complex(self.cplx2))

        # Check if attributes were copied correctly
        self.assertEqual(
            foo.get_complex_concentration(
                self.cplx1), foobar.get_complex_concentration(
                self.cplx1))
        self.assertEqual(
            bar.get_complex_concentration(
                self.cplx2), foobar.get_complex_concentration(
                self.cplx2))


@unittest.skipIf(SKIP, "skipping tests")
class TestTubeIOTest(unittest.TestCase):
    def setUp(self):
        self.t0 = objects.NuskellDomain(dtype='short', prefix='t')
        self.t1 = objects.NuskellDomain(dtype='short', prefix='t')
        self.t2 = objects.NuskellDomain(dtype='short', prefix='t')
        self.d3 = objects.NuskellDomain(dtype='long', prefix='d')
        self.d4 = objects.NuskellDomain(dtype='long', prefix='d')
        self.d5 = objects.NuskellDomain(dtype='long', prefix='d')
        self.d6 = objects.NuskellDomain(dtype='long', prefix='d')
        self.d7 = objects.NuskellDomain(dtype='long', prefix='d')
        self.t0c = ~self.t0
        self.t1c = ~self.t1
        self.t2c = ~self.t2
        self.d3c = ~self.d3
        self.d4c = ~self.d4
        self.d5c = ~self.d5
        self.d6c = ~self.d6
        self.d7c = ~self.d7

    def tearDown(self):
        objects.clear_memory()

    def test_IO_dna(self):
        # NOTE: This function writes to a file, so we need to compare it to
        # pre-written output-file. Doesn't feel necessary at this point ...
        t0 = self.t0
        t1 = self.t1
        t2 = self.t2
        d3 = self.d3
        d4 = self.d4
        d5 = self.d5
        d6 = self.d6
        d7 = self.d7
        t0c = self.t0c
        t1c = self.t1c
        t2c = self.t2c
        d3c = self.d3c
        d4c = self.d4c
        d5c = self.d5c
        d6c = self.d6c
        d7c = self.d7c

        sequence =  [ d4,  t0, '+',  d6,  t2, '+',  d3, '+',  d3,  t2, '+',  t0,  d5,  t1, '+',  t0,  d7,  t1, '+', t1c, d7c, t1c, d5c, t2c, '+', d3c, d3c, t2c, d6c, t0c, d4c, t0c]
        structure = ['(', '(', '+', '(', '(', '+', '(', '+', '(', '(', '+', '.', '(', '(', '+', '.', '(', '(', '+', ')', ')', ')', ')', ')', '+', ')', ')', ')', ')', ')', ')', '.']

        foo = objects.NuskellComplex(sequence=sequence, structure=structure)
        fooIO = objects.TestTubeIO(objects.TestTube(complexes = {foo.name: [foo, float("inf"), None, None]}))

        # fooIO.write_dnafile(sys.stdout)
        f = open(os.devnull, "w")
        fooIO.write_dnafile(f)

    def test_IO_kernel(self):
        t0 = self.t0
        t1 = self.t1
        t2 = self.t2
        d3 = self.d3
        d4 = self.d4
        d5 = self.d5
        d6 = self.d6
        d7 = self.d7
        t0c = self.t0c
        t1c = self.t1c
        t2c = self.t2c
        d3c = self.d3c
        d4c = self.d4c
        d5c = self.d5c
        d6c = self.d6c
        d7c = self.d7c

        sequence = [ d4, t0, '+', d6, t2, '+', d3, '+', d3, t2, '+', t0, d5, t1, '+', t0, d7, t1, '+', t1c, d7c, t1c, d5c, t2c, '+', d3c, d3c, t2c, d6c, t0c, d4c, t0c]
        structure = [ '(', '(', '+', '(', '(', '+', '(', '+', '(', '(', '+', '.', '(', '(', '+', '.', '(', '(', '+', ')', ')', ')', ')', ')', '+', ')', ')', ')', ')', ')', ')', '.']

        foo = objects.NuskellComplex(sequence=sequence, structure=structure)
        fooIO = objects.TestTubeIO( objects.TestTube( complexes = {foo.name: [foo, float("inf"), None, 'signal']}))

        f = open(os.devnull, "w")
        fooIO.write_pil(f)


if __name__ == '__main__':
    unittest.main()
