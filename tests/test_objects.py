
import os
import sys
import unittest
import nuskell.objects as objects

class DomainObjectTest(unittest.TestCase):
  # Test the Nuskell Version of DNAObjects
  def setUp(self):
    # preprocessing for unittesting
    pass 

  def tearDown(self):
    # clean up even if unittests failed
    objects.Domain.id_counter = 0

  def test_DomainInit(self):
    doodle = objects.Domain(list('Y'*5), prefix='doodle')

    self.assertIsInstance(doodle, objects.Domain, "doodle is a Domain")
    self.assertEqual(str(doodle), '{}'.format(doodle.name), "print Domain")
    self.assertEqual(doodle.length, 5, "Domain length")
    self.assertEqual(doodle.sequence, list('Y'*5), "Domain sequence")

    moodle = objects.Domain(list('Y'*5))
    self.assertEqual(str(moodle), '{}'.format(moodle.name), 
        "Automatic Domain Name")

  def test_ComplementDomainInit(self):
    foo = objects.Domain(list('Y'*5))

    # Conflicting Constraints
    with self.assertRaises(ValueError):
      bar = foo.get_ComplementDomain(list('R'*3))
    with self.assertRaises(ValueError):
      bar = foo.get_ComplementDomain(list('Y'*5))
    with self.assertRaises(ValueError):
      foo.update_constraints(list('R'*6))

    bar = foo.get_ComplementDomain(list('R'*5))

    moo = ~foo
    self.assertTrue(bar == moo, "bar is moo")
    self.assertTrue(bar is moo, "bar is moo")
    self.assertTrue(foo == ~bar, "foo is complement of bar")
    self.assertFalse(moo == ~bar, "moo is not complement of bar")
    self.assertFalse(foo == bar)
    self.assertTrue(bar.is_ComplementDomain, "bar is complement")
    self.assertFalse(foo.is_ComplementDomain, "foo is complement")

  def test_domains_of_domains(self):
    d1aa = objects.Domain(list('N')) 
    d1ab = objects.Domain(list('Y')) 
    d1a = objects.Domain([d1aa, d1ab]) 

    d1b = objects.Domain(list('RR'))
    d1c = objects.Domain(list('NN'))
    d1 = objects.Domain([d1a,d1b,d1c])

    self.assertIsInstance(d1, objects.Domain, "d1 is a Domain")

    for d in d1.sequence :
      self.assertIsInstance(d, objects.Domain, "Sequence of Domains")

    self.assertEqual(d1.length, 3, "Length of Domain Sequence")
    self.assertListEqual(d1.base_sequence, list('NYRRNN'))
    self.assertEqual(d1.base_length, 6)

  def test_complement_domains_of_domains_of_domains(self):
    d1aa = objects.Domain(list('N')) 
    d1ab = objects.Domain(list('Y')) 
    d1a = objects.Domain([d1aa, d1ab]) 

    d1b = objects.Domain(list('RR'))
    d1c = objects.Domain(list('NN'))
    d1 = objects.Domain([d1a,d1b,d1c])
    with self.assertRaises(NotImplementedError):
      d2 = d1.get_ComplementDomain(list('R'*6))

    with self.assertRaises(NotImplementedError):
      d1.update_constraints(list('R'*6))

  def dont_test_complementarity(self):
    # TODO: This is not properly implemented yet.
    # e.g. foo is in ~bar
    pass

class ComplexObjectTest(unittest.TestCase):

  def setUp(self):
    self.d1 = objects.Domain(list('Y'*5))
    self.d2 = objects.Domain(list('N'*5))
    self.d3 = objects.Domain(list('R'*5))
    self.d1c = self.d1.get_ComplementDomain(list('N'*5))
    self.d2c = self.d2.get_ComplementDomain(list('D'*5))
    self.d3c = self.d3.get_ComplementDomain(list('H'*5))

  def tearDown(self):
    objects.Domain.id_counter = 0
    objects.Complex.id_counter = 0


  def test_ComplexInit(self):
    #NOTE: There is no particular reason for this Error, so it might change!
    with self.assertRaises(objects.NuskellObjectError):
      foo = objects.Complex()

    foo = objects.Complex(sequence=list('RNNNY'), structure=list('(...)'))
    self.assertIsInstance(foo, objects.Complex)

    self.assertEqual(foo.sequence, list('RNNNY'))
    self.assertEqual(foo.structure, list('(...)'))
    self.assertEqual(foo.lol_sequence, [list('RNNNY')])
    self.assertEqual(foo.nucleotide_sequence, list('RNNNY'))
    self.assertEqual(foo.rotate_once, foo)
    for r in foo.rotate:
      self.assertEqual(r.sequence, list('RNNNY'))
      self.assertEqual(r.structure, list('(...)'))

  def test_ComplexDomains(self):
    foo = objects.Complex(sequence=[self.d1, self.d2, self.d3, '+', self.d1,
      '+', self.d1c, self.d3c, self.d1c, self.d2], structure=list('..(+(+))..'))

    self.assertEqual(foo.sequence, [self.d1, self.d2, self.d3, '+', self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2])
    self.assertEqual(foo.lol_sequence, [[self.d1, self.d2, self.d3], [self.d1], [self.d1c, self.d3c, self.d1c, self.d2]])
    self.assertEqual(foo.nucleotide_sequence, list('YYYYYNNNNNRRRRR+YYYYY+RRRRRYYYYYRRRRRNNNNN'))

    bar = objects.Complex(sequence=[self.d1c, self.d3c, self.d1c, self.d2, '+',
      self.d1, self.d2, self.d3, '+', self.d1], structure=list('((..+..)+)'))

    self.assertEqual(foo, bar)
    self.assertTrue(foo == bar)

  def test_rotations(self):
    foo = objects.Complex(sequence=[self.d1, self.d2, self.d3, '+', self.d1,
      '+', self.d1c, self.d3c, self.d1c, self.d2], structure=list('..(+(+))..'))
    self.assertEqual(foo.rotate_once, foo)
    self.assertEqual(foo.sequence, [self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2, '+', self.d1, self.d2, self.d3])
    self.assertEqual(foo.structure,list('(+)(..+..)'))
    self.assertEqual(foo.nucleotide_sequence, list('YYYYY+RRRRRYYYYYRRRRRNNNNN+YYYYYNNNNNRRRRR'))
    for r in foo.rotate :
      self.assertEqual(r, foo)
    self.assertEqual(foo.sequence, [self.d1, '+', self.d1c, self.d3c, self.d1c, self.d2, '+', self.d1, self.d2, self.d3])
    self.assertEqual(foo.structure,list('(+)(..+..)'))
    self.assertEqual(foo.nucleotide_sequence, list('YYYYY+RRRRRYYYYYRRRRRNNNNN+YYYYYNNNNNRRRRR'))
 
class TestTubeIOTest(unittest.TestCase):
  def setUp(self):
    self.t0 = objects.Domain(list('N'*5), prefix='t')
    self.t1 = objects.Domain(list('N'*5), prefix='t')
    self.t2 = objects.Domain(list('N'*5), prefix='t')
    self.d3 = objects.Domain(list('N'*15), prefix='d')
    self.d4 = objects.Domain(list('N'*15), prefix='d')
    self.d5 = objects.Domain(list('N'*15), prefix='d')
    self.d6 = objects.Domain(list('N'*15), prefix='d')
    self.d7 = objects.Domain(list('N'*15), prefix='d')
    self.t0c = self.t0.get_ComplementDomain(list('N'*5))
    self.t1c = self.t1.get_ComplementDomain(list('N'*5))
    self.t2c = self.t2.get_ComplementDomain(list('N'*5))
    self.d3c = self.d3.get_ComplementDomain(list('N'*15))
    self.d4c = self.d4.get_ComplementDomain(list('N'*15))
    self.d5c = self.d5.get_ComplementDomain(list('N'*15))
    self.d6c = self.d6.get_ComplementDomain(list('N'*15))
    self.d7c = self.d7.get_ComplementDomain(list('N'*15))

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

    sequence = [d4, t0, '+', d6, t2,     '+', d3,  '+', d3, t2,   '+', t0, d5, t1,    '+', t0, d7, t1,    '+', t1c, d7c, t1c, d5c, t2c, '+', d3c, d3c, t2c, d6c, t0c, d4c, t0c]
    structure =['(', '(', '+', '(', '(', '+', '(', '+', '(', '(', '+', '.', '(', '(', '+', '.', '(', '(', '+', ')', ')', ')', ')', ')', '+', ')', ')', ')', ')', ')', ')', '.']
        
    foo = objects.Complex(sequence=sequence, structure = structure)
    fooIO = objects.TestTubeIO(objects.TestTube(complexes={foo.name: (foo, float("inf"), None)}))

    #fooIO.write_dnafile(sys.stdout)
    f = open(os.devnull,"w")
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

    sequence = [d4, t0, '+', d6, t2,     '+', d3,  '+', d3, t2,   '+', t0, d5, t1,    '+', t0, d7, t1,    '+', t1c, d7c, t1c, d5c, t2c, '+', d3c, d3c, t2c, d6c, t0c, d4c, t0c]
    structure =['(', '(', '+', '(', '(', '+', '(', '+', '(', '(', '+', '.', '(', '(', '+', '.', '(', '(', '+', ')', ')', ')', ')', ')', '+', ')', ')', ')', ')', ')', ')', '.']

    foo = objects.Complex(sequence=sequence, structure = structure)
    fooIO = objects.TestTubeIO(objects.TestTube(complexes={foo.name: (foo, float("inf"), None)}))

    f = open(os.devnull,"w")
    fooIO.write_pil_kernel(f)

if __name__ == '__main__':
  unittest.main()

