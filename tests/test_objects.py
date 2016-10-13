
import unittest
import nuskell.objects as objects

class DomainObjectTest(unittest.TestCase):
  # Test the Nuskell Version of DNAObjects
  def setUp(self):
    # preprocessing for unittesting
    pass 

  def tearDown(self):
    # clean up even if unittests failed
    pass

  #def testing(self):
  #  self.assertEqual(1,1, "test 1")
  #  self.assertTrue(True, "test true")
  #  self.assertFalse(False, "test false")
  #  with self.assertRaises(Warning):
  #    ...

  def test_DomainInit(self):
    with self.assertRaises(ValueError):
      doodle = objects.Domain()

    doodle = objects.Domain(name='doodle', constraints=list('Y'*5))

    self.assertIsInstance(doodle, objects.Domain, "doodle is a Domain")
    self.assertEqual(str(doodle), 'doodle', "print Domain")
    self.assertEqual(doodle.length, 5, "Domain length")
    self.assertEqual(doodle.sequence, list('Y'*5), "Domain sequence")

    moodle = objects.Domain(constraints=list('Y'*5))
    self.assertEqual(str(moodle), 'domain_{}'.format(moodle.id), "Automatic Domain Name")

    #self.assertTrue(False)

  def test_ComplemtDomainInit(self):
    foo = objects.Domain(constraints=list('Y'*5))

    # Conflicting Constraints
    with self.assertRaises(ValueError):
      bar = foo.get_ComplementDomain(constraints=list('R'*3))
    with self.assertRaises(ValueError):
      bar = foo.get_ComplementDomain(constraints=list('Y'*5))

    bar = foo.get_ComplementDomain(constraints=list('R'*5))
    self.assertEqual(foo.id, bar.id)

    moo = ~foo
    self.assertTrue(bar == moo, "bar is moo")
    self.assertTrue(foo == ~bar, "foo is complement of bar")
    self.assertFalse(moo == ~bar, "moo is not complement of bar")

  def test_NusDomain(self):
    nus = objects.NusDomain(constraints=list('Y'*5), domaintag='toehold')
    self.assertEqual(str(nus), 't{}'.format(nus.id), "NusDomain Name")

  def test_nucleotide_constraints(self):
    pass

  def test_domains_of_domains(self):
    #print 'm', moodle.base_length
    #domain.base_length
    #domain.base_sequence
    pass

  def test_complementarity(self):
    # foo is in ~bar
    pass

if __name__ == '__main__':
  unittest.main()

