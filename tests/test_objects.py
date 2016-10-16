
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

  def test_DomainInit(self):
    with self.assertRaises(ValueError):
      doodle = objects.Domain()

    doodle = objects.Domain(name='doodle', constraints=list('Y'*5))

    self.assertIsInstance(doodle, objects.Domain, "doodle is a Domain")
    self.assertEqual(str(doodle), 'doodle', "print Domain")
    self.assertEqual(doodle.length, 5, "Domain length")
    self.assertEqual(doodle.sequence, list('Y'*5), "Domain sequence")

    moodle = objects.Domain(constraints=list('Y'*5))
    self.assertEqual(str(moodle), 'domain_{}'.format(moodle.id), 
        "Automatic Domain Name")

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
    self.assertFalse(foo == bar)

  def test_NusDomain(self):
    nus = objects.NusDomain(constraints=list('Y'*5), domaintag='toehold')
    self.assertEqual(str(nus), 't{}'.format(nus.id), "NusDomain Name")

    cus = nus.get_ComplementDomain(constraints=list('R'*5))
    self.assertEqual(cus.id, nus.id)

    mus = ~cus
    self.assertTrue(mus == nus)
    self.assertTrue(mus == ~cus)
    self.assertTrue(nus == ~cus)
    # The following line revealed an unexpected bug.
    self.assertFalse(nus == ~mus)

  def test_domains_of_domains(self):
    d1aa = objects.Domain(constraints='N') 
    d1ab = objects.Domain(constraints='Y') 

    d1a = objects.Domain(subdomains=[d1aa, d1ab]) 
    d1b = objects.Domain(constraints='RR')
    d1c = objects.Domain(constraints='NN')

    d1 = objects.Domain(subdomains=[d1a,d1b,d1c])

    self.assertIsInstance(d1, objects.Domain, "d1 is a Domain")

    for d in d1.sequence :
      self.assertIsInstance(d, objects.Domain, "Sequence of Domains")

    self.assertEqual(d1.length, 3, "Length of Domain Sequence")
    self.assertListEqual(d1.base_sequence, list('NYRRNN'))
    self.assertEqual(d1.base_length, 6)

  def test_nucleotide_constraints(self):
    pass

  def test_complementarity(self):
    # TODO: This is not properly implemented yet.
    # e.g. foo is in ~bar
    pass


if __name__ == '__main__':
  unittest.main()

