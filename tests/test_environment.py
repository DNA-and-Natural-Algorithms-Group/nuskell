import unittest
import nuskell.interpreter as ni
from nuskell.interpreter.environment import *

class TestReaction(unittest.TestCase):
  def setUp(self):
    self.r1 = Reaction(['A', 'B', 'C'], ['X'], True)
    self.r2 = Reaction([], ['B'], False)
    self.r3 = Reaction([], ['B'], True)

  def test_assignments(self):
    self.assertEqual(self.r1.reactants, ['A', 'B', 'C'], 'test_reactants 1')
    self.assertEqual(self.r2.reactants, [], 'test_reactants 2')
    self.assertEqual(self.r3.products, ['B'], 'test_products')
    self.assertEqual(self.r3.reversible, True, 'test_reversibility')


class TestDomain(unittest.TestCase):
  def setUp(self):
    # TODO: a domain with 0 length? 
    self.d0 = Domain(0, id=-25)
    self.d1 = Domain(8)
    self.d2 = Domain(4)
    self.d3 = Domain(14, id=1)

  def tearDown(self):
    """ clean up after unittest (even if not successfull), 
    removing this causes very weird behavior... 
    """
    Domain.domain_id = 0

  def test_assignments(self):
    self.assertEqual(self.d0.id, -25,'test_domain_id 0')
    self.assertEqual(self.d0.id, -25,'test_domain_id 0')
    self.assertEqual(self.d1.id, 1, 'test_domain_id 1')
    self.assertEqual(self.d2.id, 2, 'test_domain_id 2')
    self.assertEqual(self.d3.id, 1, 'test_domain_id 3')

    self.assertEqual(str(self.d1), 'd1', 'test_str 1')
    self.assertEqual(self.d1.length, 8, 'test_length 1')
    self.assertEqual(str(self.d3), 'd1', 'test_str 2')
    self.assertEqual(self.d3.length, 14, 'test_length 2')

  def test_weird_behavior(self):
    # potential misusage test-case
    print 'Xd1', self.d1, self.d1.length, self.d1.id
    print 'Xd3', self.d3, self.d3.length, self.d3.id
    self.assertEqual(self.d1 == self.d3, True, 'test_is_equal 1')

class TestEnvironment(unittest.TestCase):
  def setUp(self):
    self.env1 = Environment('env1')

  def tearDown(self):
    pass

  def test_base_level_functions(self):
    # Feel free to add base-level functions here
    self.assertEqual(len(self.env1.env), 1, 'in baselevel')
    baselevel = {'tail', 'complement', 'infty', 'short', 'long', 'unique',
        'flip', 'empty', 'rev_reactions', 'irrev_reactions'}

    self.assertEqual(set(self.env1.env[0].keys()), baselevel, 'in baselevel')

  def test_create_binding(self):
    self.env1._create_binding('testfunc', Function('ar','body'))
    self.assertEqual('testfunc' in self.env1.env[-1].keys(), True, 'in baselevel')

    #self.env1._create_binding('tail', Function('ar','value'))

class Test_funky(unittest.TestCase):
  def setUp(self):
    self.env1 = Environment('env1')

    self.t1 = Domain(5)

    domains = map(Domain, [5, 6, 7])
    dotparens = ['.', '(', ')']
    attributes = { }
    self.s1 = Structure(domains, dotparens, attributes)

  def tearDown(self):
    Domain.domain_id = 0

  def test_builtin_tail(self):
    l = [1,2,3,4]
    self.assertEqual(funky()._tail([l]), [2,3,4], 'tail testing')

  def test_builtin_complement(self):
    t1 = self.t1
    t2 = self.s1
    t3 = '('
    t4 = ')'
    t5 = 'something'
    tl = [t1, t2, t3, t4, t5]

    self.assertEqual(funky()._complement([t1]), Domain(5, id=-1), 
        'complement testing')

    # super annoying test ... abusing the __str__ function of Domains
    r2 = funky()._complement([t2])
    r2dots = list(r2.dotparens)
    self.assertEqual(r2dots, ['(',')','.'], 'complement testing 2')

    r2doml = list(r2.domains)
    r2dom1 = r2doml[0]
    r2dom2 = r2doml[1]
    r2dom3 = r2doml[2]
    r2doms = [r2dom1, r2dom2, r2dom3]

    x2dom1 = Domain(7,id=-4)
    x2dom2 = Domain(6,id=-3) 
    x2dom3 = Domain(5,id=-2)
    x2doms = [x2dom1, x2dom2, x2dom3]

    self.assertEqual(r2doms, x2doms, 'complement testing 2')

    self.assertEqual(funky()._complement([t3]), ')', 'complement testing')
    self.assertEqual(funky()._complement([t4]), '(', 'complement testing')
    self.assertEqual(funky()._complement([t5]), 'something', 'complement testing')
    
    rl = list(funky()._complement([tl]))
    rlt = [t5, t3, t4]
    self.assertEqual(rl[:3], rlt, 'complement testing')
    self.assertEqual(isinstance(rl[3], Structure), True, 'complement testing')
    self.assertEqual(isinstance(rl[4], Domain), True, 'complement testing')

  def test_flip(self):
    n = 3 # dimenstion of inner list
    l = [['a','b','c'], ['d','e','f'], ['g','h','i'], ['j','k','l']]
    r = [['a', 'd', 'g', 'j'], ['b', 'e', 'h', 'k'], ['c','f','i','l']]

    self.assertEqual(funky()._flip([l,n]), r, 'flip testing')
    print funky()._flip([l, n])


  def test_if(self):
    """ not ready yet """
    #pass
    #content = [
    #    ['where', 
    #      ['==', 
    #        ['trailer', 
    #          ['id', 'x']], 
    #        ['trailer', ['list']]]], 
    #      ['where', 
    #        ['trailer', 
    #          ['num', '0']]], 
    #      ['where', 
    #        ['+', 
    #          ['trailer', 
    #            ['num', '1']], 
    #          ['trailer', 
    #            ['id', 'len'], 
    #            ['apply', 
    #              ['where', 
    #                ['trailer', 
    #                  ['id', 'tail'], 
    #                  ['apply', 
    #                    ['where', 
    #                      ['trailer', 
    #                        ['id', 'x']]]]]]]]]]]
    #print funky()._if(self.env1, content)
    #self.env1._create_binding('tail', Function('ar','value'))

