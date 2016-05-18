import unittest
import nuskell.parser as np


class TestTLSparser(unittest.TestCase):

  def setUp(self):
    """ test individual header functions """
    self.f1 = "function len(x) = if x == [] then 0 else 1 + len(tail(x))"
    self.f2 = "function range(x) = if x == 0 then [] else range(x - 1) + [x - 1]"
    self.f3 = "function sum(x) = if len(x) == 0 then empty elseif len(x) == 1 then x[0] else x[0] + sum(tail(x))"
    self.f4 = "function reverse(x) = if x == [] then [] else reverse(tail(x)) + [x[0]] "
    self.f5 = "function rxn_degree(x, r) = if len(x) == 0 then [] elseif len(x[0].reactants) == r then [x[0]] + rxn_degree(tail(x), r) else rxn_degree(tail(x), r)"
    self.f6 = "function unirxn(x) = if len(x) == 0 then [] elseif len(x[0].reactants) == 1 then [x[0]] + unirxn(tail(x)) else unirxn(tail(x)) "
    self.f7 = "function birxn(x) = if len(x) == 0 then [] elseif len(x[0].reactants) == 2 then [x[0]] + birxn(tail(x)) else birxn(tail(x)) "
    self.f8 = "function map(f, x) = if len(x) == 0 then [] else [f(x[0])] + map(f, tail(x)) "
    self.f9 = "function map2(f, y, x) = if len(x) == 0 then [] else [f(y, x[0])] + map2(f, y, tail(x)) "

  def tearDown(self):
    """ clean up after unittest (even if not successfull), 
    removing this causes very weird behavior... 
    """
    pass

  def test_tls_parse_data(self):
    #print self.f1

    #"""
    #[[function] [len(x)] = [if] [[x] [==] [[]]] [then] [0] else 1 + len(tail(x))]
    #"""
    pf1 = [
        ['function', 
          ['id', 'len'], 
          [['id', 'x']], 
          ['if', 

            ['where', 
              ['==', 
                ['trailer', 
                  ['id', 'x']], 
                ['trailer', 
                  ['list']]]], 

            ['where', 
              ['trailer', 
                ['num', '0']]], 

            ['where', 
              ['+', 
                ['trailer', 
                  ['num', '1']], 
                ['trailer', 
                  ['id', 'len'], 
                  ['apply', 
                    ['where', 
                      ['trailer', 
                        ['id', 'tail'], 
                        ['apply', 
                          ['where', 
                            ['trailer', 
                              ['id', 'x']]]]]]]]]]]]]

    self.assertEqual(
        np.ts_parser.parse(self.f1), pf1, 'built-in function len(x)')

    #print self.f2
    #print np.ts_parser.parse(self.f2)
    pf2 = [['function', ['id', 'range'], [['id', 'x']], 
      ['if', ['where', ['==', ['trailer', ['id', 'x']], ['trailer', ['num', '0']]]], 
             ['where', ['trailer', ['list']]], 
             ['where', ['+', ['trailer', ['id', 'range'], ['apply', ['where', ['-', ['trailer', ['id', 'x']], ['trailer', ['num', '1']]]]]], ['trailer', ['list', ['where', ['-', ['trailer', ['id', 'x']], ['trailer', ['num', '1']]]]]]]]]]]
    self.assertEqual(
        np.ts_parser.parse(self.f2), pf2, 'built-in function range(x)')

if __name__ == '__main__':
  unittest.main()

