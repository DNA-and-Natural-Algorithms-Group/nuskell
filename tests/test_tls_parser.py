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

    # Figure out parsing for global statement:

    #self.m5 = "global toehold = short()" # parsed correctly!
    #self.m5 = "global [toehold, franky, fish] = short()" # parsed correctly, but may not be desired?
    self.m5 = "global [toehold, franky, fish] = [short(), long(), short()]" # parsed correctly, but may not be desired?
    #self.m5 = "global [toehold = short(); franky=long(); fish=short()]" # not paresd!

    self.m4 = """class bgate(s1, s2, l) = [ "b c d + e f g + f* e* d* c* b* a*" | "( ( ( + ( ( ~ + )  )  )  )  )  .", "h + i f*" | "~ + ~ .", "b c d" | ". . ." ]
    where { a = s1.a ; b = s1.b ; c = s1.c ; d = s2.a ; e = s2.b ; f = s2.c ; [g, h, j] = flip(map(gmac, l), 3) ; i = reverse(j) } """

    #self.m3 = "module rxn(r) = if len(r.reactants) == 1 then unimolecular(r) elseif len(r.reactants) == 2 then bimolecular(r) else print('warn me')"
    self.m3 = "module rxn(r) = if len(r.reactants) == 1 then unimolecular(r) where {void = print('y')} else print('warn me')"
    self.m2 = "module bimolecular(r) = infty(l) + infty(t) + infty(b) where [l, t, b] = bgate(r.reactants[0], r.reactants[1], r.products)"
    self.m1 = "module main(crn) = sum(map(rxn, crn)) where crn = irrev_reactions(crn)"

  def tearDown(self):
    """ clean up after unittest (even if not successfull), 
    removing this causes very weird behavior... 
    """
    pass

  def test_tls_parse_data(self):

    #"""
    #print self.f1
    #[[function] [len(x)] = if x == [] then 0 else 1 + len(tail(x))]
    #"""

    # [type=function, [name=[id:len], args=[[id:x]], body=[...]]]

    "[[function [len, [x]] = [if] [[x] [==] [[]]] [then] [0] else 1 + len(tail(x))]"
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
    pf2 = [['function', 
      ['id', 'range'], [['id', 'x']], ['if', 
        ['where', ['==',
          ['trailer', ['id', 'x']], 
          ['trailer', ['num', '0']]
        ]], 
        ['where', ['trailer', ['list']]], 
        ['where', ['+', 
          ['trailer', ['id', 'range'], 
            ['apply', ['where', ['-', ['trailer', ['id', 'x']], ['trailer',
                ['num', '1']]]]]], 
          ['trailer', ['list', ['where', ['-', ['trailer',
            ['id', 'x']], ['trailer', ['num', '1']]]]]]]]]]]

    self.assertEqual(np.ts_parser.parse(self.f2), pf2, 'built-in function range(x)')


    pf3 = [['function', 
      ['id', 'sum'], [['id', 'x']], 
      ['if', 
        ['where', ['==', 
          ['trailer', ['id', 'len'], ['apply', ['where', ['trailer', ['id', 'x']]]]], 
          ['trailer', ['num', '0']]]], 
        ['where', ['trailer', ['id', 'empty']]], 
        ['where', ['==', 
          ['trailer', ['id', 'len'], ['apply', ['where', ['trailer', ['id', 'x']]]]], 
          ['trailer', ['num', '1']]]], 
        ['where', 
          ['trailer', ['id', 'x'], 
            ['index', 
              ['where', ['trailer', ['num', '0']]]]]], 
        ['where', ['+', 
          ['trailer', ['id', 'x'], ['index', ['where', ['trailer', ['num', '0']]]]], 
          ['trailer', ['id', 'sum'], 
            ['apply', 
              ['where', 
                ['trailer', 
                  ['id', 'tail'], 
                  ['apply', 
                    ['where', 
                      ['trailer', 
                        ['id', 'x']]]]]]]]]]]]]

    self.assertEqual(np.ts_parser.parse(self.f3), pf3, 'built-in function sum(x)')


    print self.m1
    pm1 = [['module', 
      ['id', 'main'], 
      [['id', 'crn']], 
      ['where', 
        ['trailer', ['id', 'sum'], 
          ['apply', ['where', ['trailer', ['id', 'map'], 
          ['apply', ['where', ['trailer', ['id', 'rxn']]], 
            ['where', ['trailer', ['id', 'crn']]]]]]]], 
        [[['id', 'crn'], 
          ['where', ['trailer', ['id', 'irrev_reactions'], 
            ['apply', 
              ['where', ['trailer', ['id', 'crn']]]]]]]]]]]

    print self.m3
    print np.ts_parser.parse(self.m3)
    pm3 = [['module', 
      ['id', 'rxn'], 
      [['id', 'r']], 
      ['if', 
        ['where', ['==', 
          ['trailer', ['id', 'len'], ['apply', ['where', ['trailer', ['id', 'r'], ['attribute', ['id', 'reactants']]]]]], 
          ['trailer', ['num', '1']]]], 
        ['where', ['trailer', ['id', 'unimolecular'], ['apply', ['where', ['trailer', ['id', 'r']]]]]], 

        ['where', ['==', 
          ['trailer', ['id', 'len'], ['apply', ['where', ['trailer', ['id', 'r'], ['attribute', ['id', 'reactants']]]]]], 
          ['trailer', ['num', '2']]]], 
        ['where', ['trailer', ['id', 'bimolecular'], ['apply', ['where', ['trailer', ['id', 'r']]]]]], 

        ['where', ['trailer', ['id', 'r'], ['index', ['where', ['trailer', ['num', '0']]]]]]]]]

    pm3 = [['module', 
      ['id', 'rxn'], 
      [['id', 'r']], 
      ['if', 
        ['where', ['==', 
          ['trailer', ['id', 'len'], ['apply', ['where', ['trailer', ['id', 'r'], ['attribute', ['id', 'reactants']]]]]], 
          ['trailer', ['num', '1']]]], 
        ['where', 
          ['trailer', ['list', 
            ['where', ['trailer', ['id', 'unimolecular'], ['apply', ['where', ['trailer', ['id', 'r']]]]]], 
            ['where', ['trailer', ['id', 'print'], ['apply', ['quote', "'do this!'"]]]]]
        ]], 

        ['where', ['trailer', ['id', 'print'], ['apply', ['quote', "'warn me'"]]]]]]]


    print self.m5
    print np.ts_parser.parse(self.m5)


if __name__ == '__main__':
  unittest.main()

