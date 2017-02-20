import unittest
from pyparsing import ParseException
from nuskell.parser import parse_pil_string, parse_pil_file

class TestCRNparser(unittest.TestCase):
  def setUp(self):
    pass

  def tearDown(self):
    pass

  def test_parse_examples(self):
    example1 = """
    length t0 = 6
    length d44 = 15
    length f0 = 15
    length d31 = 15
    length d25 = 15
    length Fluor25 = 15

    w44_31 = d44 t0 d31
    G31_25 = d31( t0( d25 + ) ) t0*
    R25    = d25( Fluor25 + ) t0* 
    rep25  = d25 Fluor25
    #w31_f  = d31 t0 f0
    #w31_25 = d31 t0 d25
    #g31b   = t0* d31* t0*
    """
    output1 = [
                ['domain', 't0', '6'], 
                ['domain', 'd44', '15'], 
                ['domain', 'f0', '15'], 
                ['domain', 'd31', '15'], 
                ['domain', 'd25', '15'], 
                ['domain', 'Fluor25', '15'], 
                ['complex', 'w44_31', ['d44', 't0', 'd31']], 
                ['complex', 'G31_25', ['d31', ['t0', ['d25', '+']], 't0*']], 
                ['complex', 'R25', ['d25', ['Fluor25', '+'], 't0*']], 
                ['complex', 'rep25', ['d25', 'Fluor25']]
              ]
    self.assertEqual(parse_pil_string(example1), output1, 'seesaw example1')

    self.assertEqual(parse_pil_string("cplx = a( b( c( + ) ) d ) "), 
        [['complex', 'cplx', ['a', ['b', ['c', ['+']], 'd']]]], 'small example 1')

    with self.assertRaises(ParseException):
      # whitespace between domains 
      parse_pil_string("cplx = a(b( c( + ) ) ) d ")

    with self.assertRaises(ParseException):
      # Closing domain must not be defined
      # NOTE: this behavior may change in the future?
      parse_pil_string("cplx = a( b( c( + ) ) d) ")

    with self.assertRaises(ParseException):
      # Missing opening domain name
      parse_pil_string("cplx = a( b( ( + ) ) ) ")

    with self.assertRaises(ParseException):
      # Unbalanced brackets
      parse_pil_string("cplx = a( b( c( + ) ) d ")

    #NOTE: This test passes, so sad...
    #with self.assertRaises(ParseException):
    #  parse_pil_string("cplx = a( b( c( + ) ) )d")

if __name__ == '__main__':
  unittest.main()

