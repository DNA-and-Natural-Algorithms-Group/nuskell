#
#  test_input.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 1/8/11.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#

import unittest
from utils import *
from input import *
from nose.tools import *
import copy

class InputStandardTests(unittest.TestCase):
	def testStandard_SLC(self):
		enumerator = input_enum('test_files/test_input_standard_SLC.in')
		
		d1 = Domain('1', 'short')
		d1c = Domain('1', 'short', True)
		d2 = Domain('2', 'short')
		d2c = Domain('2', 'short', True)
		d3 = Domain('3', 'short')
		d3c = Domain('3', 'short', True)
		d4 = Domain('4', 'long')
		d4c = Domain('4', 'long', True)
		d5 = Domain('5', 'short')
		d5c = Domain('5', 'short', True)
		d6 = Domain('6', 'long')
		d6c = Domain('6', 'long', True)
		d7 = Domain('7', 'short')
		d7c = Domain('7', 'short', True)
		
		enum_doms = enumerator.domains[:]
		doms = [d1, d1c, d2, d2c, d3, d3c, d4, d4c, d5, d5c, d6, d6c, d7, d7c]
		
		enum_doms.sort()
		doms.sort()
		
		assert (enum_doms == doms)
		
		PS = Strand('PS', [d3c, d2c, d1c, d5, d6])
		OP = Strand('OP', [d1, d2, d3, d4])
		SP = Strand('SP', [d5, d6])
		BS = Strand('BS', [d7c, d6c, d5c, d1, d2, d3])
		Cat = Strand('Cat', [d6, d7])
		
		enum_strands = enumerator.strands[:]
		strands = [PS, OP, SP, BS, Cat]
		
		enum_strands.sort()
		strands.sort()
		
		assert (enum_strands == strands)
		
		C1 = Complex('C1', [PS, OP], [[(1, 2), (1, 1), (1, 0), None, None], [(0, 2), (0, 1), (0, 0), None]])
		C2 = Complex('C2', [SP, BS], [[(1, 2), (1, 1)], [None, (0, 1), (0, 0), None, None, None]])

		CatC = Complex('Cat', [Cat], [[None, None]])
		OPC = Complex('OP', [OP], [[None, None, None, None]])
		SPC = Complex('SP', [SP], [[None, None]])
		
		I1 = Complex('I1', [SP, Cat, BS], [[(2, 2), (2, 1)], [None, (2, 0)], [(1, 1), (0, 1), (0, 0), None, None, None]])
		I2 = Complex('I2', [SP, Cat, BS], [[(2, 2), None], [(2, 1), (2, 0)], [(1, 1), (1, 0), (0, 0), None, None, None]])
		I3 = Complex('I3', [BS, Cat], [[(1, 1), (1, 0), None, None, None, None], [(0, 1), (0, 0)]])
		I4 = Complex('I4', [BS, OP, PS, Cat], [[(3,1), (3, 0), (2, 3), None, None, None],[(2, 2), (2, 1), (2, 0), None],[(1, 2), (1, 1), (1, 0), (0, 2), None],[(0, 1), (0, 0)]])
		I5 = Complex('I5', [BS, PS, Cat], [[(2, 1), (2, 0), (1, 3), (1, 2), (1, 1), (1, 0)], [(0, 5), (0, 4), (0, 3), (0, 2), None], [(0, 1), (0, 0)]])
		I6 = Complex('I6', [BS, PS, Cat], [[(2, 1), (1, 4), (1, 3), (1, 2), (1, 1), (1, 0)], [(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)], [None, (0, 0)]])
		W = Complex('W', [BS, PS], [[None, (1, 4), (1, 3), (1, 2), (1, 1), (1, 0)], [(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)]])
		
		enum_complexes = enumerator.initial_complexes[:]
		complexes = [C1, C2, CatC, OPC, SPC, I1, I2, I3, I4, I5, I6, W]
		
		
		enum_complexes.sort()
		complexes.sort()
		
		assert (enum_complexes == complexes)

	def testStandard_3arm_junction(self):
		enumerator = input_enum('test_files/test_input_standard_3arm_junction.in')
		
		da = Domain('a', 6)		
		db = Domain('b', 6)
		dc = Domain('c', 6)
		dx = Domain('x', 6)
		dy = Domain('y', 6)
		dz = Domain('z', 6)
		dac = Domain('a', 6, True)
		dbc = Domain('b', 6, True)
		dcc = Domain('c', 6, True)
		dxc = Domain('x', 6, True)
		dyc = Domain('y', 6, True)
		dzc = Domain('z', 6, True)
		
		enum_doms = enumerator.domains[:]
		doms = [da, dac, db, dbc, dc, dcc, dx, dxc, dy, dyc, dz, dzc]
				
		enum_doms.sort()
		doms.sort()
		
		assert (enum_doms == doms)
		
		I = Strand('I', [dyc, dbc, dxc, dac])
		A = Strand('A', [da, dx, db, dy, dzc, dcc, dyc, dbc, dxc])
		B = Strand('B', [db, dy, dc, dz, dxc, dac, dzc, dcc, dyc])
		C = Strand('C', [dc, dz, da, dx, dyc, dbc, dxc, dac, dzc])
		
		enum_strands = enumerator.strands[:]
		strands = [I, A, B, C]
		
		enum_strands.sort()
		strands.sort()
		
		assert (enum_strands == strands)
		
		IC = Complex('I', [I], [[None, None, None, None]])
		AC = Complex('A', [A], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		BC = Complex('B', [B], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		CC = Complex('C', [C], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		
		IA = Complex('IA', [I, A], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), None, None, None, None, None]])
		IAB = Complex('IAB', [I, A, B], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), None, None, None, None, None]])
		IABC = Complex('IABC', [I, A, B, C], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (3, 3), (3, 2), (3, 1), (3, 0), None], [(2, 7), (2, 6), (2, 5), (2, 4), None, None, None, None, None]])
		ABC = Complex('ABC', [A, B, C], [[(2, 7), (2, 6), (2, 5), (2, 4), (1, 3), (1, 2), (1, 1), (1, 0), None], [(0, 7), (0, 6), (0, 5), (0, 4), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (0, 3), (0, 2), (0, 1), (0, 0), None]])
		
		enum_complexes = enumerator.initial_complexes[:]
		complexes = [IC, AC, BC, CC, IA, IAB, IABC, ABC]
		
		
		enum_complexes.sort()
		complexes.sort()
		
		
		assert (enum_complexes == complexes)
	
	def testPil(self):
		enumerator = input_pil('test_files/test_input_pil.pil')
		
		# domains
		da = Domain('a', 6, sequence="CTACTC")		
		db = Domain('b-seq', 8, sequence="ACATCGAN")
		dz = Domain('z', 6, sequence="TTTCCA")
		dac = Domain('a', 6, True, sequence="CTACTC")
		dbc = Domain('b-seq', 8, True, sequence="ACATCGAN")
		dzc = Domain('z', 6, True, sequence="TTTCCA")
		
		dq = Domain('q', 20, sequence="CTACTCACATCGANTTTCCA")
		dqc = Domain('q', 20, True, sequence="CTACTCACATCGANTTTCCA")

		enum_doms = enumerator.domains[:]
		doms = [da, dac, db, dbc, dz, dzc, dq, dqc]
				
		enum_doms.sort()
		doms.sort()
		
		print
		print "Domains in enumerator:"
		print enum_doms
		print "Expected:"
		print doms
		print [cmp(enum_doms[i],doms[i]) for i in xrange(len(doms)) ]
		assert (tuple(enum_doms) == tuple(doms))

		# strands
		A = Strand('A', [da, db, dz, dq])
		B = Strand('B', [dqc, dzc, dbc])
		
		enum_strands = enumerator.strands[:]
		strands = [A, B]
		
		enum_strands.sort()
		strands.sort()
		
		print
		print "Strands in enumerator:"
		print enum_strands
		print "Expected:"
		print strands
		print [enum_strands[i] == strands[i] for i in xrange(len(strands)) ]
		assert (tuple(enum_strands) == tuple(strands))
		
		# complexes
		# structure AB = A B : .(((+)))
		AB = Complex('AB', [A,B], [[None, (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1)]])

		enum_complexes = enumerator.initial_complexes[:]
		complexes = [AB]
		
		enum_complexes.sort()
		complexes.sort()
		
		print
		print "Complexes in enumerator:"
		print enum_complexes
		print "Expected:"
		print complexes
		print [enum_complexes[i] == complexes[i] for i in xrange(len(complexes)) ]
		assert (tuple(enum_complexes) == tuple(complexes))

	def testPil_venkataraman2007(self):
		enumerator = input_pil('test_files/examples/venkataraman2007/venkataraman2007.pil')
		

	def testPil_3arm_junction(self):
		enumerator = input_pil('test_files/test_input_pil_3arm_junction.pil')
		
		da = Domain('a', 6, sequence="CTACTC")		
		db = Domain('b', 6, sequence="TCCTCA")
		dc = Domain('c', 6, sequence="TTTCCA")
		dx = Domain('x', 6, sequence="CTACTC")
		dy = Domain('y', 6, sequence="TCCTCA")
		dz = Domain('z', 6, sequence="TTTCCA")
		dac = Domain('a', 6, True, sequence="CTACTC")
		dbc = Domain('b', 6, True, sequence="TCCTCA")
		dcc = Domain('c', 6, True, sequence="TTTCCA")
		dxc = Domain('x', 6, True, sequence="CTACTC")
		dyc = Domain('y', 6, True, sequence="TCCTCA")
		dzc = Domain('z', 6, True, sequence="TTTCCA")
		
		enum_doms = enumerator.domains[:]
		doms = [da, dac, db, dbc, dc, dcc, dx, dxc, dy, dyc, dz, dzc]
				
		enum_doms.sort()
		doms.sort()
		
		assert (enum_doms == doms)
		
		I = Strand('I', [dyc, dbc, dxc, dac])
		A = Strand('A', [da, dx, db, dy, dzc, dcc, dyc, dbc, dxc])
		B = Strand('B', [db, dy, dc, dz, dxc, dac, dzc, dcc, dyc])
		C = Strand('C', [dc, dz, da, dx, dyc, dbc, dxc, dac, dzc])
		
		enum_strands = enumerator.strands[:]
		strands = [I, A, B, C]
		
		enum_strands.sort()
		strands.sort()
		
		assert (enum_strands == strands)
		
		IC = Complex('I', [I], [[None, None, None, None]])
		AC = Complex('A', [A], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		BC = Complex('B', [B], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		CC = Complex('C', [C], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		
		IA = Complex('IA', [I, A], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), None, None, None, None, None]])
		IAB = Complex('IAB', [I, A, B], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), None, None, None, None, None]])
		IABC = Complex('IABC', [I, A, B, C], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (3, 3), (3, 2), (3, 1), (3, 0), None], [(2, 7), (2, 6), (2, 5), (2, 4), None, None, None, None, None]])
		ABC = Complex('ABC', [A, B, C], [[(2, 7), (2, 6), (2, 5), (2, 4), (1, 3), (1, 2), (1, 1), (1, 0), None], [(0, 7), (0, 6), (0, 5), (0, 4), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (0, 3), (0, 2), (0, 1), (0, 0), None]])
		
		enum_complexes = enumerator.initial_complexes[:]
		complexes = [IC, AC, BC, CC, IA, IAB, IABC, ABC]
		
		
		enum_complexes.sort()
		complexes.sort()
		
		
		assert (enum_complexes == complexes)
	
	
	def testErrors(self):
		
		def testDuplicateDomain():
			enum = input_enum('test_files/test_input_errors/test_input_duplicate_domain.in')
		assert_raises(Exception, testDuplicateDomain)
		
		def testDuplicateStrand():
			enum = input_enum('test_files/test_input_errors/test_input_duplicate_strand.in')	
		assert_raises(Exception, testDuplicateStrand)
		
		def testDuplicateComplex():
			enum = input_enum('test_files/test_input_errors/test_input_duplicate_complex.in')	
		assert_raises(Exception, testDuplicateComplex)

		def testMissingDomain():
			enum = input_enum('test_files/test_input_errors/test_input_missing_domain.in')	
		assert_raises(Exception, testMissingDomain)

		def testMissingStrand():
			enum = input_enum('test_files/test_input_errors/test_input_missing_strand.in')	
		assert_raises(Exception, testMissingStrand)
		
		def testComplexError():
			enum = input_enum('test_files/test_input_errors/test_input_size_mismatch.in')
		assert_raises(Exception, testComplexError)
		
		enum = input_enum('test_files/test_input_errors/test_input_warnings.in')	
		assert (enum != None)
			
	
class InputKernel(unittest.TestCase):

	def test_parse_kernel(self):
		def listify(x):
			return [i if isinstance(i, basestring) else listify(i) for i in x] 

		assert listify(parse_kernel("a")) == ['a']
		assert listify(parse_kernel("a*")) == ['a*']
		assert listify(parse_kernel("a^*")) == ['a^*']
		assert listify(parse_kernel("1^*")) == ['1^*']
		assert listify(parse_kernel("foo^*")) == ['foo^*']
		assert listify(parse_kernel("a()")) == ['a', []]
		assert listify(parse_kernel("a(+)")) == ['a', ['+']]
		assert listify(parse_kernel("a+b")) == ['a', '+', 'b']
		assert listify(parse_kernel("a(b +)")) == ['a', ['b', '+']]
		assert listify(parse_kernel("a b(c) d")) == ['a', 'b', ['c'], 'd']
		assert listify(parse_kernel("a(b c) d^ e(f)")) == ['a', ['b', 'c'], 'd^', 'e', ['f']]
		assert listify(parse_kernel("a b c*(d + e) f a")) == ['a', 'b', 'c*', ['d', '+', 'e'], 'f', 'a']

		# parse_kernel("^foo")
		assert_raises(Exception, lambda: parse_kernel("^foo"))

		# parse_kernel("a(")
		assert_raises(Exception, lambda: parse_kernel("a("))
		assert_raises(Exception, lambda: parse_kernel("b)"))
		assert_raises(Exception, lambda: parse_kernel("a())"))

	def test_kernel_helpers(self):
		assert parse_identifier("a") == ("a", 1, 'long')
		assert parse_identifier("b^") == ("b", 1, 'short')
		assert parse_identifier("c*") == ("c", -1, 'long')
		assert parse_identifier("d^*") == ("d", -1, 'short')

		domains = {}
		a = auto_domain("a", 1, domains)
		ap = auto_domain("a*", 1, domains)
		assert a.can_pair(ap)
		assert auto_domain("a*", 1, domains) == ap
		assert auto_domain("a", -1, domains) == ap

		bp = auto_domain("b*", 1, domains)
		b = auto_domain("b", 1, domains)
		assert auto_domain("b*", -1, domains) == b
		assert auto_domain("b", 1, domains) == b

		t = auto_domain("t", 1, domains)
		# assert_raises(SystemExit, lambda: auto_domain("t^", 1, domains))

	def test_from_kernel(self):
		(kdomains, kstrands, kcomplexes) = from_kernel(["FancyComplexName = a( b(c + d) e) f"])

		domains = { 
			'a' : Domain('a', 12, is_complement=False, sequence='None'),
			'a*' : Domain('a', 12, is_complement=True, sequence='None'),
			'b' : Domain('b', 12, is_complement=False, sequence='None'),
			'b*' : Domain('b', 12, is_complement=True, sequence='None'),
			'c' : Domain('c', 12, is_complement=False, sequence='None'),
			'c*' : Domain('c', 12, is_complement=True, sequence='None'),
			'd' : Domain('d', 12, is_complement=False, sequence='None'),
			'd*' : Domain('d', 12, is_complement=True, sequence='None'),
			'e' : Domain('e', 12, is_complement=False, sequence='None'),
			'e*' : Domain('e', 12, is_complement=True, sequence='None'),
			'f' : Domain('f', 12, is_complement=False, sequence='None'),
			'f*' : Domain('f', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(kdomains.values())

		# Strands 
		strands = { 
			'a_b_c' : Strand('a_b_c', [domains['a'], domains['b'], domains['c']]),
			'd_b*_e_a*_f' : Strand('d_b*_e_a*_f', [domains['d'], domains['b*'], domains['e'], domains['a*'], domains['f']])
		}
		assert set(strands.values()) == set(kstrands.values())

		# Complexes 
		complexes = { 
			'FancyComplexName' : Complex('FancyComplexName', [strands['a_b_c'], strands['d_b*_e_a*_f']], [[(1, 3), (1, 1), None], [None, (0, 1), None, (0, 0), None]])
		}
		assert complexes == kcomplexes

	def test_kernel_1(self):
		enumerator = input_pil('test_files/test_input_kernel_1.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'a' : Domain('a', 12, is_complement=False, sequence='None'),
			'a*' : Domain('a', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'a' : Strand('a', [domains['a']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'1' : Complex('1', [strands['a']], [[None]]),
			'2' : Complex('2', [strands['a']], [[None]])
		}
		assert set(complexes.values()) == set(enumerator.complexes)
		

	def test_kernel_2(self):
		enumerator = input_pil('test_files/test_input_kernel_2.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'a' : Domain('a', 12, is_complement=False, sequence='None'),
			'a*' : Domain('a', 12, is_complement=True, sequence='None'),
			'b' : Domain('b', 12, is_complement=False, sequence='None'),
			'b*' : Domain('b', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'a_b' : Strand('a_b', [domains['a'], domains['b']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'1' : Complex('1', [strands['a_b']], [[None, None]])
		}
		assert set(complexes.values()) == set(enumerator.complexes)


	def test_kernel_3(self):
		enumerator = input_pil('test_files/test_input_kernel_3.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'a' : Domain('a', 12, is_complement=False, sequence='None'),
			'a*' : Domain('a', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'a' : Strand('a', [domains['a']]),
			'a*' : Strand('a*', [domains['a*']]),
			'a_a*' : Strand('a_a*', [domains['a'], domains['a*']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'1' : Complex('1', [strands['a_a*']], [[(0, 1), (0, 0)]]),
			'2' : Complex('2', [strands['a'], strands['a*']], [[(1, 0)], [(0, 0)]]),
			'3' : Complex('3', [strands['a'], strands['a*']], [[(1, 0)], [(0, 0)]])
		}
		assert set(complexes.values()) == set(enumerator.complexes)


	def test_kernel_4(self):
		enumerator = input_pil('test_files/test_input_kernel_4.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'a' : Domain('a', 12, is_complement=False, sequence='None'),
			'a*' : Domain('a', 12, is_complement=True, sequence='None'),
			'b' : Domain('b', 12, is_complement=False, sequence='None'),
			'b*' : Domain('b', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'a' : Strand('a', [domains['a']]),
			'a*_b_a' : Strand('a*_b_a', [domains['a*'], domains['b'], domains['a']]),
			'a_b' : Strand('a_b', [domains['a'], domains['b']]),
			'a_b_a*' : Strand('a_b_a*', [domains['a'], domains['b'], domains['a*']]),
			'b' : Strand('b', [domains['b']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'1' : Complex('1', [strands['a_b']], [[None, None]]),
			'2' : Complex('2', [strands['a'], strands['b']], [[None], [None]]),
			'3' : Complex('3', [strands['a_b_a*']], [[(0, 2), None, (0, 0)]]),
			'4' : Complex('4', [strands['a*_b_a']], [[(0, 2), None, (0, 0)]])
		}
		assert set(complexes.values()) == set(enumerator.complexes)




	def test_kernel_5(self):
		enumerator = input_pil('test_files/test_input_kernel_5.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'a' : Domain('a', 12, is_complement=False, sequence='None'),
			'a*' : Domain('a', 12, is_complement=True, sequence='None'),
			'b' : Domain('b', 12, is_complement=False, sequence='None'),
			'b*' : Domain('b', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'a*_b_a' : Strand('a*_b_a', [domains['a*'], domains['b'], domains['a']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'1' : Complex('1', [strands['a*_b_a']], [[(0, 2), None, (0, 0)]])
		}
		assert set(complexes.values()) == set(enumerator.complexes)



	def test_kernel_6(self):
		enumerator = input_pil('test_files/test_input_kernel_6.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'a' : Domain('a', 12, is_complement=False, sequence='None'),
			'a*' : Domain('a', 12, is_complement=True, sequence='None'),
			'b' : Domain('b', 12, is_complement=False, sequence='None'),
			'b*' : Domain('b', 12, is_complement=True, sequence='None'),
			'c' : Domain('c', 12, is_complement=False, sequence='None'),
			'c*' : Domain('c', 12, is_complement=True, sequence='None'),
			'd' : Domain('d', 12, is_complement=False, sequence='None'),
			'd*' : Domain('d', 12, is_complement=True, sequence='None'),
			'e' : Domain('e', 12, is_complement=False, sequence='None'),
			'e*' : Domain('e', 12, is_complement=True, sequence='None'),
			'f' : Domain('f', 12, is_complement=False, sequence='None'),
			'f*' : Domain('f', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'a_b_c*_d' : Strand('a_b_c*_d', [domains['a'], domains['b'], domains['c*'], domains['d']]),
			'e_c_f_a' : Strand('e_c_f_a', [domains['e'], domains['c'], domains['f'], domains['a']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'1' : Complex('1', [strands['a_b_c*_d'], strands['e_c_f_a']], [[None, None, (1, 1), None], [None, (0, 2), None, None]])
		}
		assert set(complexes.values()) == set(enumerator.complexes)


	def test_kernel_7(self):
		enumerator = input_pil('test_files/test_input_kernel_7.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'a' : Domain('a', 15, is_complement=False, sequence='N'),
			'a*' : Domain('a', 15, is_complement=True, sequence='N'),
			'b' : Domain('b', 12, is_complement=False, sequence='None'),
			'b*' : Domain('b', 12, is_complement=True, sequence='None'),
			'c' : Domain('c', 12, is_complement=False, sequence='None'),
			'c*' : Domain('c', 12, is_complement=True, sequence='None'),
			'd' : Domain('d', 12, is_complement=False, sequence='None'),
			'd*' : Domain('d', 12, is_complement=True, sequence='None'),
			'e' : Domain('e', 12, is_complement=False, sequence='None'),
			'e*' : Domain('e', 12, is_complement=True, sequence='None'),
			'f' : Domain('f', 12, is_complement=False, sequence='None'),
			'f*' : Domain('f', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'a_b_c*_d' : Strand('a_b_c*_d', [domains['a'], domains['b'], domains['c*'], domains['d']]),
			'e_c_f_a' : Strand('e_c_f_a', [domains['e'], domains['c'], domains['f'], domains['a']]),
			's' : Strand('s', [domains['a'], domains['b'], domains['c'], domains['d'], domains['e'], domains['f']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'1' : Complex('1', [strands['a_b_c*_d'], strands['e_c_f_a']], [[None, None, (1, 1), None], [None, (0, 2), None, None]])
		}
		assert set(complexes.values()) == set(enumerator.complexes)


	def test_kernel_8(self):
		enumerator = input_pil('test_files/test_input_kernel_8.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'2' : Domain('2', 12, is_complement=False, sequence='None'),
			'2*' : Domain('2', 12, is_complement=True, sequence='None'),
			'3' : Domain('3', 12, is_complement=False, sequence='None'),
			'3*' : Domain('3', 12, is_complement=True, sequence='None'),
			'a' : Domain('a', 12, is_complement=False, sequence='None'),
			'a*' : Domain('a', 12, is_complement=True, sequence='None'),
			'b' : Domain('b', 12, is_complement=False, sequence='None'),
			'b*' : Domain('b', 12, is_complement=True, sequence='None'),
			'c' : Domain('c', 12, is_complement=False, sequence='None'),
			'c*' : Domain('c', 12, is_complement=True, sequence='None'),
			't' : Domain('t', 6, is_complement=False, sequence='None'),
			't*' : Domain('t', 6, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'2_3' : Strand('2_3', [domains['2'], domains['3']]),
			'3*_a*' : Strand('3*_a*', [domains['3*'], domains['a*']]),
			'a*_b*_c_b_a_2*_t*' : Strand('a*_b*_c_b_a_2*_t*', [domains['a*'], domains['b*'], domains['c'], domains['b'], domains['a'], domains['2*'], domains['t*']]),
			't_2_3' : Strand('t_2_3', [domains['t'], domains['2'], domains['3']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'1' : Complex('1', [strands['t_2_3']], [[None, None, None]]),
			'2' : Complex('2', [strands['2_3'], strands['3*_a*'], strands['a*_b*_c_b_a_2*_t*']], [[(2, 5), (1, 0)], [(0, 1), None], [(2, 4), (2, 3), None, (2, 1), (2, 0), (0, 0), None]])
		}
		assert set(complexes.values()) == set(enumerator.complexes)

	def test_kernel_9(self):
		enumerator = input_pil('test_files/test_input_kernel_9.pil')
		enumerator.dry_run()

		# Domains 
		domains = { 
			'2' : Domain('2', 12, is_complement=False, sequence='None'),
			'2*' : Domain('2', 12, is_complement=True, sequence='None'),
			'3' : Domain('3', 12, is_complement=False, sequence='None'),
			'3*' : Domain('3', 12, is_complement=True, sequence='None')
		}
		assert set(domains.values()) == set(enumerator.domains)

		# Strands 
		strands = { 
			'2_3' : Strand('2_3', [domains['2'], domains['3']])
		}
		assert set(strands.values()) == set(enumerator.strands)

		# Complexes 
		complexes = { 
			'2' : Complex('2', [strands['2_3']], [[None, None, None]]),
			'3' : Complex('2', [strands['2_3']], [[None, None, None]]),
			'4' : Complex('2', [strands['2_3']], [[None, None, None]]),
			'5' : Complex('2', [strands['2_3']], [[None, None, None]])
		}
		complexes['2'].concentration = 5e-7
		complexes['3'].concentration = 5e-7
		complexes['4'].concentration = 5e-7
		complexes['5'].concentration = 7e-13

		out_domains, out_strands, out_complexes = index_parts(enumerator)
		assert complexes['2'].concentration == out_complexes['2'].concentration
		assert complexes['3'].concentration == out_complexes['3'].concentration
		assert complexes['4'].concentration == out_complexes['4'].concentration
		assert complexes['5'].concentration == out_complexes['5'].concentration
