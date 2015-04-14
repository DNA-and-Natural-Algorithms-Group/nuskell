#
#  test_utils.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/21/10.
#

import unittest
from utils import *
from nose.tools import *
import copy
import input

class MiscTests(unittest.TestCase):
	def testWrap(self):
		assert wrap(0,3) == 0
		assert wrap(2,3) == 2
		assert wrap(3,3) == 0
		assert wrap(4,3) == 1
		assert wrap(-1,3) == 2
		assert wrap(-2,3) == 1

	def testNaturalSort(self):
		assert natural_sort(['c10','b','c2','c1','a']) == ['a','b','c1','c2','c10']
		
	def testParseDotParen(self):
		#						012345678
		assert parse_dot_paren('(((...)))') == [[(0,8),(0,7),(0,6),None,None,None,(0,2),(0,1),(0,0)]]
		
		#						0	1
		#						012 012345
		assert parse_dot_paren('(((+...)))') == [[(1,5),(1,4),(1,3)],[None,None,None,(0,2),(0,1),(0,0)]]
		
		#						0    1
		#						0123 01234
		assert parse_dot_paren('((((+))).)') == [[(1,4),(1,2),(1,1),(1,0)],[(0,3),(0,2),(0,1),None,(0,0)]]

	def testParseConcentration(self):
		assert parse_concentration('0.05 uM') == 5e-8
		assert parse_concentration('17e-4 M') == 17e-4
		assert parse_concentration('2e-5pM') == 2e-17
		assert parse_concentration('1pM') == 1e-12

	def testParseParameters(self):
		assert parse_parameters('[1nt]') == { 'concentration': None }
		assert parse_parameters('[@1uM]') == { 'concentration': 1e-6 }
		assert parse_parameters('1nt@1pM') == { 'concentration': 1e-12 }

	def testFormatSI(self):
		assert format_si(1e-6) == (1, 'u')
		assert format_si(5e-9) == (5, 'n')
		assert format_si(50e-12) == (50, 'p')

class LoopTests(unittest.TestCase):
	def testLoop(self):
		(domains, strands, complexes) = input.from_kernel(["C = 1 2 3() + 4"])
		
		parts = []
		domain_length = 0
		complex = complexes['C']

		# skip the closing helical domain
		for loc in [(0,0), (0,1), (0,2), None, (1,0)]:
			if loc is None:
				parts.append(None)
				continue

			domain = complex.get_domain(loc)
			parts.append( (domain, complex.get_structure(loc), loc) )
			
			# we'll naively assume all domains have the same length for the 
			# sake of this example
			domain_length = len(domain)

		# construct a loop from parts
		loop = Loop(parts)

		assert loop.bases == domain_length*3
		assert loop.stems == 1
		assert loop.is_open == True

class DomainTests(unittest.TestCase):
	def setUp(self):
		self.d1 = Domain('d1', 5)
		self.d2 = Domain('d2', 4, True)
		self.d3 = Domain('d3', 'short')
		self.d4 = Domain('d4', 'long')
		self.d5 = Domain('d5', 6, sequence='ATGCGA')
		self.d6  = Domain('d6', 6, is_complement=True, sequence='ATGCGA')
		self.d6a = Domain('d6', 6, is_complement=False, sequence='ATGCGA')
		self.d6b = Domain('d6',7, sequence='ATGCGA')
		
		
		self.d7 = Domain('d1', 5, False)
		self.d8 = Domain('d2', 4)
		
	def testRepr(self):
		assert self.d1.__repr__() == "Domain(d1)"
		assert self.d2.__repr__() == "Domain(d2*)"
		assert self.d3.__repr__() == "Domain(d3)"
		assert self.d4.__repr__() == "Domain(d4)"
		assert self.d5.__repr__() == "Domain(d5)"
		assert self.d6.__repr__() == "Domain(d6*)"
	
	def testIdentity(self):
		assert self.d1.identity == 'd1'
		assert self.d2.identity == 'd2'
		assert self.d7.identity == 'd1'
		assert self.d6.identity == 'd6'
		def assnId(self):
			self.d1.identity = 'd2'
		assert_raises(AttributeError, assnId, self)
	
	def testName(self):
		assert self.d1.name == 'd1'
		assert self.d2.name == 'd2*'
		def assnName(self):
			self.d1.name = 'd2'
		assert_raises(AttributeError, assnName, self)
	
	def testEq(self):
		assert not (self.d1 == self.d2)
		assert self.d1 == self.d7
		assert not (self.d2 == self.d8)
	
	def testCmp(self):
		assert self.d1 < self.d2 # different name
		assert self.d5 < self.d6 # different name
		assert self.d6 > self.d5 # different name
		assert self.d6 > self.d6a # same name, different complement
		assert self.d6a < self.d6 # same name, different complement
		assert self.d6a < self.d6b # same name, different length
		
	def testStr(self):
		assert self.d1.__str__() == self.d1.name
	
	def testLength(self):
		assert self.d1.length == 5
		assert len(self.d1) == 5
		
		assert self.d2.length == 4
		assert len(self.d2) == 4
		
		assert self.d3.length == SHORT_DOMAIN_LENGTH
		assert self.d4.length == LONG_DOMAIN_LENGTH
		def assnLen(self):
			self.d1.length = 1
		assert_raises(AttributeError, assnLen, self)
	
	def testSequence(self):
		assert self.d1.sequence == None
		assert self.d5.sequence == 'ATGCGA'
		assert self.d6.sequence == 'TCGCAT'
		def assnSeq(self):
			self.d1.sequence = 'ATCG'
		assert_raises(AttributeError, assnSeq, self)
	
	def testIsComplement(self):
		assert not self.d1.is_complement
		assert self.d2.is_complement
		assert not self.d7.is_complement
		def assnIC(self):
			self.d1.is_complement = True
		assert_raises(AttributeError, assnIC, self)
			
# TODO: maybe should replace this with something published...
def setUpSLC(self):
		self.domains = {}
		self.domains['1'] = Domain('1', 'short')
		self.domains['1*'] = Domain('1', 'short', True)
		self.domains['2'] = Domain('2', 'short')
		self.domains['2*'] = Domain('2', 'short', True)
		self.domains['3'] = Domain('3', 'short')
		self.domains['3*'] = Domain('3', 'short', True)
		self.domains['4'] = Domain('4', 'long')
		self.domains['4*'] = Domain('4', 'long', True)
		self.domains['5'] = Domain('5', 'short')
		self.domains['5*'] = Domain('5', 'short', True)
		self.domains['6'] = Domain('6', 'long')
		self.domains['6*'] = Domain('6', 'long', True)
		self.domains['7'] = Domain('7', 'short')
		self.domains['7*'] = Domain('7', 'short', True)
		
		self.strands = {}
		self.strands['PS'] = Strand('PS', [self.domains['3*'], self.domains['2*'], self.domains['1*'], self.domains['5'], self.domains['6']])
		self.strands['OP'] = Strand('OP', [self.domains['1'], self.domains['2'], self.domains['3'], self.domains['4']])
		self.strands['SP'] = Strand('SP', [self.domains['5'], self.domains['6']])
		self.strands['BS'] = Strand('BS', [self.domains['7*'], self.domains['6*'], self.domains['5*'], self.domains['1'], self.domains['2'], self.domains['3']])
		self.strands['Cat'] = Strand('Cat', [self.domains['6'], self.domains['7']])
		
		self.complexes = {}
		self.complexes['C1'] = Complex('C1', [self.strands['PS'], self.strands['OP']], [[(1, 2), (1, 1), (1, 0), None, None], [(0, 2), (0, 1), (0, 0), None]])
		self.complexes['I1'] = Complex('I1', [self.strands['SP'], self.strands['Cat'], self.strands['BS']], [[(2, 2), (2, 1)], [None, (2, 0)], [(1, 1), (0, 1), (0, 0), None, None, None]])
		self.complexes['Cat'] = Complex('Cat', [self.strands['Cat']], [[None, None]])
		self.complexes['I4'] = Complex('I4', [self.strands['BS'], self.strands['OP'], self.strands['PS'], self.strands['Cat']], [[(3,1), (3, 0), (2, 3), None, None, None],[(2, 2), (2, 1), (2, 0), None],[(1, 2), (1, 1), (1, 0), (0, 2), None],[(0, 1), (0, 0)]])
		self.complexes['I3'] = Complex('I3', [self.strands['BS'], self.strands['Cat']], [[(1, 1), (1, 0), None, None, None, None], [(0, 1), (0, 0)]])
		
class StrandTests(unittest.TestCase):
	def setUp(self):
		setUpSLC(self)

	def testEq(self):
		assert not self.strands['PS'] == self.strands['OP']
		PSfake = Strand('PS', [self.domains['3*'], self.domains['2'], self.domains['1*'], self.domains['5'], self.domains['6']])
		assert not self.strands['PS'] == PSfake
		PS = Strand('PS', [self.domains['3*'], self.domains['2*'], self.domains['1*'], self.domains['5'], self.domains['6']])
		assert self.strands['PS'] == PS
	
	def testHash(self):
		assert not hash(self.strands['PS']) == hash(self.strands['OP'])
		PSfake = Strand('PS', [self.domains['3*'], self.domains['2'], self.domains['1*'], self.domains['5'], self.domains['6']])
		assert not hash(self.strands['PS']) == hash(PSfake)
		PS = Strand('PS', [self.domains['3*'], self.domains['2*'], self.domains['1*'], self.domains['5'], self.domains['6']])
		assert hash(self.strands['PS']) == hash(PS)
		
	def testName(self):
		assert self.strands['PS'].name == 'PS'
		def assnName(self):
			self.strands['PS'].name = 'OP'
		assert_raises(AttributeError, assnName, self)
		
	def testDomains(self):
		assert self.strands['PS'].domains == [self.domains['3*'], self.domains['2*'], self.domains['1*'], self.domains['5'], self.domains['6']]
		def assnDoms(self):
			self.strands['PS'].domains = []
		assert_raises(AttributeError, assnDoms, self)
		
	def testLength(self):
		assert self.strands['PS'].length == 5
		assert self.strands['Cat'].length == 2
		def assnLength(self):
			self.strands['PS'].length = []
		assert_raises(AttributeError, assnLength, self)

class ComplexTests(unittest.TestCase):
	def setUp(self):
		setUpSLC(self)
		
		self.complexes['C2'] = Complex( 'C2', [self.strands['OP'], self.strands['PS']], [[None, (1, 1), (1, 0), None], [(0, 2), (0, 1), None, None, None]] )

		# not a feasible structure; just for testing __hash__ing
		self.complexes['C3'] = Complex( 'C3', [self.strands['OP'], self.strands['OP']], [[None, (1, 1), (1, 0), None], [(0, 2), (0, 1), None, None, None]] )
		
	
	def testConstructor(self):
		assert self.complexes['C1'].strands == [self.strands['OP'], self.strands['PS']]
		assert self.complexes['C1'].structure == [[(1, 2), (1, 1), (1, 0), None], [(0, 2), (0, 1), (0, 0), None, None]]

		assert self.complexes['I1'].strands == [self.strands['BS'], self.strands['SP'], self.strands['Cat']]
		assert self.complexes['I1'].structure == [[(2, 1), (1, 1), (1, 0), None, None, None], [(0, 2), (0, 1)], [None, (0, 0)]]
		assert self.complexes['Cat'].strands == [self.strands['Cat']]
		assert self.complexes['Cat'].structure == [[None, None]]
		assert self.complexes['I4'].strands == [self.strands['BS'], self.strands['OP'], self.strands['PS'], self.strands['Cat']]
						
	def testEq(self):
		assert not self.complexes['C1'] == self.complexes['I1']		
		assert self.complexes['C1'] == copy.deepcopy(self.complexes['C1'])
		c = copy.deepcopy(self.complexes['C1']).rotate_strands()		
		assert not self.complexes['C1'] == c
	
	def testHash(self):
		# copy before hash is computed
		c = copy.deepcopy(self.complexes['C1']) #.rotate_strands()	
		assert hash(self.complexes['C1']) == hash(c)

		assert not hash(self.complexes['C1']) == hash(self.complexes['I1'])
		assert hash(self.complexes['C1']) == hash(copy.deepcopy(self.complexes['C1']))
		
		# same strands, different structure
		assert not hash(self.complexes['C1']) == hash(self.complexes['C2'])
	
		# different strands, same structure
		assert not hash(self.complexes['C3']) == hash(self.complexes['C2'])
	
	def testGetDomain(self):
		print self.complexes['C1'].get_domain((0,0))
		assert self.complexes['C1'].get_domain((0,0)) == self.domains['1']
		
		print self.complexes['C1'].get_domain((0,1))
		assert self.complexes['C1'].get_domain((0,1)) == self.domains['2']
		
		print self.complexes['C1'].get_domain((1,0))
		assert self.complexes['C1'].get_domain((1,0)) == self.domains['3*']
		
	def testName(self):
		assert self.complexes['C1'].name == 'C1'
		def assnName(self):
			self.complexes['C1'].name = 'C2'
		assert_raises(AttributeError, assnName, self)
		
	def testStrands(self):
		assert self.complexes['C1'].strands == [self.strands['OP'], self.strands['PS']]
		def assnStrands(self):
			self.complexes['C1'].strands = []
		assert_raises(AttributeError, assnStrands, self)
	
	def testStrandIndex(self):
		assert self.complexes['C1'].strand_index('OP') == 0
		assert self.complexes['C1'].strand_index('PS') == 1
		assert self.complexes['C1'].strand_index('nonexistent') == -1	
	
	def testStructure(self):
		assert self.complexes['C1'].structure == [[(1, 2), (1, 1), (1, 0), None], [(0, 2), (0, 1), (0, 0), None, None]]
		def assnStructure(self):
			self.complexes['C1'].structure = []
		assert_raises(AttributeError, assnStructure, self)

	def testAvailableDomains(self):
		C1doms = [(self.domains['4'], 0, 3), (self.domains['5'], 1, 3), (self.domains['6'], 1, 4)]
		C1doms.sort(key=lambda dom: dom[0].name)
		assert self.complexes['C1'].available_domains == C1doms
		
		I1doms = [(self.domains['1'], 0, 3), (self.domains['2'], 0, 4), (self.domains['3'], 0, 5), (self.domains['6'], 2, 0)]		
		I1doms.sort(key=lambda dom: dom[0].name)		
		assert self.complexes['I1'].available_domains == I1doms
		
		t1 = Strand('t1', [self.domains['1'], self.domains['2'], self.domains['1*']])
		ct1 = Complex('ct1', [t1], [[(0, 2), None, (0, 0)]])
		assert ct1.available_domains == []
		
		t2 = Strand('t2', [self.domains['1'], self.domains['2'], self.domains['3']])
		t3 = Strand('t3', [self.domains['3*'], self.domains['2'], self.domains['1*']])
		ct2 = Complex('ct2', [t2, t3], [[(1, 2), None, (1, 0)], [(0, 2), None, (0, 0)]])
		assert ct2.available_domains == []

	
	def testAvailableDomains2(self):
		from input import input_enum

		# Example from 3-arm junction

		
		three_arm_nodal_enum = input_enum('test_files/examples/3-arm-junction.enum')
		enum = three_arm_nodal_enum
		
		domains = {}
		strands = {}
		complexes = {}
		
		for domain in enum.domains:
			domains[domain.name] = domain
		
		for strand in enum.strands:
			strands[strand.name] = strand
		
		for complex in enum.initial_complexes:
			complexes[complex.name] = complex
		
		C24 = Complex('24',[strands['n2'], strands['n3']], [[None, None, (1, 1), (1, 0)], [(0, 3), (0, 2), None, None]])
		assert sorted(C24.available_domains) == sorted([(domains['1'], 0, 0), (domains['1*'], 1, 2), (domains['2'], 0, 1), (domains['3*'], 1, 3)])

		C24_ = Complex('C24_',[strands['n2'], strands['n3']],[[(1, 2), None, (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 0), None]])
		assert C24_.available_domains == [(domains['3*'], 1, 3)]
		
		# Original test case
		#rxns = bind11(C24)
		#assert rxns[0].products[0] == C24_
		#assert rxns[0].products[0].available_domains == [(domains['3*'], 1, 3)]
														
	def testRotateStrands(self):	
		c2 = self.complexes['I1'].rotate_strands()
		c3 = c2.rotate_strands()
		assert c3.strands == [self.strands['Cat'], self.strands['BS'], self.strands['SP']]
		assert c3.structure == [[None, (1, 0)], [(0, 1), (2, 1), (2, 0), None, None, None], [(1, 2), (1, 1)]]
		c4 = self.complexes['Cat'].rotate_strands()
		assert c4.strands == [self.strands['Cat']]
		assert c4.structure == [[None, None]]
	
	def testRotateStrands2(self):
		from input import input_enum
		self.biggate_enum = input_enum('test_files/examples/sarma2010/biggate.in');
		(domains,strands,complexes) = index_parts(self.biggate_enum)
	
		c16 = Complex('16',[strands['a1'],strands['a1'],strands['b1']],parse_dot_paren('(..+(((+))))'))
		c38 = Complex('38',[strands['a1'],strands['b1'],strands['a1']],parse_dot_paren('(((+)))(+)..'))
		assert c16 == c38
	
	def testRotateLocation(self):
		i1 = self.complexes['I1']
		c3 = Complex('C3', [self.strands['SP'], self.strands['Cat'], self.strands['BS']], [[(2, 2), (2, 1)], [None, (2, 0)], [(1, 1), (0, 1), (0, 0), None, None, None]])
		assert c3.strands == [self.strands['BS'], self.strands['SP'], self.strands['Cat']]
		assert c3.structure == [[(2, 1), (1, 1), (1, 0), None, None, None], [(0, 2), (0, 1)], [None, (0, 0)]]
		assert c3 == self.complexes['I1']
		assert c3._rotations == -2


		assert c3.rotate_location((2,1)) == (0,1)
		assert c3.rotate_location((0,0)) == (1,0)
		assert c3.rotate_location((1,0)) == (2,0)

		# assert c3.strands == [self.strands['Cat'], self.strands['BS'], self.strands['SP']]
		# assert c3.structure == 
		# c4 = self.complexes['Cat'].rotate_strands()
		# assert c4.strands == [self.strands['Cat']]
		# assert c4.structure == [[None, None]]

	def testDotParenString(self):
		str = self.complexes['I1'].dot_paren_string()
		assert str == "(((...+))+.)"
		
		str = self.complexes['I4'].dot_paren_string()
		assert str == "(((...+(((.+)))).+))"
	
	def testKernelString(self):
		assert self.complexes['I4'].kernel_string() == "7*( 6*( 5*( 1 2 3 + 1( 2( 3( 4 + ) ) ) ) 6 + ) )"
		
	def testCheckStructure(self):
		#                 0                 1                 2                 3                  4
		s1 = Strand('S1',[self.domains['1'],self.domains['2'],self.domains['3'],self.domains['1*'],self.domains['2*']])
		c1 = Complex('C1',[s1],[(0,4),(0,3),None,(0,1),(0,0)])
		def checkStruct1():
			c1.check_structure()
		assert_raises(Exception,checkStruct1)
	
		c2 = Complex('C2',[s1],[(0,4),None,None,None,None])
		def checkStruct2():
			c2.check_structure()
		assert_raises(Exception,checkStruct2)
	
		assert self.complexes['C1'].check_structure()
		assert self.complexes['C2'].check_structure()
		assert self.complexes['Cat'].check_structure()		
	
class RestingStateTests(unittest.TestCase):
	def setUp(self):
		setUpSLC(self)
		self.rs = RestingState('RS1', [self.complexes['C1'], self.complexes['Cat'], self.complexes['I1']])
		self.rs1 = RestingState('RS2', [self.complexes['C1'], \
									Complex('123',self.complexes['Cat'].strands,self.complexes['Cat'].structure), \
									Complex('456',self.complexes['I1'].strands,self.complexes['I1'].structure)])
	
	def testConstructor(self):
		assert self.rs._name == 'RS1'
				
		comp = sorted([self.complexes['C1'], self.complexes['Cat'], self.complexes['I1']])
		
		assert self.rs._complexes == comp
		
		
	def testName(self):
		assert self.rs.name == 'RS1'
		def assnName(self):
			self.rs.name = 'RS2'
			
		assert_raises(AttributeError, assnName, self)
		
	def testStr(self):
		print str(self.rs)
		print str(self.rs1)
		
		assert str(self.rs) == 'I1'
		assert str(self.rs1) == 'C1'
		
	def testRepr(self):
		print repr(self.rs)
		assert repr(self.rs) == "RestingState(\"RS1\", %s)" % str(self.rs.complexes)
		
	def testComplexes(self):	
		self.rs.complexes[0] = None
		
		comp = sorted([self.complexes['C1'], self.complexes['Cat'], self.complexes['I1']])
		
		assert self.rs.complexes == comp
		
		def assnComplex(self):
			self.rs.complexes = []
		
		assert_raises(AttributeError, assnComplex, self)