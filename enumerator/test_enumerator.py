#
#  test_enumerator.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 2/6/11.
#

import unittest
from utils import *
from reactions import ReactionPathway
from input import input_enum
from enumerator import *
from nose.tools import *
import copy

class EnumeratorTests(unittest.TestCase):
	def setUp(self):
		self.SLC_enumerator = input_enum('test_files/test_input_standard_SLC.in')
		self.domains = {}
		self.strands = {}
		self.complexes = {}		
		
#		for domain in self.SLC_enumerator.domains:
#			self.domains[domain.name] = domain
#		
#		for strand in self.SLC_enumerator.strands:
#			self.strands[strand.name] = strand
#		
#		for complex in self.SLC_enumerator.initial_complexes:
#			self.complexes[complex.name] = complex
		
		(self.domains,self.strands,self.complexes) = index_parts(self.SLC_enumerator)
		
		
	
		self.three_arm_enumerator = input_enum('test_files/test_input_standard_3arm_junction.in')
		self.domains2 = {}
		self.strands2 = {}
		self.complexes2 = {}		
		
		(self.domains2,self.strands2,self.complexes2) = index_parts(self.three_arm_enumerator)
		
#		for domain in self.three_arm_enumerator.domains:
#			self.domains2[domain.name] = domain
#		
#		for strand in self.three_arm_enumerator.strands:
#			self.strands2[strand.name] = strand
#		
#		for complex in self.three_arm_enumerator.initial_complexes:
#			self.complexes2[complex.name] = complex
			
	def testDomains(self):
		exp_domains = [Domain('1', 'short'), Domain('1', 'short', True), Domain('2', 'short'), Domain('2', 'short', True), Domain('3', 'short'), Domain('3', 'short', True), Domain('4', 'long'), Domain('4', 'long', True), Domain('5', 'short'), Domain('5', 'short', True), Domain('6', 'long'), Domain('6', 'long', True), Domain('7', 'short'), Domain('7', 'short', True)]
		act_domains = self.SLC_enumerator.domains
		
		exp_domains.sort()
		act_domains.sort()
		
		
		assert exp_domains == act_domains
		act_domains[0] = Domain('99', 'short')
		
		assert exp_domains != act_domains
		
		def assnDomains(self):
			self.SLC_enumerator.domains = []
			
		assert_raises(AttributeError, assnDomains, self)
		
	def testStrands(self):
		exp_strands = [Strand('PS', [self.domains['3*'], self.domains['2*'], self.domains['1*'], self.domains['5'], self.domains['6']]), Strand('OP', [self.domains['1'], self.domains['2'], self.domains['3'], self.domains['4']]), Strand('SP', [self.domains['5'], self.domains['6']]), Strand('BS', [self.domains['7*'], self.domains['6*'], self.domains['5*'], self.domains['1'], self.domains['2'], self.domains['3']]), Strand('Cat', [self.domains['6'], self.domains['7']])]
		act_strands = self.SLC_enumerator.strands
		
		exp_strands.sort()
		act_strands.sort()
		
		assert exp_strands == act_strands
		
		act_strands[0] = Strand('A', [self.domains['1']])
		
		assert exp_strands != act_strands
		
		def assnStrands(self):
			self.SLC_enumerator.domains = []
			
		assert_raises(AttributeError, assnStrands, self)
		
	def testInitialComplexes(self):
		PS = self.strands['PS']
		OP = self.strands['OP']
		Cat = self.strands['Cat']
		SP = self.strands['SP']
		BS = self.strands['BS']
	
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

		exp_complexes = [C1, C2, CatC, OPC, SPC, I1, I2, I3, I4, I5, I6, W]
		act_complexes = self.SLC_enumerator.initial_complexes
		
		exp_complexes.sort()
		act_complexes.sort()
		
		assert exp_complexes == act_complexes
		
		act_complexes[0] = Complex('A', [Cat], [[None, None]])
		
		assert exp_complexes != act_complexes
		
		def assnComplexes(self):
			self.SLC_enumerator.initial_complexes = []
			
		assert_raises(AttributeError, assnComplexes, self)
		
	def testSegmentNeighborhood1(self):
		enum = Enumerator(self.SLC_enumerator._domains, self.SLC_enumerator._strands, [self.complexes['Cat']])
		res = enum.segment_neighborhood(enum.initial_complexes, [])
	
		assert res == {'resting_states': [RestingState('0', enum.initial_complexes)], 'resting_state_complexes': enum.initial_complexes, 'transient_state_complexes': []}
				
	
	def testSegmentNeighborhood2(self):
		enum = self.SLC_enumerator
		complex_set = [self.complexes['Cat'], self.complexes['C2'], self.complexes['I1'], self.complexes['I2'], self.complexes['I3'], self.complexes['SP']]
		reaction_set = []
		# multi-molecular reactions never appear in segment neighborhood arguments
		#reaction_set.append(ReactionPathway('bind21', [self.complexes['Cat'], self.complexes['C2']], [self.complexes['I1']]))
		#reaction_set.append(ReactionPathway('bind21', [self.complexes['I3'], self.complexes['SP']], [self.complexes['I2']]))
		reaction_set.append(ReactionPathway('open', [self.complexes['I1']], [self.complexes['Cat'], self.complexes['C2']]))
		reaction_set.append(ReactionPathway('open', [self.complexes['I2']], [self.complexes['I3'], self.complexes['SP']]))
		reaction_set.append(ReactionPathway('branch_3way', [self.complexes['I1']], [self.complexes['I2']]))
		reaction_set.append(ReactionPathway('branch_3way', [self.complexes['I2']], [self.complexes['I1']]))
	
		res = enum.segment_neighborhood(complex_set, reaction_set)
		exp = {'resting_states': sorted([
										 RestingState('0', [self.complexes['Cat']]),
										 RestingState('1', [self.complexes['C2']]),
										 RestingState('2', [self.complexes['SP']]),
										 RestingState('3', [self.complexes['I3']])
										 ]), 
			'resting_state_complexes': sorted([self.complexes['Cat'], self.complexes['C2'], self.complexes['SP'], self.complexes['I3']]), 
			'transient_state_complexes': sorted([self.complexes['I1'], self.complexes['I2']])}
	
		assert res == exp
																									 
									
	def testSegmentNeighborhood3(self):
		enum = self.three_arm_enumerator
		complex_set = [self.complexes2['IABC'], self.complexes2['I'], self.complexes2['ABC']]
		reaction_set = []
		reaction_set.append(ReactionPathway('branch_3way', [self.complexes2['IABC']], [self.complexes2['I'], self.complexes2['IABC']]))
		
		res = enum.segment_neighborhood(complex_set, reaction_set)
		exp = {
			'resting_states': sorted([
									 RestingState('0', [self.complexes2['I']]),
									 RestingState('1', [self.complexes2['ABC']]) 
									  ]),
			'resting_state_complexes': sorted([self.complexes2['I'], self.complexes2['ABC']]),
			'transient_state_complexes': sorted([self.complexes2['IABC']])
			}
		assert res == exp
																									 
	def testProcessNeighborhood1(self):
		# regression
		pass
		# enum = Enumerator(self.SLC_enumerator._domains, self.SLC_enumerator._strands, [self.complexes['Cat']])
		# enum._E = []
		# enum._N = []
		# enum._S = []
		# enum._T = []
		# enum._B = []
		# enum._resting_states = []
		# enum._reactions = []
		# enum._complexes = []
		# enum.process_neighborhood(self.complexes['Cat'])
		# assert enum._S == [self.complexes['Cat']]
			
	def testProcessNeighborhood2(self):
		# regression
		pass
		# enum = Enumerator(self.SLC_enumerator._domains, self.SLC_enumerator._strands, [self.complexes['I1']])
		# enum._N = []
		# enum._S = []
		# enum._T = []
		# enum._E = []
		# enum._B = enum.initial_complexes
		# enum._resting_states = []
		# enum._reactions = []
		# enum._complexes = []
		# enum.process_neighborhood(self.complexes['I1'])
		# assert enum._S == sorted([self.complexes['Cat'], self.complexes['C2'], self.complexes['SP'], self.complexes['I3']])
		# assert enum._T == sorted([self.complexes['I1'], self.complexes['I2']])
		# assert enum._resting_states == sorted([
		# 									   RestingState('0', [self.complexes['Cat']]),
		# 									   RestingState('1', [self.complexes['C2']]),
		# 									   RestingState('2', [self.complexes['SP']]),
		# 									   RestingState('3', [self.complexes['I3']])
		# 									   ])
		# reaction_set = []
		# reaction_set.append(ReactionPathway('open', [self.complexes['I1']], [self.complexes['Cat'], self.complexes['C2']]))
		# reaction_set.append(ReactionPathway('open', [self.complexes['I2']], [self.complexes['I3'], self.complexes['SP']]))
		# reaction_set.append(ReactionPathway('branch_3way', [self.complexes['I1']], [self.complexes['I2']]))
		# reaction_set.append(ReactionPathway('branch_3way', [self.complexes['I2']], [self.complexes['I1']]))
		# reaction_set.sort()
		
		# assert sorted(enum._reactions) == reaction_set
		# assert enum._N == []
	
	def testProcessNeighborhood3(self):
		# regression
		pass
		# enum = Enumerator(self.three_arm_enumerator._domains, self.three_arm_enumerator._strands, [self.complexes2['IABC']])
		# enum._N = []
		# enum._S = []
		# enum._T = []
		# enum._E = []
		# enum._B = enum.initial_complexes
		# enum._resting_states = []
		# enum._reactions = []
		# enum._complexes = []
		# enum.process_neighborhood(self.complexes2['IABC'])
		# assert enum._S == sorted([self.complexes2['I'], self.complexes2['ABC']])
		# assert enum._T == sorted([self.complexes2['IABC']])
		# assert enum._resting_states	== sorted([
		# 									   RestingState('0', [self.complexes2['I']]),
		# 									   RestingState('1', [self.complexes2['ABC']]) 
		# 									   ])
		# reaction_set = []
		# reaction_set.append(ReactionPathway('branch_3way', [self.complexes2['IABC']], [self.complexes2['I'], self.complexes2['ABC']]))
		# assert sorted(enum._reactions) == sorted(reaction_set)


	def testEnumeration1(self):
		enum = input_enum('test_files/test_input_standard_simple.in')
		strands = {}
		for strand in enum.strands:
			strands[strand.name] = strand
		
		expected_complexes = []
		for complex in enum.initial_complexes:
			expected_complexes.append(complex)
			
		new_complex = Complex('C3', [strands['S1'], strands['S2']], [[(1, 0), None],[(0, 0)]])
		expected_complexes.append(new_complex)
		enum.enumerate()
		resting_complexes = enum.resting_complexes
		
		expected_complexes.sort()
		resting_complexes.sort()



		assert expected_complexes == resting_complexes
		
		
	def testEnumeration2(self):
		enum = Enumerator(self.SLC_enumerator._domains, self.SLC_enumerator._strands, [self.complexes['I1']])
		enum.enumerate()
		
		exp_initial_complexes = [self.complexes['I1']]
		assert exp_initial_complexes == enum.initial_complexes
			
		reaction_set = []
		reaction_set.append(ReactionPathway('bind21', [self.complexes['Cat'], self.complexes['C2']], [self.complexes['I1']]))
		reaction_set.append(ReactionPathway('bind21', [self.complexes['I3'], self.complexes['SP']], [self.complexes['I2']]))
		reaction_set.append(ReactionPathway('open', [self.complexes['I1']], [self.complexes['Cat'], self.complexes['C2']]))
		reaction_set.append(ReactionPathway('open', [self.complexes['I2']], [self.complexes['I3'], self.complexes['SP']]))
		reaction_set.append(ReactionPathway('branch_3way', [self.complexes['I1']], [self.complexes['I2']]))
		reaction_set.append(ReactionPathway('branch_3way', [self.complexes['I2']], [self.complexes['I1']]))
		reaction_set.sort()	
			
		exp_reactions = reaction_set
		res_reactions = sorted(enum.reactions)
			
		assert exp_reactions == res_reactions
			
		
		exp_resting_states = sorted([
									RestingState('0', [self.complexes['Cat']]),
									RestingState('1', [self.complexes['C2']]),
									RestingState('2', [self.complexes['SP']]),
									RestingState('3', [self.complexes['I3']])
									])
			
		assert exp_resting_states == enum.resting_states
			
		
		exp_complexes = sorted([self.complexes['I1'], self.complexes['I2'], self.complexes['Cat'], self.complexes['C2'], self.complexes['SP'], self.complexes['I3']])
		res_complexes = sorted(enum.complexes)
		
		assert res_complexes == exp_complexes
			
		exp_resting_complexes = sorted([self.complexes['Cat'], self.complexes['C2'], self.complexes['SP'], self.complexes['I3']])
		res_resting_complexes = sorted(enum.resting_complexes)
		assert res_resting_complexes == exp_resting_complexes
			
		exp_transient_complexes = sorted([self.complexes['I1'], self.complexes['I2']])
		res_transient_complexes = sorted(enum.transient_complexes)
		assert res_transient_complexes == exp_transient_complexes

	def testEnumeration3(self):
		enum = Enumerator(self.SLC_enumerator.domains, self.SLC_enumerator.strands, [self.complexes['Cat'], self.complexes['C1'], self.complexes['C2']])
		enum.enumerate()
		
		# shortcut because of large number of uses
		c = self.complexes

		exp_initial_complexes = sorted([self.complexes['Cat'], self.complexes['C1'], self.complexes['C2']])
		res_initial_complexes = sorted(enum.initial_complexes)
		
		assert exp_initial_complexes == res_initial_complexes

		reaction_set = []
		reaction_set.append(ReactionPathway('bind21', [c['Cat'], c['C2']], [c['I1']]))
		reaction_set.append(ReactionPathway('open', [c['I1']], [c['Cat'], c['C2']]))
		reaction_set.append(ReactionPathway('branch_3way', [c['I1']], [c['I2']]))
		reaction_set.append(ReactionPathway('branch_3way', [c['I2']], [c['I1']]))
		reaction_set.append(ReactionPathway('open', [c['I2']], [c['SP'], c['I3']]))
		reaction_set.append(ReactionPathway('bind21', [c['SP'], c['I3']], [c['I2']]))
		reaction_set.append(ReactionPathway('bind21', [c['I3'], c['C1']], [c['I4']]))
		reaction_set.append(ReactionPathway('open', [c['I4']], [c['I3'], c['C1']]))
		reaction_set.append(ReactionPathway('branch_3way', [c['I4']], [c['OP'], c['I5']]))

		I6 = copy.deepcopy(c['I5'])
		Cat_index = I6.strand_index('Cat')
		BS_index = I6.strand_index('BS')
		PS_index = I6.strand_index('PS')

		I6.structure[Cat_index][0] = None
		I6.structure[BS_index][1] = (PS_index, 4)
		I6.structure[PS_index][4] = (BS_index, 1)
		
		self.complexes['I6'] = I6
		I6._name = 'I6'
				
		reaction_set.append(ReactionPathway('branch_3way', [c['I5']], [c['I6']]))
		reaction_set.append(ReactionPathway('branch_3way', [c['I6']], [c['I5']]))
		reaction_set.append(ReactionPathway('open', [c['I6']], [c['W'], c['Cat']]))
		reaction_set.append(ReactionPathway('bind21', [c['W'], c['Cat']], [c['I6']]))

		I7 = copy.deepcopy(c['I4'])
		Cat_index = I7.strand_index('Cat')
		BS_index = I7.strand_index('BS')
		PS_index = I7.strand_index('PS')
	
		I7.structure[Cat_index][0] = None
		I7.structure[BS_index][1] = (PS_index, 4)
		I7.structure[PS_index][4] = (BS_index, 1)
	
		self.complexes['I7'] = I7
		I7._name = 'I7'
		
		reaction_set.append(ReactionPathway('branch_3way', [c['I4']], [c['I7']]))
		reaction_set.append(ReactionPathway('branch_3way', [c['I7']], [c['I4']]))
		reaction_set.append(ReactionPathway('branch_3way', [c['I7']], [c['OP'], c['I6']]))
		
		
		I8 = Complex('I8', [self.strands['BS'], self.strands['OP'], self.strands['PS']], 
							[[None, (2, 4), (2, 3), None, None, None], [(2, 2), (2, 1), (2, 0), None], 
							 [(1, 2), (1, 1), (1, 0), (0, 2), (0, 1)]])
		self.complexes['I8'] = I8
		I8._name = 'I8'

		reaction_set.append(ReactionPathway('open', [c['I7']], [c['I8'], c['Cat']]))
		reaction_set.append(ReactionPathway('branch_3way', [c['I8']], [c['OP'], c['W']]))
		
		
		exp_reactions = sorted(reaction_set)
		res_reactions = sorted(enum.reactions)
		
		assert exp_reactions == res_reactions

		exp_resting_states = sorted([RestingState('0', [self.complexes['Cat']]), 
							  RestingState('1', [self.complexes['C2']]),
							  RestingState('2', [self.complexes['SP']]),
							  RestingState('3', [self.complexes['I3']]),
							  RestingState('4', [self.complexes['C1']]),
							  RestingState('5', [self.complexes['OP']]),
							  RestingState('6', [self.complexes['W']])])
		res_resting_states = sorted(enum.resting_states)

		assert exp_resting_states == res_resting_states

		exp_resting_complexes = sorted([self.complexes['Cat'], self.complexes['C2'], self.complexes['SP'], self.complexes['I3'], self.complexes['C1'], self.complexes['OP'], self.complexes['W']])

		assert sorted(enum.resting_complexes) == exp_resting_complexes

		exp_transient_complexes = sorted([self.complexes['I1'], self.complexes['I2'], 
										  self.complexes['I4'], self.complexes['I5'], self.complexes['I6'],
										  self.complexes['I7'], self.complexes['I8']])
		res_transient_complexes	= sorted(enum.transient_complexes)

			
		assert res_transient_complexes == exp_transient_complexes

		assert sorted(exp_transient_complexes + exp_resting_complexes) == sorted(enum.complexes)
		
#	def testEnumeration4(self):
#		self.tet_enumerator = input_enum('test_files/examples/sadowski-tetrahedron.enum')
#		self.tet_enumerator.enumerate()
	
	def testEnumeration4(self):
		self.three_arm_nodal_enum = input_enum('test_files/examples/3-arm-junction.enum')
		self.three_arm_nodal_enum.enumerate()
	
	def testEnumeration5(self):
		self.seesaw_enumerator = input_enum('test_files/examples/seesaw/seesaw.enum')
		self.seesaw_enumerator.enumerate()

	def testEnumeration6(self):
		self.bounded_dendrimer = input_enum('test_files/examples/bounded-dendrimer.enum')

		self.bounded_dendrimer.MAX_COMPLEX_SIZE = 15	
		self.bounded_dendrimer.MAX_REACTION_COUNT = 1000
		self.bounded_dendrimer.MAX_COMPLEX_COUNT = 200
		self.bounded_dendrimer.RELEASE_CUTOFF = 8

		self.bounded_dendrimer.enumerate()
		
	def testEnumerationPolymer(self):
		# Test polymer detection
				
		
		# We're going to shrink some of these constants to trigger the exception
		# Test that too many reactions triggers exception
		polymer_enum = self.polymer_enum = input_enum('test_files/test_input_standard_polymer.in')
		polymer_enum.MAX_REACTION_COUNT = 10
		polymer_enum.enumerate()
		print "%d Complexes" % len(polymer_enum.complexes)
		print "%d Reactions" % len(polymer_enum.reactions)
		
		# The failure condition stops if .reactions > MAX_REACTION_COUNT, so we test for that (rather 
		# than that .reactions < MAX_REACTION_COUNT, which is not guaranteed.
		assert(len(polymer_enum.reactions) >= polymer_enum.MAX_REACTION_COUNT)
		
		# Now test that too many complexes also causes the error
		polymer_enum = self.polymer_enum = input_enum('test_files/test_input_standard_polymer.in')
		polymer_enum.MAX_COMPLEX_COUNT = 10
		polymer_enum.enumerate()
		print "%d Complexes" % len(polymer_enum.complexes)
		print "%d Reactions" % len(polymer_enum.reactions)
		# We're not examining len(polymer_enum.complexes) because that doesn't include ._S, which *is* 
		# tested for in the failure mode.
		assert((len(polymer_enum.complexes) + len(polymer_enum._S)) >= polymer_enum.MAX_COMPLEX_COUNT)
		
		complexes = polymer_enum._E + polymer_enum._T + polymer_enum._S
		assert max([len(c.strands) for c in complexes]) <= polymer_enum.MAX_COMPLEX_SIZE
		
		# Now we want to make sure that no reactions in the enumerator point to complexes that weren't in the list
		undefined_complexes = []
		
		for reaction in polymer_enum.reactions:
			for product in reaction.products:
				if not (product in polymer_enum.complexes):
					undefined_complexes.append(product)
					print "Reaction: %s, Product: %s" % (repr(reaction), repr(product))
		
		#undefined_complexes = [product for reaction in self.polymer_enum.reactions for product in reaction.products if not (product in complexes) ]
		
		print "Undefined complexes:"
		print undefined_complexes
		assert len(undefined_complexes) == 0
	