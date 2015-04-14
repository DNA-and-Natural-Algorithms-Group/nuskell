#
#  test_output.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 2/6/11.
#

import unittest
from utils import *
from enumerator import *
from nose.tools import *
from input import *
from output import *
import copy
import filecmp

class OutputTests(unittest.TestCase):
	# Disable output tests until all other functionality is working!
	
	# __test__ = False
	
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

		self.SLC_enumerator_reduced = Enumerator(self.SLC_enumerator.domains, self.SLC_enumerator.strands, [self.complexes['Cat'], self.complexes['C1'], self.complexes['C2']])
			
		self.three_arm_enumerator = input_enum('test_files/test_input_standard_3arm_junction.in')
		
				
		for complex in self.three_arm_enumerator.initial_complexes:
			self.complexes[complex.name] = complex
				
		self.three_arm_enumerator_reduced = Enumerator(self.three_arm_enumerator.domains, self.three_arm_enumerator.strands, [self.complexes['I'], self.complexes['A'], self.complexes['B'], self.complexes['C']])

		self.simple_enumerator = input_enum('test_files/test_input_standard_simple.in')
				
		self.simple2_enumerator = input_enum('test_files/test_input_standard_simple2.in')
		
		
		self.SLC_enumerator.enumerate()
		self.SLC_enumerator_reduced.enumerate()
		# self.three_arm_enumerator.enumerate()
		self.three_arm_enumerator_reduced.enumerate()
		self.simple_enumerator.enumerate()
		self.simple2_enumerator.enumerate()
		
	def testOutputLegacy(self):
		output_legacy(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_Legacy.out', output_condensed=False)

	def testOutputLegacy2(self):
		output_legacy(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_Legacy2.out', output_condensed=False)
		
	def testOutputLegacy3(self):
		output_legacy(self.simple2_enumerator, 'test_files/test_output/test_output_Legacy3.out', output_condensed=False)
	
	# def testOutputFullGraph(self):
	# 	output_full_graph(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_FullGraph.out')
	
	# def testOutputFullGraph2(self):
	# 	output_full_graph(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_FullGraph2.out')
	
	# def testOutputFullGraph3(self):
	# 	output_full_graph(self.simple2_enumerator, 'test_files/test_output/test_output_FullGraph3.out')

	# def testOutputFullGraph4(self):
	# 	output_full_graph(self.SLC_enumerator, 'test_files/test_output/test_output_FullGraph4.out')
		
	# def testOutputCondensedGraph(self):
	# 	output_condensed_graph(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_CondensedGraph.out')
	
	# def testOutputCondensedGraph2(self):
	# 	output_condensed_graph(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_CondensedGraph2.out')
	
	# def testOutputCondensedGraph3(self):
	# 	output_condensed_graph(self.simple2_enumerator, 'test_files/test_output/test_output_CondensedGraph3.out')
	
	# def testOutputCondensedGraph4(self):
	# 	output_condensed_graph(self.SLC_enumerator, 'test_files/test_output/test_output_CondensedGraph4.out')
		
	# def testOutputGraph(self):
	# 	output_graph(self.simple2_enumerator,'test_files/test_output/test_output_Graph.out',output_condensed=True)
	# 	output_graph(self.simple2_enumerator,'test_files/test_output/test_output_GraphCondensed.out',output_condensed=False)


	# PIL output
	def testOutputPIL(self):
		output_pil(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_FullGraph.pil',output_condensed=False)
	
	def testOutputPIL2(self):
		output_pil(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_FullGraph2.pil',output_condensed=False)
	
	def testOutputPIL3(self):
		output_pil(self.simple2_enumerator, 'test_files/test_output/test_output_FullGraph3.pil',output_condensed=False)
	
	def testOutputPIL4(self):
		output_pil(self.SLC_enumerator, 'test_files/test_output/test_output_FullGraph4.pil',output_condensed=False)
			

	# Condensed PIL output	
	def testOutputCondensedPIL(self):
		output_pil(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_CondensedGraph.pil',output_condensed=True)
	
	def testOutputCondensedPIL2(self):
		output_pil(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_CondensedGraph2.pil',output_condensed=True)
	
	def testOutputCondensedPIL3(self):
		output_pil(self.simple2_enumerator, 'test_files/test_output/test_output_CondensedGraph3.pil',output_condensed=True)
	
	def testOutputCondensedPIL4(self):
		output_pil(self.SLC_enumerator, 'test_files/test_output/test_output_CondensedGraph4.pil',output_condensed=True)
	

	# CRN
	def testOutputCRN(self):
		output_crn(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_FullGraph.crn',output_condensed=False)
		filecmp.cmp('test_files/test_output/test_output_FullGraph.crn','test_files/expected_output/test_output_FullGraph.crn')

	def testOutputCRN2(self):
		output_crn(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_FullGraph2.crn',output_condensed=False)
		filecmp.cmp('test_files/test_output/test_output_FullGraph2.crn','test_files/expected_output/test_output_FullGraph2.crn')
	
	def testOutputCRN3(self):
		output_crn(self.simple2_enumerator, 'test_files/test_output/test_output_FullGraph3.crn',output_condensed=False)
		filecmp.cmp('test_files/test_output/test_output_FullGraph3.crn','test_files/expected_output/test_output_FullGraph3.crn')
	
	def testOutputCRN4(self):
		output_crn(self.SLC_enumerator, 'test_files/test_output/test_output_FullGraph4.crn',output_condensed=False)
		filecmp.cmp('test_files/test_output/test_output_FullGraph4.crn','test_files/expected_output/test_output_FullGraph4.crn')


	def testOutputCRN_condensed(self):
		output_crn(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_CondensedGraph.crn',output_condensed=True)
		filecmp.cmp('test_files/test_output/test_output_CondensedGraph.crn','test_files/expected_output/test_output_CondensedGraph.crn')

	def testOutputCRN_condensed2(self):
		output_crn(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_CondensedGraph2.crn',output_condensed=True)
		filecmp.cmp('test_files/test_output/test_output_CondensedGraph2.crn','test_files/expected_output/test_output_CondensedGraph2.crn')
	
	def testOutputCRN_condensed3(self):
		output_crn(self.simple2_enumerator, 'test_files/test_output/test_output_CondensedGraph3.crn',output_condensed=True)
		filecmp.cmp('test_files/test_output/test_output_CondensedGraph3.crn','test_files/expected_output/test_output_CondensedGraph3.crn')
	
	def testOutputCRN_condensed4(self):
		output_crn(self.SLC_enumerator, 'test_files/test_output/test_output_CondensedGraph4.crn',output_condensed=True)
		filecmp.cmp('test_files/test_output/test_output_CondensedGraph4.crn','test_files/expected_output/test_output_CondensedGraph4.crn')

	# K
	def testOutputK(self):
		output_k(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_FullGraph.k',output_condensed=False)
		filecmp.cmp('test_files/test_output/test_output_FullGraph.k','test_files/expected_output/test_output_FullGraph.k')

	def testOutputK2(self):
		output_k(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_FullGraph2.k',output_condensed=False)
		filecmp.cmp('test_files/test_output/test_output_FullGraph2.k','test_files/expected_output/test_output_FullGraph2.k')
	
	def testOutputK3(self):
		output_k(self.simple2_enumerator, 'test_files/test_output/test_output_FullGraph3.k',output_condensed=False)
		filecmp.cmp('test_files/test_output/test_output_FullGraph3.k','test_files/expected_output/test_output_FullGraph3.k')
	
	def testOutputK4(self):
		output_k(self.SLC_enumerator, 'test_files/test_output/test_output_FullGraph4.k',output_condensed=False)
		filecmp.cmp('test_files/test_output/test_output_FullGraph4.k','test_files/expected_output/test_output_FullGraph4.k')


	def testOutputK_condensed(self):
		output_k(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_CondensedGraph.k',output_condensed=True)
		filecmp.cmp('test_files/test_output/test_output_CondensedGraph.k','test_files/expected_output/test_output_CondensedGraph.k')

	def testOutputK_condensed2(self):
		output_k(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_CondensedGraph2.k',output_condensed=True)
		filecmp.cmp('test_files/test_output/test_output_CondensedGraph2.k','test_files/expected_output/test_output_CondensedGraph2.k')
	
	def testOutputK_condensed3(self):
		output_k(self.simple2_enumerator, 'test_files/test_output/test_output_CondensedGraph3.k',output_condensed=True)
		filecmp.cmp('test_files/test_output/test_output_CondensedGraph3.k','test_files/expected_output/test_output_CondensedGraph3.k')
	
	def testOutputK_condensed4(self):
		output_k(self.SLC_enumerator, 'test_files/test_output/test_output_CondensedGraph4.k',output_condensed=True)
		filecmp.cmp('test_files/test_output/test_output_CondensedGraph4.k','test_files/expected_output/test_output_CondensedGraph4.k')



	# JSON output
		
	def testOutputJSON(self):
		output_json(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_JSON.json')
		enumerator = load_json('test_files/test_output/test_output_JSON.json')
		assert enumerator == self.SLC_enumerator_reduced
		
	def testOutputJSON2(self):
		output_json(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_JSON2.json')
		enumerator = load_json('test_files/test_output/test_output_JSON2.json')
		assert enumerator == self.three_arm_enumerator_reduced
		
	def testOutputJSONCondensed(self):
		output_json(self.SLC_enumerator_reduced, \
				'test_files/test_output/test_output_JSONCondensed.json',output_condensed=True)
		enumerator = load_json('test_files/test_output/test_output_JSON.json')
		assert enumerator == self.SLC_enumerator_reduced
		
	def testOutputJSONCondensed2(self):
		output_json(self.three_arm_enumerator_reduced, \
				'test_files/test_output/test_output_JSONCondensed2.json',output_condensed=True)
		enumerator = load_json('test_files/test_output/test_output_JSON2.json')
		assert enumerator == self.three_arm_enumerator_reduced
			
	def testJSONInputOutputLoop(self):

		output_json(self.SLC_enumerator_reduced, 'test_files/test_output/testJSONInputOutputLoop.json')
		enumerator = load_json('test_files/test_output/testJSONInputOutputLoop.json')
		assert enumerator == self.SLC_enumerator_reduced

	def testOutputSBML(self):
		output_sbml(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_SBML.sbml')
		
	def testOutputSBML2(self):
		output_sbml(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_SBML2.sbml')
		
	def testOutputSBMLCondensed(self):
		output_sbml(self.SLC_enumerator_reduced, 'test_files/test_output/test_output_SBMLCondensed.sbml',output_condensed=True)
		
	def testOutputSBML2Condensed(self):
		output_sbml(self.three_arm_enumerator_reduced, 'test_files/test_output/test_output_SBML2Condensed.sbml',output_condensed=True)
			
	def testOutputLegacyCondensed(self):
		output_legacy(self.simple_enumerator, 'test_files/test_output/test_output_LegacyCondensed.out', output_condensed=True)
		
		
	def testCondenseRestingStates(self):
		# for enum in [self.SLC_enumerator, self.three_arm_enumerator]:
		for enum in [self.SLC_enumerator]:
			enum.enumerate()
			condensed = condense_resting_states(enum)
				
			assert (sorted(condensed['resting_states']) == \
					sorted(enum.resting_states))
			
			# Test that transient complexes have been eliminated from complexes
			resting_complexes = set([complex for state in condensed['resting_states'] for complex in state.complexes])
			predicted_resting_complexes = set(enum.complexes) - set(enum.transient_complexes)
			
			print resting_complexes
			print predicted_resting_complexes
			
			assert (resting_complexes == predicted_resting_complexes)
			
			# Test transient complexes have been removed from reactions
			transients = set(enum.transient_complexes)
			for reaction in condensed['reactions']:
				for reactant in reaction.reactants:
					assert reactant not in transients
					
				for product in reaction.products:
					assert product not in transients
					
			# Test that no reactions are duplicated
			assert (sorted(condensed['reactions']) == sorted(list(set(condensed['reactions']))))
	
	def testCondenseRestingStates1(self):
		condensed = condense_resting_states(self.simple2_enumerator)
		
		enum = self.simple2_enumerator
#		domains = {}
#		strands = {}
#		complexes = {}
#		
#		for domain in enum.domains:
#			domains[domain.name] = domain
#		
#		for strand in enum.strands:
#			strands[strand.name] = strand
#		
#		for complex in enum.initial_complexes:
#			complexes[complex.name] = complex
		
		(domains,strands,complexes) = index_parts(enum)
		
#		# End-state Complexes 
#		structure RC1 = S1 + S3 : ((+))
#		structure RC2 = S2 : .
#		structure C1 = S1 + S2 : .(+)
#		structure C2 = S3 : ..
#		
#		# Resting-state sets 
#		state RS_RC1 = RC1
#		state RS_RC2 = RC2
#		state RS_C1 = C1
#		state RS_C2 = C2
#		
#		# Condensed Reactions 
#		kinetic C1 + C2 -> RC1 + RC2 

		RC1 = Complex('RC1',[strands['S1'], strands['S3']],parse_dot_paren('((+))'))
		RC2 = Complex('RC2',[strands['S2']],parse_dot_paren('.'))
		
		RS_RC1 = RestingState('RS_RC1',[RC1])
		RS_RC2 = RestingState('RS_RC2',[RC2])
		RS_C1 = RestingState('RS_C1',[complexes['C1']])
		RS_C2 = RestingState('RS_C2',[complexes['C2']])
		
		expected_reactions = [ReactionPathway('condensed',[RS_C1,RS_C2],[RS_RC1,RS_RC2])]
		
		print "Expected reactions"
		print expected_reactions
		print "Actual reactions"
		print condensed['reactions']
		
		
		assert condensed['reactions'] == expected_reactions
		assert sorted([RS_RC1,RS_RC2,RS_C1,RS_C2]) == sorted(condensed['resting_states'])
				
				
	def testCondenseRestingStates2(self):
		enum = self.SLC_enumerator
		condensed = condense_resting_states(enum)
		
#		domains = {}
#		strands = {}
#		complexes = {}
#		
#		for domain in enum.domains:
#			domains[domain.name] = domain
#		
#		for strand in enum.strands:
#			strands[strand.name] = strand
#		
#		for complex in enum.initial_complexes:
#			complexes[complex.name] = complex

		(domains,strands,complexes) = index_parts(enum)
	
#		state RS_C1 = C1
#		state RS_C2 = C2
#		state RS_Cat = Cat
#		state RS_I3 = I3
#		state RS_OP = OP
#		state RS_SP = SP
#		state RS_W = W

#		kinetic I3 + C1 -> W + OP + Cat 
#		kinetic I3 + SP -> C2 + Cat 
#		kinetic C2 + Cat -> I3 + SP

		RS_C1 = RestingState('RS_C1',[complexes['C1']])
		RS_C2 = RestingState('RS_C2',[complexes['C2']])
		RS_Cat = RestingState('RS_Cat',[complexes['Cat']])
		RS_I3 = RestingState('RS_I3',[complexes['I3']])
		RS_OP = RestingState('RS_OP',[complexes['OP']])
		RS_SP = RestingState('RS_SP',[complexes['SP']])
		RS_W = RestingState('RS_W',[complexes['W']])
		
		expected_reactions = [ReactionPathway('condensed',[RS_I3, RS_C1],[RS_W,RS_OP,RS_Cat]), \
								ReactionPathway('condensed',[RS_I3, RS_SP],[RS_C2,RS_Cat]),
								ReactionPathway('condensed',[RS_C2, RS_Cat],[RS_I3,RS_SP])]
		
		print "Expected reactions"
		print expected_reactions
		print "Actual reactions"
		print condensed['reactions']
		
		assert sorted(condensed['reactions']) == sorted(expected_reactions)
		assert sorted([RS_C1,RS_C2,RS_Cat,RS_I3,RS_OP,RS_SP,RS_W]) == sorted(condensed['resting_states'])

