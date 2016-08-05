#
#  output.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/21/10.
#

import copy
import json
import subprocess
import collections

import nuskell.include.peppercorn.utils as utils
import nuskell.include.peppercorn.reactions as reactions
from nuskell.include.peppercorn.reactions import ReactionPathway, auto_name
from nuskell.include.peppercorn.condense import condense_resting_states



def output_legacy(enumerator, filename, output_condensed = False, output_rates = False, condense_options = {}):
	"""
	Legacy text-based output scheme similar to that used in the original 
	enumerator. Designed to be simultaneously parsable and human-readable.
	Supports output of condensed graph in addition to the full graph. Does
	not support strands.
	"""
	
	def write_complex(output_file,complex):
		output_file.write(str(complex) + "\n")
		names = []
		for strand in complex.strands:
			names.append(" ".join(map(str,strand.domains)))
		strands_string = " + ".join(names)
		output_file.write(strands_string + "\n")
		output_file.write(str(complex.dot_paren_string()) + "\n")
		output_file.write("\n")
	
	def write_reaction(output_file,reaction):
		reactants = map(str,reaction.reactants)
		products = map(str,reaction.products)
		reac_string_list = [" + ".join(reactants),"->"," + ".join(products),"\n"]
		reac_string = ' '.join(reac_string_list)
		output_file.write(reac_string)
	
#	def write_reaction(output_file,reaction):
#		reactants = reaction.reactants
#		products = reaction.products
#		reac_string_list = [reactants[0].name]
#		for reactant in reactants[1:]:
#			reac_string_list.append(" + " + reactant.name)
#		reac_string_list.append(" -> ")
#		reac_string_list.append(products[0].name)
#		for product in products[1:]:
#			reac_string_list.append(" + " + product.name)
#		reac_string_list.append("\n")
#		reac_string = ''.join(reac_string_list)
#		output_file.write(reac_string)
		
	complexes = enumerator.complexes
	transient_complexes = enumerator.transient_complexes
	resting_complexes = enumerator.resting_complexes
	reactions = enumerator.reactions
	resting_states = enumerator.resting_states
	
	output_file = open(filename, 'w')
	output_file.write("###### Enumerated Output ######\n")
	output_file.write("\n\n# Domains \n")
	for domain in sorted(enumerator.domains):
		if(not domain.is_complement):
			output_file.write("sequence " + domain.name + " = : " + str(domain.length) + "\n")
	
	output_file.write("\n\n# End-state Complexes \n")
	for complex in sorted(resting_complexes):
		write_complex(output_file,complex)
		
	# Not part of the original output; omitting
#	output_file.write("###############################\n")
#	output_file.write("\n\n# Resting-state sets \n")
#	for resting_state in sorted(resting_states):
#		output_file.write("state " + resting_state.name + " = : ")
#		output_file.write(str(resting_state.complexes[0]))
#		for complex in resting_state.complexes[1:]:
#			output_file.write(" %s" % complex)
#		output_file.write("\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Fast (Transition) Complexes \n")
	for complex in sorted(transient_complexes):
		write_complex(output_file,complex)
	output_file.write("###############################\n")
	output_file.write("\n\n# Reactions \n")
	for reaction in sorted(reactions):
		write_reaction(output_file,reaction)
		
	if (output_condensed):
		output_file.write("###############################\n")
		output_file.write("\n\n# Condensed Reactions \n")
		condensed = condense_resting_states(enumerator, **condense_options)
		new_reactions = condensed['reactions']
		for reaction in sorted(new_reactions):
			write_reaction(output_file,reaction)
				
	output_file.close()

output_legacy.supports_condensed = True		

def output_pil(enumerator, filename, output_condensed = False, output_rates = True, condense_options = {}):
	"""
	Text-based output using the Pepper Intermediate Language (PIL)
	"""
	
	def write_complex(output_file,complex):
		params = ""
		if complex.concentration is not None:
			params = "[@ %g %sM] " % utils.format_si(complex.concentration)
		output_file.write("structure " + params + str(complex) + " = ")
		names = map(lambda strand: strand.name, complex.strands)
		strands_string = " + ".join(names)
		output_file.write(strands_string + " : ")
		output_file.write(str(complex.dot_paren_string()) + "\n")
	
	def write_reaction(output_file,reaction):
		reactants = map(str,reaction.reactants)
		products = map(str,reaction.products)

		if output_rates:
			rate_units = "/M" * (reaction.arity[0]-1) + "/s"
			rate_const = "[%12.8g %s]" % (reaction.rate(), rate_units) 
		else: rate_const = ""

		reac_string_list = ["kinetic",rate_const," + ".join(reactants),"->"," + ".join(products),"\n"]
		reac_string = ' '.join(reac_string_list)
		output_file.write(reac_string)
		
	complexes = enumerator.complexes
	transient_complexes = enumerator.transient_complexes
	resting_complexes = enumerator.resting_complexes
	reactions = enumerator.reactions
	resting_states = enumerator.resting_states
	
	output_file = open(filename, 'w')
	output_file.write("###### Enumerated Output ######\n")
	output_file.write("\n# Domains \n")
	
	def seq(dom):
		if(dom.sequence != None):
			return dom.sequence
		else:
			return "N" * len(dom)
	
	for domain in utils.natural_sort(enumerator.domains):
		if(not domain.is_complement):
			output_file.write("sequence " + domain.name + " = " + seq(domain) + " : " + str(domain.length) + "\n")
	
	output_file.write("\n# Strands \n")
	for strand in utils.natural_sort(enumerator.strands):
		output_file.write("strand " + strand.name + " = " + \
						" ".join(map(lambda dom: dom.name, strand.domains)) + "\n")
	
	output_file.write("\n# Resting-state Complexes \n")
	for complex in utils.natural_sort(resting_complexes):
		write_complex(output_file,complex)
		
	

	
	if (output_condensed):
		condensed = condense_resting_states(enumerator, **condense_options)

		output_file.write("\n# Resting-state sets \n")
		resting_states = condensed['resting_states']
		for resting_state in utils.natural_sort(resting_states):
			output_file.write("# state " + str(resting_state) + " = { " + " ".join(map(str,resting_state.complexes)) + " }\n")

		output_file.write("\n# Condensed Reactions \n")
		new_reactions = condensed['reactions']
		for reaction in sorted(new_reactions):
			write_reaction(output_file,reaction)
	else:	
		output_file.write("\n# Transient Complexes \n")
		for complex in utils.natural_sort(transient_complexes):
			write_complex(output_file,complex)
		output_file.write("\n# Detailed Reactions \n")
		for reaction in sorted(reactions): #utils.natural_sort(reactions):
			write_reaction(output_file,reaction)
		
			
	output_file.close()


def output_kernel(enumerator, filename, output_condensed = False, output_rates = True, condense_options = {}):
	"""
	Text-based output using the Pepper Intermediate Language (PIL)
	"""
	
	def write_complex(output_file,complex):
		params = ""
		if complex.concentration is not None:
			params = "[@ %g %sM] " % utils.format_si(complex.concentration)
		output_file.write(str(complex) + params + " = " + complex.kernel_string() + "\n")
	
	def write_reaction(output_file,reaction):
		if output_rates:
			rate_units = "/M" * (reaction.arity[0]-1) + "/s"
			rate_const = "[%f %s]" % (reaction.rate(), rate_units) 
		else: rate_const = ""
		reac_string_list = [rate_const,reaction.kernel_string()]
		reac_string = ' '.join(reac_string_list)

		reactants = map(str,reaction.reactants)
		products = map(str,reaction.products)
		prefix = " ".join(["# ",str(reaction)," + ".join(reactants),"->"," + ".join(products)])

		output_file.write(prefix + "\n" + reac_string + "\n")
		
	complexes = enumerator.complexes
	transient_complexes = enumerator.transient_complexes
	resting_complexes = enumerator.resting_complexes
	reactions = enumerator.reactions
	resting_states = enumerator.resting_states
	
	output_file = open(filename, 'w')
	output_file.write("###### Enumerated Output ######\n")
	output_file.write("\n# Domains \n")
	
	def seq(dom):
		if(dom.sequence != None):
			return dom.sequence
		else:
			return "N" * len(dom)
	
	for domain in utils.natural_sort(enumerator.domains):
		if(not domain.is_complement):
			output_file.write("sequence " + domain.name + " = " + seq(domain) + " : " + str(domain.length) + "\n")
	
	output_file.write("\n# Strands \n")
	for strand in utils.natural_sort(enumerator.strands):
		output_file.write("strand " + strand.name + " = " + \
						" ".join(map(lambda dom: dom.name, strand.domains)) + "\n")
	
	output_file.write("\n# Resting-state Complexes \n")
	for complex in utils.natural_sort(resting_complexes):
		write_complex(output_file,complex)

	
	if (output_condensed):
		condensed = condense_resting_states(enumerator, **condense_options)

		output_file.write("\n# Resting-state sets \n")
		resting_states = condensed['resting_states']
		for resting_state in utils.natural_sort(resting_states):
			output_file.write("# state " + str(resting_state) + " = { " + " ".join(map(str,resting_state.complexes)) + " }\n")

		output_file.write("\n# Condensed Reactions \n")
		new_reactions = condensed['reactions']
		for reaction in sorted(new_reactions):
			write_reaction(output_file,reaction)
	else:	
		output_file.write("\n# Transient Complexes \n")
		for complex in utils.natural_sort(transient_complexes):
			write_complex(output_file,complex)
		output_file.write("\n# Detailed Reactions \n")
		for reaction in sorted(reactions): #utils.natural_sort(reactions):
			write_reaction(output_file,reaction)
		
			
	output_file.close()


def output_json(enumerator, filename, output_condensed = False, output_rates = True, condense_options = {}):
	"""
	JSON-based output schema intended to be easily machine parsable. Uses
	python's JSON serialization libraries.
	"""
	
	def serializeComplex(complex):
		temp_strands = []
		for strand in complex.strands:
			temp_strands.append(strand.name)
		return { 
			'name':complex.name,
			'strands': temp_strands,
			'structure': complex.structure,
			'dot-paren': complex.dot_paren_string(),
			'dot-paren-full': complex.dot_paren_string_full(),
		}
	
	def serializeReaction(reaction):
		# temp_reactants = []
		# for reactant in reaction.reactants:
		# 	temp_reactants.append(reactant.name)
		# temp_products = []
		# for product in reaction.products:
		# 	temp_products.append(product.name)
		reactants = map(str,reaction.reactants)
		products = map(str,reaction.products)
		return {
			"name":reaction.name,
			"reactants": reactants,  #temp_reactants,
			"products": products #temp_products
		}
	
	def serializeDomain(domain):
		temp_domain = {
			"name": domain.name,
			"length": domain.length,
			"is_complement": domain.is_complement
		}

		# temp_domain = {}
		# temp_domain['name'] = domain.name
		# temp_domain['length'] = domain.length
		# temp_domain['is_complement'] = domain.is_complement
		if domain.sequence != None:
			temp_domain['sequence'] = domain.sequence
		return temp_domain
		
	def serializeStrand(strand):
		return {
			"name": strand.name,
			"domains": [domain.name for domain in strand.domains]
		}
		# temp_strand = {}
		# temp_strand['name'] = strand.name
		# temp_domains = []
		# for domain in strand.domains:
		# 	temp_domains.append(domain.name)
		# temp_strand['domains'] = temp_domains
		# return temp_strand
		
	def serializeRestingState(resting_state):
		return {
			"name": str(resting_state),
			"complexes": [complex.name for complex in resting_state.complexes]
		}
		# temp_resting_state = {}
		# temp_complexes = []
		# for complex in resting_state.complexes:
		# 	temp_complexes.append(complex.name)
		# temp_resting_state['name'] = resting_state.name
		# temp_resting_state['complexes'] = temp_complexes
		# return temp_resting_state
		
	object_out = {}
	
	object_out['domains'] = map(serializeDomain,enumerator.domains)
	object_out['strands'] = map(serializeStrand,enumerator.strands)
	object_out['resting_complexes'] = map(serializeComplex,enumerator.resting_complexes)
	object_out['transient_complexes'] = map(serializeComplex,enumerator.transient_complexes)
	object_out['resting_states'] = map(serializeRestingState,enumerator.resting_states)
	object_out['initial_complexes'] = map(serializeComplex,enumerator.initial_complexes)
	object_out['reactions'] = map(serializeReaction,enumerator.reactions)
	
	if output_condensed:
		condensed = condense_resting_states(enumerator, **condense_options)
		object_out['condensed_reactions'] = map(serializeReaction,condensed['reactions'])
		
	fout = open(filename, 'w')
	json.dump(object_out, fout, indent=4)
	fout.close()
 
output_json.supports_condensed = True

def output_graph(enumerator, filename, output_condensed = False, output_rates = False):
	if not output_condensed:
		output_full_graph(enumerator, filename)
	else:
		output_condensed_graph(enumerator, filename)
	
output_graph.supports_condensed = True

def output_dotfile(enumerator, filename, output_condensed = False, output_rates = False):
	if not output_condensed:
		output_full_dotfile(enumerator, filename)
	else:
		output_condensed_dotfile(enumerator, filename)
	
output_dotfile.supports_condensed = True

def output_full_dotfile(enumerator, filename):
	fout = open(filename + ".dot", "w")
	fout.write("digraph G {\n")
	fout.write('size="7,10"\n')
	fout.write('page="8.5,11"\n')
	fout.write('node[width=0.25,height=0.375,fontsize=9]\n')
	
	# We need to determine the complex clusters for the graph. Complexes with
	# the same cyclic permutation of strands are placed in the same cluster.
	
	strand_cyclic_permutations = []
	clusters = []
	
	# We loop through all complexes
	for complex in enumerator.transient_complexes:
		complex._resting = False
		flag = False
		
		# Check to see if we've already seen this strand ordering
		for (i, perm) in enumerate(strand_cyclic_permutations):
			if perm == complex.strands:
				clusters[i].append(complex)
				flag = True
				break
				
		# If not, we add it
		if not flag:
			strand_cyclic_permutations.append(complex.strands)
			clusters.append([complex])
			
	for complex in enumerator.resting_complexes:
		complex._resting = True
		flag = False
		
		# Check to see if we've already seen this strand ordering
		for (i, perm) in enumerate(strand_cyclic_permutations):
			if perm == complex.strands:
				clusters[i].append(complex)
				flag = True
				break
				
		# If not, we add it
		if not flag:
			strand_cyclic_permutations.append(complex.strands)
			clusters.append([complex])
	
	# We now draw the clusters on the graph
	for (i, cluster) in enumerate(clusters):
		fout.write("subgraph cluster%d {\n" % i)
		strands = [cluster[0].strands[0].name]
		for strand in cluster[0].strands[1:]:
			strands.append(" + ")
			strands.append(strand.name)
		strands_string = ''.join(strands)
		fout.write('label="%s"\n' % strands_string)
		fout.write('fontsize=6\n')
		for complex in cluster:
			extra_params = ""
			if complex._resting:
				extra_params = ",style=filled,color=gold1"
			fout.write('%s [label="%s: %s"%s];\n' % (str(complex), 
													str(complex), 
													complex.dot_paren_string(), 
													extra_params
													))
		fout.write("}\n")
		
	# We now draw the reactions. If there is one product and one reagent, then
	# we just draw an edge between the two. Otherwise we create a reaction node.
	for (i, reaction) in enumerate(enumerator.reactions):
		if (len(reaction.products) == 1) and (len(reaction.reactants) == 1):
			fout.write("%s -> %s\n" % (str(reaction.reactants[0]),
									   str(reaction.products[0])))
		else:
			reaction_label = "R_%d" % i
			# We label unimolecular reactions blue, and other reactions red
			reaction_color = "red"
			if len(reaction.reactants) == 1:
				reaction_color = "blue"
			fout.write('%s [label="",shape=circle,height=0.12,width=0.12,fontsize=1,style=filled,color=%s];\n' %
							(reaction_label, reaction_color))
			
			# We now make all the edges needed
			for reagent in reaction.reactants:
				fout.write("%s -> %s\n" % (str(reagent), reaction_label))
				
			for product in reaction.products:
				fout.write("%s -> %s\n" % (reaction_label, str(product)))
				
				
	fout.write("}\n")
	fout.close()
	


def output_full_graph(enumerator, filename):
	"""
	Products graphical output representing the full reaction graph, with all
	reactions and complexes. Transient and resting states are differentiated
	by color.
	"""
	
	output_full_dotfile(enumerator,filename)

	# Create the output file.
	# TODO: make 'pdf' configurable
	subprocess.call(["dot", "-O", "-Teps", "%s.dot" % filename])


def output_condensed_dotfile(enumerator, filename):
	fout = open(filename + ".dot", "w")
	fout.write("digraph G {\n")
	fout.write('size="7,10"\n')
	fout.write('page="8.5,11"\n')
	fout.write('node[width=0.25,height=0.375,fontsize=9]\n')
	
	condensed_graph = condense_resting_states(enumerator, **condense_options)
	
	# We loop through all resting states, drawing them on the graph
	for state in condensed_graph['resting_states']:
		fout.write('%s [label="%s"]\n' % (str(state), str(state)))
	
	# We now add reactions. We always create a reaction node.
	for (i, reaction) in enumerate(condensed_graph['reactions']):
		reaction_label = "R_%d" % i
		fout.write('%s [label="",shape=circle,height=0.12,width=0.12,fontsize=1,style=filled,color=red];\n' % reaction_label)
		
		for reagent in reaction.reactants:
			fout.write("%s -> %s\n" % (str(reagent), reaction_label))
			
		for product in reaction.products:
			fout.write("%s -> %s\n" % (reaction_label, str(product)))
			
	fout.write("}\n")
	fout.close()

def output_condensed_graph(enumerator, filename):
	"""
	Products graphical output representing the condensed reaction graph, with
	condensed reactions and resting states aggregated into single nodes.
	"""
	
	output_condensed_dotfile(enumerator,filename)
	
	# Create the output file.
	# TODO: make 'pdf' configurable
	subprocess.call(["dot", "-O", "-Teps", "%s.dot" % filename])
	

def output_sbml(enumerator,filename, output_condensed = False, output_rates = True, condense_options = {}):
	# # default initial concentration of all species is 100 nM
	# initial_concentration = 10e-7 

	import xml.dom.minidom
	header = '<?xml version="1.0" encoding="UTF-8"?>'
	out = [header,
		'<sbml level="2" version="3" xmlns="http://www.sbml.org/sbml/level2/version3">',
		'<model name="%s">' % filename,
		'<listOfUnitDefinitions>',
            '<unitDefinition id="per_second">',
                '<listOfUnits>',
                    '<unit kind="second" exponent="-1"/>',
               ' </listOfUnits>',
            '</unitDefinition>',
            '<unitDefinition id="litre_per_mole_per_second">',
                '<listOfUnits>',
                    '<unit kind="mole"   exponent="-1"/>',
                    '<unit kind="litre"  exponent="1"/>',
                    '<unit kind="second" exponent="-1"/>',
                '</listOfUnits>',
            '</unitDefinition>',
        '</listOfUnitDefinitions>',
        '<listOfCompartments>',
			'<compartment id="reaction" size="1e-3" />',
		'</listOfCompartments>',
		'<listOfSpecies>']
	
	if(output_condensed):
		condensed = condense_resting_states(enumerator, **condense_options)
		complexes = condensed['resting_states']
		reactions = condensed['reactions']
	else:
		complexes = enumerator.complexes
		reactions = enumerator.reactions
	
	def id(species):
		return "s_"+species.name

	# build elements for each species
	if(output_condensed):
		for resting_state in complexes:
			# is_initial = any(c in enumerator.initial_complexes for c in resting_state.complexes)
			initial_concentration = sum((c.concentration if c.concentration is not None else 0.0) for c in resting_state.complexes)
			out.append('<species compartment="reaction" id="%(id)s" name="%(name)s" initialConcentration="%(initial)g"/>' \
				% {"name": resting_state.name, "id": id(resting_state), "initial": initial_concentration })
	else:
		for complex in complexes:
			# is_initial = (complex in enumerator.initial_complexes)
			initial_concentration = complex.concentration if complex.concentration is not None else 0.0
			out.append('<species compartment="reaction" id="%(id)s" name="%(name)s" initialConcentration="%(initial)g"/>' \
				% {"name": complex.name, "id": id(complex), "initial": initial_concentration })
	
	out += ['</listOfSpecies>','<listOfReactions>']

	# list reactions
	for (i, reaction) in enumerate(reactions):
		out += ['<reaction id="r_%d" reversible="false">' % i,
                '<listOfReactants>'] + \
					['<speciesReference species="%s"/>' % id(species) for species in reaction.reactants] + \
				['</listOfReactants>',
                '<listOfProducts>'] + \
					['<speciesReference species="%s"/>' % id(species) for species in reaction.products]	+ \
				['</listOfProducts>']

		# unimolecular rate constants have units 1/s, bimolecular rate 
		# constants have units 1/M/s
		units = 'per_second' if reaction.arity[0] == 1 else 'litre_per_mole_per_second'

		out += ['<kineticLaw>',
			'<math xmlns="http://www.w3.org/1998/Math/MathML">',
				'<apply>',
					'<times />',
					'<ci>k</ci>'] + \
					['<ci>'+id(s)+'</ci>' for s in reaction.reactants] + \
				['</apply>',
			'</math>',
			'<listOfParameters>',
	            '<parameter id="k"  value="%.10f" units="%s"/>' % (reaction.rate(), units), 
	        '</listOfParameters>',
		'</kineticLaw>',
		'</reaction>']

	out.extend(['</listOfReactions>','</model>','</sbml>']);

	
	doc = xml.dom.minidom.parseString("".join(out))
	
	fout = open(filename, "w")
	fout.write(header+'\n'+doc.documentElement.toprettyxml(indent="\t"))
	fout.close()
	

def output_crn(enumerator, filename, output_condensed = False, output_rates = True, condense_options = {}):
	output_file = open(filename, 'w')

	def write_reaction(output_file,reaction):
		reactants = map(str,reaction.reactants)
		products = map(str,reaction.products)
		reac_string_list = [" + ".join(reactants),"->"," + ".join(products),"\n"]
		reac_string = ' '.join(reac_string_list)
		output_file.write(reac_string)

	reactions = enumerator.reactions
	if (output_condensed):
		condensed = condense_resting_states(enumerator, **condense_options)
		reactions = condensed['reactions']

	for reaction in sorted(reactions): #utils.natural_sort(reactions):
		write_reaction(output_file,reaction)

	output_file.close()

def output_k(enumerator, filename, output_condensed = False, output_rates = True, condense_options = {}):
	output_file = open(filename, 'w')

	def write_reaction(output_file,reaction,reversible=True):
		reactants = sorted([ str(count)+'*'+str(x) if count > 1 else str(x) \
			for (x, count) in collections.Counter(sorted(reaction.reactants)).iteritems() ])
		products = sorted([ str(count)+'*'+str(x) if count > 1 else str(x) \
			for (x, count) in collections.Counter(sorted(reaction.products)).iteritems() ])

		if reversible: arrow = '<->'
		else: arrow = '->'

		rxn = [" + ".join(reactants),arrow," + ".join(products),"\n"]
		output_file.write(' '.join(rxn))

	reactions = enumerator.reactions
	if (output_condensed):
		condensed = condense_resting_states(enumerator, **condense_options)
		reactions = condensed['reactions']

	for reaction in sorted(reactions): #utils.natural_sort(reactions):
		write_reaction(output_file,reaction,True)

	output_file.close()

# warning: I tried to name this output_test_case, but nosetests kept trying to 
# run it as a test because magic. 
# https://nose.readthedocs.org/en/latest/writing_tests.html#test-modules
def output_case(enumerator, filename, output_condensed = False, condense_options = {}):
	def tab(x):
		return "\t" * x

	of = open(filename, 'w')
	of.write("def test_%s(enumerator):\n" % filename)

	# Domains
	of.write(tab(1) + "# Domains \n")
	of.write(tab(1) + "domains = { \n")
	lines = []
	for domain in utils.natural_sort(enumerator.domains):
		#  name, length, is_complement=False, sequence=None
		lines.append(tab(2) + "'%s' : Domain('%s', %d, is_complement=%s, sequence='%s')" \
			% (domain.name, domain.identity, domain.length, domain.is_complement, domain.sequence))

	of.write(",\n".join(lines) + "\n")
	of.write(tab(1) + "}\n")
	of.write(tab(1) + "assert set(domains.values()) == set(enumerator.domains)\n\n")

	# Strands
	of.write(tab(1) + "# Strands \n")
	of.write(tab(1) + "strands = { \n")
	lines = []
	for strand in utils.natural_sort(enumerator.strands):
		doms = ", ".join("domains['%s']" % dom.name for dom in strand.domains)
		lines.append(tab(2) + "'%s' : Strand('%s', [%s])" % (strand.name, strand.name, doms))
	of.write(",\n".join(lines) + "\n")
	of.write(tab(1) + "}\n")
	of.write(tab(1) + "assert set(strands.values()) == set(enumerator.strands)\n\n")

	# Complexes
	of.write(tab(1) + "# Complexes \n")
	of.write(tab(1) + "complexes = { \n")
	lines = []
	for complex in utils.natural_sort(enumerator.complexes):
		strands = ", ".join("strands['%s']" % strand.name for strand in complex.strands)
		lines.append(tab(2) + "'%s' : Complex('%s', [%s], %r)" % (complex.name, complex.name, strands, complex.structure))
	of.write(",\n".join(lines) + "\n")
	of.write(tab(1) + "}\n")
	of.write(tab(1) + "assert set(complexes.values()) == set(enumerator.complexes)\n\n")

	# Reactions
	of.write(tab(1) + "# Reactions \n")
	of.write(tab(1) + "reactions = { \n")
	lines = []
	for reaction in utils.natural_sort(enumerator.reactions):
		reactants = ", ".join("complexes['%s']" % complex.name for complex in reaction.reactants)
		products  = ", ".join("complexes['%s']" % complex.name for complex in reaction.products)
		lines.append(tab(2) + "ReactionPathway('%s', [%s], [%s])" % (reaction.name, reactants, products))
	of.write(",\n".join(lines) + "\n")
	of.write(tab(1) + "}\n")
	of.write(tab(1) + "assert set(reactions) == set(enumerator.reactions)\n\n")


	if(output_condensed):
		condensed = condense_resting_states(enumerator, **condense_options)
		
		# Resting states
		of.write(tab(1) + "# Resting states \n")
		of.write(tab(1) + "resting_states = { \n")
		lines = []
		for rs in utils.natural_sort(condensed['resting_states']):
			complexes = ", ".join("complexes['%s']" % complex.name for complex in rs.complexes)
			lines.append(tab(2) + "'%s' : RestingState('%s', [%s])" % (rs.name, rs.name, complexes))
	
		of.write(",\n".join(lines) + "\n")
		of.write(tab(1) + "}\n")
		of.write(tab(1) + "assert set(resting_states.values()) == set(condensed['resting_states'])\n\n")


		# Condensed reactions
		of.write(tab(1) + "# Condensed Reactions \n")
		of.write(tab(1) + "condensed_reactions = { \n")
		lines = []
		for reaction in utils.natural_sort(condensed['reactions']):
			reactants = ", ".join("resting_states['%s']" % species.name for species in reaction.reactants)
			products  = ", ".join("resting_states['%s']" % species.name for species in reaction.products)
			lines.append(tab(2) + "ReactionPathway('%s', [%s], [%s])" % (reaction.name, reactants, products))
		of.write(",\n".join(lines) + "\n")
		of.write(tab(1) + "}\n")
		of.write(tab(1) + "assert set(condensed_reactions) == set(condensed['reactions'])\n\n")

	of.close()





text_output_functions = {
	# 'standard': output_legacy,
	'legacy': output_legacy,
	'pil': output_pil,
	'kernel': output_kernel,
	'json': output_json,
	'enjs': output_json,
	'sbml': output_sbml,
	'crn': output_crn,
	'test': output_case
}

graph_output_functions = {
	'graph': output_graph
}
