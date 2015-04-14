#
#  input.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/26/10.
#

import logging
import json
import utils
from utils import *
from enumerator import Enumerator
from reactions import ReactionPathway, get_auto_name
import reactions
import re

def input_enum(filename):
	"""
	Initializes and returns an enumerator from a standard input file.
	"""
	fin = open(filename, 'r')
	
	domains = {}
	strands = {}
	complexes = {}
	
	line_counter = 1
	line = fin.readline()
	
	# We loop over all the lines in the file
	while (line != ""):
		line = line.strip()
		
		# This was an empty line
		if line == "":
			line = fin.readline()
			continue
			
		# This was a comment
		elif line.startswith("#"):
			line = fin.readline()
			continue
			
		# This is the beginning of a domain
		elif line.startswith("domain"):
			# e.g.: 
			#       "domain a : 6"
			# parts: 0      1 2 3
			
			parts = line.split()
						
			domain_name = parts[1]
			if domain_name in domains:
				logging.error("Duplicate domain name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not re.match(r'\w+$',domain_name):
				logging.warn("Non-alphanumeric domain name %s encountered in input line %d"
								% (domain_name, line_counter))
			
			# The domain length could be either short or long or it could be
			# an exact number
			domain_length = parts[3]
			if not ((domain_length == 'short') or (domain_length == 'long')):
				domain_length = int(domain_length)
				if domain_length <= 0:
					logging.warn("Domain of length %d found in input line %d" 
									% (domain_length, line_counter))
				
			# Check to see if a sequence is specified
			if len(parts) > 4:
				domain_sequence = parts[4]
			else:
				domain_sequence = None
			
			# Create the new domains
			new_dom = Domain(domain_name, domain_length, 
							 sequence=domain_sequence)
			new_dom_comp = Domain(domain_name, domain_length, 
								  sequence=domain_sequence,
								  is_complement=True)
			
			domains[domain_name] = new_dom
			domains["%s*" % domain_name] = new_dom_comp
		
		# This is the beginning of a strand	
		elif line.startswith("strand"):
			# e.g.: 
			#       "strand A : a x b y z* c* y* b* x*"
			# parts: 0      1 2 3 4 5 6 ...
			
			parts = line.split()
			
			strand_name = parts[1]
			if strand_name in strands:
				logging.error("Duplicate strand name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not re.match(r'\w+$',strand_name):
				logging.warn("Non-alphanumeric strand name %s encountered in input line %d"
								% (strand_name, line_counter))
			
			strand_doms = []
			for domain_name in parts[3:]:
				if not domain_name in domains:
					logging.error("Invalid domain name %s encountered in input line %d"
									% (domain_name, line_counter))
					raise Exception()
								
				strand_doms.append(domains[domain_name])
			
			if len(strand_doms) == 0:
				logging.warn("0-length strand encountered in input line %d")			
				
			new_strand = Strand(strand_name, strand_doms)
			strands[strand_name] = new_strand
			
		# This is the beginning of a complex	
		elif line.startswith("complex"):
			# e.g.:
			# complex A :\n
			# A\n				<- strands_line
			# .(((..)))\n		<- structure_line
			
			
			parts = line.split()
			
			complex_name = parts[1]
			if complex_name in complexes:
				logging.error("Duplicate complex name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not re.match(r'\w+$',complex_name):
				logging.warn("Non-alphanumeric complex name %s encountered in input line %d"
								% (complex_name, line_counter))
			
			complex_strands = []
			
			strands_line = fin.readline()
			strands_line = strands_line.strip()
			strands_line_parts = strands_line.split()
			for strand_name in strands_line_parts:
				if not strand_name in strands:
					logging.error("Invalid strand name %s encountered in input line %d"
									% (strand_name, line_counter))
					raise Exception()
				else:
					complex_strands.append(strands[strand_name])
				
			structure_line = fin.readline()
			structure_line = structure_line.strip()
			
			complex_structure = parse_dot_paren(structure_line)
			
			
			struct_length = sum(map(len,complex_structure))
			domains_length = sum(map(len,complex_strands))
			if(struct_length != domains_length):
				logging.error("Complex %(name)s has %(doms)d domains but structure size %(struct_length)d. (structure was %(struct)s)"
								% {"name":complex_name,"doms":domains_length,"struct_length":struct_length, "struct":structure_line})
				raise Exception()
			
			complex = Complex(complex_name, complex_strands, complex_structure)
			complex.check_structure()
			complexes[complex_name] = complex			
			
		else:
			logging.error("Unexpected characters encountered in input line %d"
							% line_counter)
			raise Exception()
		line = fin.readline()
		line_counter += 1
		
	domains = domains.values()
	strands = strands.values()
	complexes = complexes.values()
	
	enumerator = Enumerator(domains, strands, complexes)
	return enumerator


def parse_kernel(line):
	from pyparsing import Group, Forward, Word, Combine, Literal, Optional, Suppress, ZeroOrMore, OneOrMore, StringEnd, delimitedList, nestedExpr, alphanums
	# letter = (* any alphanumeric character *)
	# identifier = letter, {letter}, [*]
	# pattern = expression, {space, expression}
	# expression = domain | loop | wildcard | break
	# domain = identifier
	# loop = identifier, "(", [pattern] ,")"
	# wildcard = "?" | "_" | "~" (* ... etc.*)
	# break = "+"

	identifier = Word(alphanums+"_-").setName("identifier")
	domain = Combine(identifier + Optional(Literal("^")) + Optional(Literal("*"))).setName("domain")
	sbreak = Literal("+").setName("strand break")
	wildcard = Literal("?").setName("wildcard")
	
	pattern = Forward()
	expression = Forward()
	
	loop = (domain + Suppress("(") + Group(Optional(pattern)) + Suppress(")")).setName("loop")
	expression <<  (loop | wildcard | sbreak | domain)
	pattern << OneOrMore(expression)

	rule = pattern + StringEnd()

	return rule.parseString(line)


def parse_identifier(identifier):
	"""
	Parse an identifier (e.g. `a^*`) to figure out the name, 
	polarity (+1 or -1) and the length ('long' or 'short').

	Returns: (name, polarity, length)
	"""

	polarity = 1
	length = 'long'
	if identifier[-1] == '*':
		polarity = -1
		identifier = identifier[:-1]
	if identifier[-1] == '^':
		identifier = identifier[:-1]
		length = 'short'
	return (identifier, polarity, length)

def auto_domain(name, polarity, domains):
	"""
	Finds or automatically generates a domain and/or its complement.
	"""

	# figure name, polarity, length
	identity, plrt, length = parse_identifier(name)
	polarity *= plrt

	# handle complements
	identifier = identity + ('*' if polarity == -1 else '')

	# search for existing domain
	if identifier in domains: 
		dom = domains[identifier]
		# if (length == 'short') != (dom.length == utils.SHORT_DOMAIN_LENGTH):
		# 	error(("Domain '%s' should be %s, but there is a already a domain %s that is not. "+ \
		# 		"I assume this was a mistake; please give all domains '%s' a caret (^), "+ \
		# 		"or remove the caret from all domain") % (name, length, dom.name, dom.name))
		return dom

	# generate new domain
	else:
		if identity in domains:
			raise "Error!"
		if "*" in identity:
			raise "Error!"

		sequence = 'N' * utils.resolve_length(length)

		new_dom = Domain(identity, length, sequence=sequence)
		new_dom_comp = Domain(identity, length,
							  is_complement=True, sequence=sequence)
		domains[identity] = new_dom
		domains[identity+'*'] = new_dom_comp

		# return domain or complement depending on requested polarity
		return new_dom if polarity == 1 else new_dom_comp

def auto_strand(doms, strands, structures_to_strands):
	"""
	Finds or automatically generates a strand from a list of 
	domains.
	"""

	# look up whether any strand with this structure exists
	tdoms = tuple(doms)
	if tdoms in structures_to_strands:
		return structures_to_strands[tdoms]
	
	# if not, make up a name
	else:
		auto_name = initial_auto_name = "_".join(d.name for d in doms)
		
		# if another strand exists with this name but not this 
		# structure, generate a new, uglier name 
		if auto_name in strands:
			auto_name += "_%d" 
			index = 2
			while (auto_name % index) in strands:
				index += 1
			auto_name = auto_name % index

			# TODO: warn about this
			warning("Auto-generated strand name %s already taken by a strand with different structure. Auto-generated strand will be named %s." % (initial_auto_name, auto_name))

		# generate new strand object
		strand = Strand(auto_name, doms)
		strands[auto_name] = strand
		structures_to_strands[tuple(doms)] = strand
		return strand

def auto_complex(strands, structure, complexes, name=None):
	"""
	Generates a complex from a list of strands and a structure
	"""
	if name is None:
		name = get_auto_name()
	complex = Complex(name, strands, structure)
	complexes[complex.name] = complex
	return complex

def resolve_kernel(lines, domains, strands, structures_to_strands, complexes):
	def resolve(parts, strands, structure):
		"""
		Recursively parse a list of `parts` generated by parse_kernel,
		populating the `strands` and `structure` list
		"""
		for part in parts:

			# strand break
			if part == '+':
				strands.append([])
				structure.append([])
			# ignore
			elif part == '(':
				pass
			elif part == ')':
				pass

			# unpaired domain
			elif isinstance(part, basestring):
				strands[-1].append(auto_domain(part, 1, domains))
				structure[-1].append(None)

			# paired domain
			else:
				# remember opening domain (last domain of last strand)
				last_dom = strands[-1][-1]
				last_index = (len(strands)-1, len(strands[-1])-1 )
				
				# resolve loop
				resolve(part, strands, structure)

				# add closing domain
				strands[-1].append( auto_domain(str(last_dom), -1, domains) )
				structure[-1].append( last_index )
				structure[last_index[0]][last_index[1]] = (len(strands)-1, len(strands[-1])-1 )
				# I was sad that this didn't work without numpy...
				# structure[last_index] = ...

	for line in lines:
		cname = None
		parameters = None
		parts = re.match(r"(\w+)\s*(\[[^\]]+\])?\s*=\s*(.*)", line)
		if parts is not None:
			cname, parameters, line = parts.groups()

		# parse parameters
		if parameters is None:
			parameters = ""
		params = utils.parse_parameters(parameters)

		# parse the structure portion
		kparts = parse_kernel(line)

		stack = []
		kstrands = [[]]
		kstructure = [[]]

		resolve(kparts, kstrands, kstructure)

		# build strands
		kstrands = [auto_strand(doms, strands, structures_to_strands) for doms in kstrands]

		# build complex
		complex = auto_complex(kstrands, kstructure, complexes, name=cname)

		# check structure is valid
		complex.check_structure()

		# apply parameters
		if params['concentration'] is not None:
			complex.concentration = params['concentration']

def from_kernel(lines):
	# split string into lines if necessary
	if isinstance(lines, basestring):
		lines = lines.split("\n")

	# remove blank lines
	lines = filter(None, lines)

	domains = {}
	strands = {}
	structures_to_strands = {}
	complexes = {}
	resolve_kernel(lines, domains, strands, structures_to_strands, complexes)
	return (domains, strands, complexes)

def enum_from_kernel(lines):
	(domains, strands, complexes) = from_kernel(lines)
	return Enumerator(domains, strands, complexes)

def input_pil(filename):
	"""
	Initializes and returns an enumerator from an input file in the Pepper Intermediate Language (PIL)
	"""
	fin = open(filename, 'r')
	
	domains = {}
	strands = {}
	complexes = {}

	# maps domain-wise strand structures to auto-generated strand names
	structures_to_strands = {}

	# We loop over all the lines in the file
	for (line_counter, line) in enumerate(fin,start=1):
		line = line.strip()
		
		# This was an empty line
		if line == "":
			continue
			
		# This was a comment
		elif line.startswith("#"):
			continue
		
		elif line.startswith("length"):
			# e.g.:
			#       "length a = 6"
			# parts:        0   1 
			parts = re.match(r"length\s*([\w-]+)\s*=\s*(\d+)\s*",line)
			if parts == None:
				logging.error("Invalid syntax on input line %d"
						% line_counter)
				logging.error(line)
				raise Exception()

			domain_name, domain_length = parts.groups()
			if domain_name in domains:
				logging.error("Duplicate domain name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not re.match(r'[\w-]+$',domain_name):
				logging.warn("Non-alphanumeric domain name %s encountered in input line %d"
								% (domain_name, line_counter))

			domain_length = int(domain_length)
			domain_sequence = "N" * domain_length

			# Create the new domains
			new_dom = Domain(domain_name, domain_length, 
							 sequence=domain_sequence)
			new_dom_comp = Domain(domain_name, domain_length, 
								  sequence=domain_sequence,
								  is_complement=True)
			
			domains[domain_name] = new_dom
			domains["%s*" % domain_name] = new_dom_comp

		# This is the beginning of a domain
		elif line.startswith("sequence"):
			# e.g.: 
			#       "sequence a = 6 : 6"
			# parts: 0        1 2 3
			
			#                  sequence    a         =   NNNNN   :     6
			parts = re.match(r"sequence\s*([\w-]+)\s*=\s*(\w+)\s*:?\s*(\d?)\s*",line)
			if parts == None:
				logging.error("Invalid syntax on input line %d"
							% line_counter)
				logging.error(line)
				raise Exception()

			domain_name, domain_sequence, length = parts.groups()
			if domain_name in domains:
				logging.error("Duplicate domain name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not re.match(r'[\w-]+$',domain_name):
				logging.warn("Non-alphanumeric domain name %s encountered in input line %d"
								% (domain_name, line_counter))
			
			# The sequence specification
			# domain_sequence = parts[1]
			domain_length = len(domain_sequence);
			
			# Create the new domains
			new_dom = Domain(domain_name, domain_length, 
							 sequence=domain_sequence)
			new_dom_comp = Domain(domain_name, domain_length, 
								  sequence=domain_sequence,
								  is_complement=True)
			
			domains[domain_name] = new_dom
			domains["%s*" % domain_name] = new_dom_comp
		
		elif line.startswith("sup-sequence"):
			# e.g.
			#       "sequence a = b c d e : 6"
			#                 0   1         2
			parts = re.match(r"sup-sequence\s*([\w-]+)\s*=\s*((?:[\w-]+\s*)+):?(\d?)",line)
			if parts == None:
				logging.error("Invalid syntax on input line %d"
							% line_counter)
				logging.error(line)
				raise Exception()

			domain_name, sequence_names, length = parts.groups()

			# domain name
			if domain_name in domains:
				logging.error("Duplicate domain name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not re.match(r'[\w-]+$',domain_name):
				logging.warn("Non-alphanumeric domain name %s encountered in input line %d"
								% (domain_name, line_counter))


			# subsequences
			sequence = ""
			for sequence_name in sequence_names.split():
				sequence_name = sequence_name.strip()
				if sequence_name == "": continue
				
				# make sure each subsequence is defined
				if not sequence_name in domains:
					logging.error("Unknown domain name '%s' in super-sequence on input line %d"
								% (sequence_name,line_counter))
					logging.error(line)
					raise Exception()

				# build up the full sequence
				sequence += domains[sequence_name].sequence

			# check for correctness
			if length:
				if int(length) != len(sequence):
					logging.error("Sequence length for super-sequence %s is %d, not equal to expected value %d on input line %d"
								% (domain_name, len(sequence), length, line_counter))
					raise Exception()

			# The sequence specification
			domain_sequence = sequence
			domain_length = len(sequence)

			# Create the new domains
			new_dom = Domain(domain_name, domain_length, 
							 sequence=domain_sequence)
			new_dom_comp = Domain(domain_name, domain_length, 
								  sequence=domain_sequence,
								  is_complement=True)
			
			domains[domain_name] = new_dom
			domains["%s*" % domain_name] = new_dom_comp

		elif line.startswith("equal"):
			parts = line.split()
			if len(parts) < 3:
				logging.error("'equal' statement does not specify at least 2 domains on input line %d"
								% line_counter)
				logging.error(line)
				raise Exception()

			source_domain_name = parts[1]
			target_domain_names = parts[2:]

			if not source_domain_name in domains:
				logging.error("Unknown domain name '%s' in 'equals' statement on input line %d"
								% (source_domain_name,line_counter))
				logging.error(line)
				raise Exception()

			source_domain = domains[source_domain_name]


			for target_domain_name in target_domain_names:
				new_dom = Domain(target_domain_name, len(source_domain), 
							 sequence=source_domain.sequence)
				new_dom_comp = Domain(target_domain_name, len(source_domain), 
									  sequence=source_domain.sequence,
									  is_complement=True)

				domains[target_domain_name] = new_dom
				domains["%s*" % target_domain_name] = new_dom_comp


		# This is the beginning of a strand	
		elif line.startswith("strand"):
			# e.g.: 
			#       "strand A = a x b y z* c* y* b* x*"
			# parts: 0      1 2 3 4 5 6 ...
			
			parts = re.match(r"strand\s*([\w-]+)\s*=\s*((?:[\w*-]+\s*)+):?(\d?)",line)
			if parts == None:
				logging.error("Invalid syntax on input line %d"
							% line_counter)
				logging.error(line)
				raise Exception()
			
			strand_name, strand_dom_names, length = parts.groups()

			if strand_name in strands:
				logging.error("Duplicate strand name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not re.match(r'\w+$',strand_name):
				logging.warn("Non-alphanumeric strand name %s encountered in input line %d"
								% (strand_name, line_counter))
			
			strand_doms = []
			for domain_name in filter(None,strand_dom_names.split()):
				if not domain_name in domains:
					logging.error("Invalid domain name %s encountered in input line %d"
									% (domain_name, line_counter))
					logging.error(line)
					print domains
					raise Exception()
								
				strand_doms.append(domains[domain_name])
			
			if len(strand_doms) == 0:
				logging.warn("0-length strand encountered in input line %d")			
				
			new_strand = Strand(strand_name, strand_doms)
			strands[strand_name] = new_strand
			
		# This is the beginning of a complex	
		elif line.startswith("structure"):
			# parse `structure` line:
			# e.g.:
			# structure A = S1 : .(((..)))
			#                  structure    [  1nt  ]      name      =   s1 s2 s3 + s4         : ....((+))...((..))....
			parts = re.match(r"structure\s+(\[[^\]]+\])?\s*([\w-]+)\s*=\s*((?:[\w-]+\s*\+?\s*)+):\s*([().+\s]+)",line)

			if parts == None:

				# parse `structure` line:
				# e.g.:
				# structure A = S1 : .(((..)))
				#                  structure    name      =   s1 s2 s3 + s4         : ....((+))...((..))....
				parts = re.match(r"structure\s+([\w-]+)\s*=\s*((?:[\w-]+\s*\+?\s*)+):\s*([().+\s]+)",line)

				if parts == None:

					logging.error("Invalid syntax on input line %d"
								% line_counter)
					logging.error(line)
					raise Exception()
				else:
					complex_name, strands_line, structure_line = parts.groups()
					parameters = ""
			else:
				parameters, complex_name, strands_line, structure_line = parts.groups() 
			
			# parse parameters
			if parameters is None:
				parameters = ""
			params = utils.parse_parameters(parameters)

			# check for duplicate complex name
			if complex_name in complexes:
				logging.error("Duplicate complex name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			# check for non-alphanumeric complex name
			if not re.match(r'\w+$',complex_name):
				logging.warn("Non-alphanumeric complex name %s encountered in input line %d"
								% (complex_name, line_counter))
			
			# get strand names, allowing optional '+' characters
			complex_strands = []
			strands_line_parts = [name for name in strands_line.split() if name != "+"]
			for strand_name in strands_line_parts:
				if not strand_name in strands:
					logging.error("Invalid strand name %s encountered in input line %d"
									% (strand_name, line_counter))
					raise Exception()
				else:
					complex_strands.append(strands[strand_name])

			# parse dot-paren structure, then do some horrible magic to guess 
			# if it's basewise or segmentwise... 
			complex_structure = parse_dot_paren(structure_line)
			struct_length = sum(map(len,complex_structure))
			domains_length = sum(map(len,complex_strands)) #sum([ len(d) for c in complex_strands for d in c.domains ])
						
			if(struct_length > domains_length):
				complex_structure = parse_basewise_dot_paren(structure_line, complex_strands)						
				struct_length = sum(map(len,complex_structure))
				if(struct_length != domains_length):
					logging.error("Complex %(name)s has %(doms)d domains but structure size %(struct_length)d. (structure was '%(struct)s')"
								% {"name":complex_name,"doms":domains_length,"struct_length":struct_length, "struct":structure_line})
					raise Exception()
			
			elif(struct_length != domains_length):

				logging.error("Complex %(name)s has %(doms)d domains but structure size %(struct_length)d. (structure was '%(struct)s')"
								% {"name":complex_name,"doms":domains_length,"struct_length":struct_length, "struct":structure_line})
				raise Exception()
			
			complex = Complex(complex_name, complex_strands, complex_structure)
			complex.check_structure()

			# apply parameters
			if params['concentration'] is not None:
				complex.concentration = params['concentration']

			complexes[complex_name] = complex		
		elif line.startswith("kinetic"):
			continue	
		elif line.strip() == "":	
			continue
		else:
			try:
				resolve_kernel([line], domains, strands, structures_to_strands, complexes)
			except Exception as e:
				logging.error("Unexpected characters encountered in input line %d; tried to parse as Kernel statement but got error: %s"
								% (line_counter, str(e)))
				raise Exception()
			
		# line = fin.readline()
		# line_counter += 1
		
	fin.close()
	domains = domains.values()
	strands = strands.values()
	complexes = complexes.values()
	
	enumerator = Enumerator(domains, strands, complexes)
	return enumerator

			
def load_json(filename):
	"""
	Loads a saved enumerator from a JSON output file at filename.
	"""
	
	fin = open(filename, 'r')
	saved = json.load(fin)
	saved_domains = saved['domains']
	
	domains = {}
	for saved_domain in saved_domains:
		if not 'sequence' in saved_domain:
			saved_domain['sequence'] = None
		if (saved_domain['is_complement']):
			saved_domain['name'] = saved_domain['name'][:-1]
		new_dom = Domain(saved_domain['name'], saved_domain['length'], is_complement=saved_domain['is_complement'], sequence=saved_domain['sequence'])
		domains[new_dom.name] = new_dom
	
	saved_strands = saved['strands']
	
	strands = {}
	for saved_strand in saved_strands:
		doms = []
		for domain in saved_strand['domains']:
			doms.append(domains[domain])
		
		new_strand = Strand(saved_strand['name'], doms)
		strands[saved_strand['name']] = new_strand
	
	complexes = {}
	resting_complexes = {}
	
	saved_resting_complexes = saved['resting_complexes']
	for saved_complex in saved_resting_complexes:
		c_strands = []
		for strand in saved_complex['strands']:
			c_strands.append(strands[strand])
		new_structure = []
		for strand in saved_complex['structure']:
			new_strand = []
			for tup in strand:
				if (tup == None):
					new_strand.append(None)
				else:
					new_strand.append(tuple(tup))
			new_structure.append(new_strand)
		new_complex = Complex(saved_complex['name'], c_strands, new_structure)
		resting_complexes[saved_complex['name']] = new_complex
		complexes[saved_complex['name']] = new_complex
	
	transient_complexes = {}	
	
	saved_transient_complexes = saved['transient_complexes']
	for saved_complex in saved_transient_complexes:
		c_strands = []
		for strand in saved_complex['strands']:
			c_strands.append(strands[strand])
		new_structure = []
		for strand in saved_complex['structure']:
			new_strand = []
			for tup in strand:
				if (tup == None):
					new_strand.append(None)
				else:
					new_strand.append(tuple(tup))
			new_structure.append(new_strand)
		new_complex = Complex(saved_complex['name'], c_strands, new_structure)
		transient_complexes[saved_complex['name']] = new_complex
		complexes[saved_complex['name']] = new_complex
	
	saved_reactions = saved['reactions']
	
	reactions = []
	
	for saved_reaction in saved_reactions:
		reactants = []
		for reactant in saved_reaction['reactants']:
			reactants.append(complexes[reactant])
		products = []
		for product in saved_reaction['products']:
			products.append(complexes[product])
			
		reaction = ReactionPathway(saved_reaction['name'], reactants, products)
		reactions.append(reaction)
	
	resting_states = []
	
	for resting_state in saved['resting_states']:
		comps = []
		for complex in resting_state['complexes']:
			comps.append(complexes[complex])
		resting_states.append(RestingState(resting_state['name'], comps))

	initial_complexes = {}
	for saved_complex in saved['initial_complexes']:
		c_strands = []
		for strand in saved_complex['strands']:
			c_strands.append(strands[strand])
		new_structure = []
		for strand in saved_complex['structure']:
			new_strand = []
			for tup in strand:
				if (tup == None):
					new_strand.append(None)
				else:
					new_strand.append(tuple(tup))
			new_structure.append(new_strand)
		new_complex = Complex(saved_complex['name'], c_strands, new_structure)
		initial_complexes[saved_complex['name']] = new_complex
		
			
	enumerator = Enumerator(domains.values(), strands.values(), initial_complexes.values())
	enumerator._complexes = complexes.values()
	enumerator._resting_states = resting_states
	enumerator._transient_complexes = transient_complexes.values()
	enumerator._resting_complexes = resting_complexes.values()
	enumerator._reactions = reactions

	
	return enumerator



text_input_functions = {
						'enum': input_enum,
						'pil': input_pil
					  }
					  
load_input_functions = {
						 'json': load_json
					   }