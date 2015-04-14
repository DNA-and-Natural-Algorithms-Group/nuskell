#
#  reactions.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/18/2010.
#  Modifications by Casey Grun and Erik Winfree 8/15/2014.

import copy
import utils
from utils import *

auto_name = 0
def get_auto_name():
	"""
	Returns a new unique name
	"""
	global auto_name
	auto_name += 1
	return str(auto_name)

RELEASE_CUTOFF_1_1 = 6
RELEASE_CUTOFF_1_N = 6
"""
Gives the maximum length of a duplex in nucleotides that should be considered 
reversibly bound; that is, helices longer than RELEASE_CUTOFF will never be
unbound by the `open` reaction function.
"""

REJECT_REMOTE = False
"""
If True, discards 3-way and 4-way remote toehold branch migration reactions
"""


# If true, 3 way branch migrations are always greedy
UNZIP = True
"""
If True, 3-way branch migrations obey "Maximum helix at a time" semantics
"""
# UNZIP = False
LEGACY_UNZIP = True
# LEGACY_UNZIP = False

class ReactionPathway(object):
	"""
	Represents a reaction node on a reaction graph. Has a list of reactants
	and a list of products, along with a name field to denote the type of
	reaction involved.
	"""
	
	def __init__(self, name, reactants, products):
		"""
		Default constructor. Takes a name, a list of reactants (Complexes) and
		a list of products (Complexes).
		"""
		self._name = name
		self._reactants = reactants
		self._reactants.sort()
		self._products = products
		self._products.sort()
		self._const = 1

		for x in self._reactants:
			if isinstance(x,Complex):
				assert x.check_structure()
				assert x.check_pseudoknots()
				assert x.check_connected()

		for x in self._products:
			if isinstance(x,Complex):
				assert x.check_structure() 
				assert x.check_pseudoknots()
				assert x.check_connected()
		
	def __repr__(self):
		return self.full_string()

	def __str__(self):
		return self.name

	def full_string(self):
		# return "ReactionPathway(%s): %s -> %s" % (self.name, str(self.reactants), str(self.products))
		return "ReactionPathway(\"%s\",%s,%s)" % (self.name, str(self.reactants), str(self.products))

	@property
	def name(self):
		"""
		Gives the name of the move type (reaction function) that generated the ReactionPathway
		"""
		return self._name
		
	@property
	def reactants(self):
		"""
		Gives the list of reactants of the ReactionPathway
		"""
		return self._reactants
		
	@property
	def products(self):
		"""
		Gives the list of products of the ReactionPathway
		"""
		return self._products
		
	@property
	def arity(self):
		"""
		Gives a pair containing the number of reactants and the number of products in the reaction
		"""
		return (len(self._reactants),len(self._products))	
	
	def __eq__(self, other):
		"""
		Compares two ReactionPathway objects. ReactionPathway objects are equal if they have the same name,
		reactants, and products.
		"""
		return (self.name == other.name) and \
			   (self.reactants == other.reactants) and \
			   (self.products == other.products)

	def __hash__(self):
		return hash(self.name) + hash(frozenset(self.reactants)) + hash(frozenset(self.products))

	def __cmp__(self, other):
		"""
		ReactionPathway objects are sorted by name, then by reactants, then by products.
		"""
		out = cmp(self.name, other.name)
		if (out != 0):
			return out
			
		out = cmp(self.reactants, other.reactants)
		if (out != 0):
			return out
			
		return cmp(self.products, other.products)
			   
	def normalize(self):
		"""
		Ensures that complexes appear on only one side of the reaction by
		removing them evenly from both sides until only one side has any.
		"""
		for reactant in self.reactants:
			while (reactant in self.reactants and \
				   reactant in self.products):
				self.reactants.remove(reactant)
				self.products.remove(reactant)
	
	def rate(self):
		"""
		Gives the rate constant for this reaction
		"""
		return self._const	

	def simple_string(self):
		return \
		"  +  ".join(str(r) for r in self.reactants) + \
		"  ->  " + \
		"  +  ".join(str(p) for p in self.products) 		

	def kernel_string(self):
		return \
		"  +  ".join(r.kernel_string() for r in self.reactants) + \
		"  ->  " + \
		"  +  ".join(p.kernel_string() for p in self.products) 

# def unpack_loop(loop):
# 	bases = 0
# 	stems = 0
# 	is_open = False
# 	strand = None
# 	for (dom, struct, loc) in loop:
# 		if struct is None:
# 			bases += 1

# 			if strand is not None and loc[1] != strand:
# 				is_open = True

# 			strand = loc[1]

# 		elif struct is not None:
# 			stems += 1

# 			if strand is not None and loc[1] != strand:
# 				is_open = True

# 			strand = struct[0]


				
# Rate constant formulae
# ----------------------------------------------------------------------------


def opening_rate(length):
	"""
	Rate constant formula for opening a duplex of a given `length`.
	"""
	# use k_open = k_hybrid exp( (length * dG_bp + dG_assoc) / RT )
	# where k_hybrid = 3x10^6 /M/s   from Zhang&Winfree 2009 and Srinivas et al 2013
	#       dG_bp = -1.7 kcal/mol
	#       dG_assoc = +1.9 kcal/mol
	#       R = 0.001987 kcal/mol/K
	#       T = (273.15 + 25) K
	# return 7.41e7 * (0.0567 ** length)
	#
	# instead, use k_hybrid = L * 3 * 10^5, which matches the above for L=10.
	# this is to be consistent with the bimolecular binding rate.
	return length * 7.41e6 * (0.0567 ** length)


def polymer_link_length(before, after):
	"""Effective length estimate for (ss+ds) linkers between two domains, one or both of which may be open."""
	L_stem = 2.0/0.43  # rough equivalent number of single-stranded nucleotides to span a stem
	if not before.is_open:
		L_before = 1 + before.bases + L_stem * before.stems
	if not after.is_open:
		L_after = 1 + after.bases + L_stem * after.stems
	if not after.is_open and not before.is_open:   # for both closed & open cases, assume shorter length matters most
		return min(L_before,L_after)
	if not before.is_open:                         
		return L_before    
	if not after.is_open:                          
		return L_after
	assert False, "should not have reached this case -- how can both sides be open?"
	raw_input("bad bad bad -- computing polymer lengths in disconnected complex!")
	
def polymer_link_rate(linker_length):
	"""Unimolecular hairpin closing rate, as a function of effective linker length. """
	# a = 2.5e7    # per second; fit from data in Bonnet et al 1998 only
	a = 1e6        # per second; Kuznetsov et al 2008, Nayak et al 2012, Tsukanov et al 2013ab, all say at least 10x slower
	b = -2.5                   # fit from data in Bonnet et al 1998, consistent with Kuznetsov et al 2008
	k = a * linker_length**b   # hairpin closing adapted (simplified) from data in Bonnet et al 1998, modified as above
	return min(k,33000)        # hairpins shorter than about 4 nt can't close any faster.

def binding_rate(length, before, after):
	"""
	Rate constant formula for unimolecular binding of a domain of the given length.  
	Could be zippering, hairpin closing, bubble closing, bulge formation, multiloop formation,
	depending on the flanking loops, which may be open or closed.  
	"""
	if not before.is_open and before.stems==1 and before.bases==0:
		if not after.is_open and after.stems==1 and after.bases==0:
			return 1e4  # bubble closing rate from Altan-Bonnet 2003
		return (1e6)/length  # zippering from Wetmur&Davidson 1968, Gueron&Leroy 1995, Srinivas et al 2013, low end
	if not after.is_open and after.stems==1 and after.bases==0:
		return (1e6)/length  # zippering from Wetmur&Davidson 1968, Gueron&Leroy 1995, Srinivas et al 2013, low end

	L = polymer_link_length(before, after)  # bulge closing assumed to be similar to faster of two hairpin closings
	return polymer_link_rate(L)             # hairpin closing adapted (simplified) from data in Bonnet et al 1998

# Diagram for 3-way branch migration, general case.
# Loops could be listed either 5'->3' or 3'->5, but they always go from invading domain to bound stem (non-inclusive).
#
#                before
#                _______  x (bound domain)
#               /       \____
# (invading) x |         ____
#               \_______/ x*
#                 
#                 after

def show_loops(before, after, message):
	"""
	Debugging help: show (partial) loops returned by find_on_loop().  ! indicates a stem, | indicates open loop break.
	"""
	print "before: [ ",
	for step in before.parts:
		print " | " if step==None else step[0].name+("!" if step[1] is not None else "")+" ",
	print " ] is_open = %r, stems = %d, bases = %d" % (before.is_open, before.stems, before.bases)
	print "after: [ ",
	for step in after.parts:
		print " | " if step==None else step[0].name+("!" if step[1] is not None else "")+" ",
	print " ] is_open = %r, stems = %d, bases = %d" % (after.is_open, after.stems, after.bases)
	raw_input(message)


def branch_3way_remote_rate(length, before, after):
	"""
	Rate constant formula for 3-way branch migration, possibly with a remote toehold
	"""
	# step = 0.1e-3  # 100us branch migration step time from Srinivas et al 2013 (not relevant)
	# k_init = k_uni exp(-dGsp / RT) with k_uni = 7.5e7, dGsp = 7.3 kcal/mol, T = 273.15 + 25 K, R = 0.001987 kcal/mol/K
	init = 3.0e-3  # sec, = 1/k_init from Srinivas et al 2013

	# show_loops(before, after, "...before & after loops for 3-way branch migration...")

	if not after.is_open and after.stems==1 and after.bases==0: 	# "standard" 3-way bm initiation (plus "before" being closed)
		return 1.0 / init / length  # each initiation has probability 1/length of succeeding.  how long it takes doesn't matter.

	# show_loops(before, after, "run_tests.py should not have any remote toeholds for 3-way branch migration")
	
	# consider a slowdown analogous to Genot et al 2011 (remote) supp info derivation
	L = polymer_link_length(before, after)  # bulge closing assumed to be similar to faster of two hairpin closings
	ratio = polymer_link_rate(L) / (1e6)    # how much slower than our (admittedly slow) zippering rate is this?
	return ratio / init / length            # we slow down initiation and thus success probability (note: ratio < 1/30)


def branch_4way_remote_rate(length, before, after):
	"""
	Rate constant formula for 4-way branch migration, possibly with a remote toehold
	"""
	# rates recalculated from Nadine Dabby, PhD Thesis, 2013, based on assuming last 6 bp dissociate faster than 4way bm step
	open_step   = 107   # sec, = 1/k_first  (this is for open 4-way case only)
	closed_step = 3.0   # sec, = 1/k_bm     (this is used for initiating closed 4-way case; consistent with Panyutin&Hsieh 1993)

	# open_step = 200 # fudge !

	# show_loops(before, after, "before & after loops for 4-way branch migration")

	if not before.is_open and not after.is_open:
		init = closed_step
		if before.bases==0 and before.stems==1 and after.bases==0 and after.stems==1:
			return 1.0 / init / length   # closed & ready-to-rock-and-roll 4 way initiation
	if before.is_open:
		init = open_step
		if after.bases==0 and after.stems==1:
			return 1.0 / init / length   # we care about probability of taking this path, not actual time
	if after.is_open:
		init = open_step
		if before.bases==0 and before.stems==1:
			return 1.0 / init / length   # we care about probability of taking this path, not actual time

	# show_loops(before, after, "run_tests.py should not have any remote toeholds for 4-way branch migration")

	# consider a slowdown analogous to Genot et al 2011 (remote) supp info derivation
	L = polymer_link_length(before, after)  # bulge closing assumed to be similar to faster of two hairpin closings
	ratio = polymer_link_rate(L) / (1e6)    # how much slower than our (admittedly slow) zippering rate is this?
	return ratio / init / length            # we slow down initiation and thus success probability (note: ratio < 1/30)

def bimolecular_binding_rate(length):
	"""
	Rate constant formula for bimolecular association (binding).
	"""
	# use k_hybrid = 3x10^6 /M/s   from Zhang&Winfree 2009
	# return 3.0e6
	#
	# instead, use k_hybrid = L * 3 * 10^5, which matches the above for L=10.
	# see Wetmur 1976 review, and Srinivas et al 2013 AEL model.
	# another motivation is to have binding rate approx = if a domain is divided into two domains.
	return length * 3e5


# Reaction functions
# ----------------------------------------------------------------------------

def bind11(reactant):
	"""
	Returns a list of reaction pathways which can be produced by 1-1 binding 
	reactions of the argument complex. The 1-1 binding reaction is the 
	hybridization of two complementary unpaired domains within a single complex 
	to produce a single unpseudoknotted product complex.
	"""
	reactions = []
	structure = reactant.structure
	
	# We iterate through all the domains
	for (strand_index, strand) in enumerate(reactant.strands):
		for (domain_index, domain) in enumerate(strand.domains):

			# The displacing domain must be free			
			if (structure[strand_index][domain_index] != None):
				continue

			d1 = strand.domains[domain_index]
			loc1 = (strand_index, domain_index)

			# search both directions around the loop for a bound domain that
			# has the same sequence (and therefore can be displaced)
			locs = find_on_loop(reactant, loc1, -1, \
				# lambda dom, struct, loc: struct is None and dom.can_pair(d1)) + \
				lambda dom1, struct1, loc1, dom2, struct2, loc2: struct1 is None and struct2 is None and dom2.can_pair(dom1)) + \
			find_on_loop(reactant, loc1, +1, \
				lambda dom1, struct1, loc1, dom2, struct2, loc2: struct1 is None and struct2 is None and dom2.can_pair(dom1))

			# build products
			for (loc1s, loc2s, before, after) in locs:
				product = do_bind11(reactant, loc1s.locs, loc2s.locs)
				reaction = ReactionPathway('bind11', [reactant], [product])

				# length of invading domain
				length = len(loc1s)

				# calculate reaction constant
				reaction._const = binding_rate(length, before, after)

				reactions.append(reaction)

			# for (loc2, before, after) in locs:
			# 	product = do_bind11(reactant, loc1, loc2)		
			# 	reaction = ReactionPathway('bind11', [reactant], [product])

			# 	# length of invading domain
			# 	length = len(d1)

			# 	# calculate reaction constant
			# 	reaction._const = binding_rate(length, before, after)
			# 	reactions.append(reaction)

	output = reactions

	# remove any duplicate reactions
	output = sorted(list(set(output)))
	return output

def do_single_bind11(reactant, loc1, loc2):
	new_structure = copy.deepcopy(reactant.structure)
	new_structure[loc1[0]][loc1[1]] = loc2
	new_structure[loc2[0]][loc2[1]] = loc1
	product = Complex(get_auto_name(),  #reactant.name + "(" + str(loc1) + "+" + str(loc2) + ")",
					  reactant.strands, new_structure)
	return product

def do_bind11(reactant, loc1s, loc2s):
	product = reactant
	for loc1, loc2 in zip(loc1s,loc2s):
		product = do_single_bind11(product, loc1, loc2)
	return product

# def do_bind11(reactant, loc1, loc2):
# 	new_structure = copy.deepcopy(reactant.structure)
# 	new_structure[loc1[0]][loc1[1]] = loc2
# 	new_structure[loc2[0]][loc2[1]] = loc1
# 	product = Complex(get_auto_name(),  #reactant.name + "(" + str(loc1) + "+" + str(loc2) + ")",
# 					  reactant.strands, new_structure)
# 	return product

def bind21(reactant1, reactant2):
	"""
	Returns a list of reaction pathways which can be produced by 2-1 binding 
	reactions of the argument complexes. The 2-1 binding reaction is the 
	hybridization of two complementary unpaired domains, each in a different 
	complex, to produce a single, unpseudoknotted product complex containing 
	all of the strands contained in either of the original complexes.
	"""
	r1_doms = reactant1.available_domains
	r2_doms = reactant2.available_domains
	
	interaction_list = []
	
	# TODO: significant optimization likely possible here because the lists
	# are sorted
	
	reactions = []
	
	# Iterate through all the free domains in reactant1
	for (i, (dom1, strand_num1, dom_num1)) in enumerate(r1_doms):
		# For each, find any domains in reactant2 that could bind
		for (j, (dom2, strand_num2, dom_num2)) in enumerate(r2_doms):
			# If it can pair, this is one possible reaction (this kind of
			# reaction cannot possibly produce a pseudoknotted structure)
			if (dom1.can_pair(dom2)):

				# combine the two complexes into one, but do not perform
				# the association
				reactions.append(combine_complexes_21(
									 reactant1, (strand_num1, dom_num1), 
								     reactant2, (strand_num2, dom_num2)))

									 
	output = []
	for complex, location1, location2 in reactions:

		assert complex.triple(*location1) is not None  
		assert complex.triple(*location2) is not None

		# build "before" and "after" loop structures
		out = find_on_loop(complex, location1, 1, 
			lambda dom1, struct1, loc1, dom2, struct2, loc2: loc1 == location1 and loc2 == location2 )

		[(loc1s, loc2s, before, after)] = out

		# zipper for max-helix semantics
		if UNZIP and not LEGACY_UNZIP:
			(loc1s, loc2s, before, after) = zipper(complex, location1, location2, before.parts, after.parts, 1, 
				lambda dom1, struct1, loc1, dom2, struct2, loc2: struct1 is None and struct2 is None and dom1.can_pair(dom2))

		product = do_bind11(complex, loc1s.locs, loc2s.locs)

		reaction = ReactionPathway('bind21', [reactant1, reactant2], 
								      [product])

		length = len(loc1s)
		reaction._const = bimolecular_binding_rate(length)

		output.append(reaction)
	
	return output
	
	
def find_external_strand_break(complex, location):
	"""
	Takes a complex and a location (strand index, domain index). This location
	refers to a domain in the `complex` which is on an external loop. Returns the 
	index of the last strand _before_ the strand break which would put `location` 
	on an external loop.
	
	Finds the location of a strand break on an external loop.  Used to 
	determine where to split a complex when merging it with another complex at
	`location`.
	"""
	
	strand_num = location[0]
	domain_num = location[1]
	
	stop_flag = False
	search_dom_index = domain_num + 1
	search_strand_index = strand_num
	insertion_strand_index = None
	
	# We will first search to the 'right', looking for the first external
	# strand break
	while insertion_strand_index == None:
		# First check to see if we've run off the end of a strand
		# in which case the external break is between this strand and the
		# next one
		if (search_dom_index == len(complex.strands[search_strand_index])):
			insertion_strand_index = search_strand_index
		else:
			paired_dom = complex.structure[search_strand_index]\
										  [search_dom_index]
			if (paired_dom == None):			
				# If the current domain is unpaired, move to the next domain
				# to the right
				search_dom_index += 1
			
			# elif ((paired_dom[0] > search_strand_index) or \
			# 		 ((paired_dom[0] == search_strand_index) and \
			# 		  (paired_dom[1] > search_dom_index))):
			elif (paired_dom > (search_strand_index,search_dom_index)):
				# If the current domain is paired to something which is to the
				# right of the current search domain, then we jump past that
				# loop
				search_strand_index = paired_dom[0]
				search_dom_index = paired_dom[1] + 1
			else:
				# Otherwise the current domain is paired to something to the
				# left of the current search domain, which means that the
				# external strand break must be to the left
				break
				
	search_dom_index = domain_num - 1
	search_strand_index = strand_num
	
	# We now search to the 'left' if we haven't found the break yet
	while insertion_strand_index == None:
		# First check to see if we've run off the end of a strand, in which
		# case the external break is between the previous strand and this one
		if (search_dom_index == -1):
			insertion_strand_index = search_strand_index - 1
		else:
			paired_dom = complex.structure[search_strand_index] \
										  [search_dom_index]
			if (paired_dom == None):			
				# If the current domain is unpaired, move to the next domain
				# to the left
				search_dom_index -= 1
			# elif ((paired_dom[0] < search_strand_index) or 
			# 		 ((paired_dom[0] == search_strand_index) and
			# 		  (paired_dom[1] < search_dom_index))):
			elif (paired_dom < (search_strand_index,search_dom_index)):
				# If the current domain is paired to something which is to the
				# left of the current search domain, then we jump past that
				# loop
				search_strand_index = paired_dom[0]
				search_dom_index = paired_dom[1] - 1
			else:
				# Otherwise the current domain is paired to something to the
				# left of the current search domain, which means that the
				# external strand break must be to the right; however we
				# checked that, so there must have been an error...
				raise \
				    Exception("Unexpected error in find_external_strand_break")
	
	return insertion_strand_index


def combine_complexes_21(complex1, location1, complex2, location2):
	"""
	Combines two complexes to form one complex, binding the domain in
	complex1 at location1 to the domain in complex2 at location2.
	
	Returns the new complex.
	"""

	# Remember locations on original complexes
	loc1 = location1
	loc2 = location2
	
	# First we need to find the external strand breaks where we will be
	# splitting the complexes
	
	insertion_index_1 = find_external_strand_break(complex1, location1)
	insertion_index_2 = find_external_strand_break(complex2, location2)
		
	# We then find the four parts, 1-4, which we will stick together
	# in order to create our final complex
	
	# If this condition is true, we will need to cut up complex1
	if insertion_index_1 >= 0:
		d1 = complex1.strands[0:(insertion_index_1+1)]
		d4 = complex1.strands[(insertion_index_1+1):]
		s1 = copy.deepcopy(complex1.structure[0:(insertion_index_1+1)])
		s4 = copy.deepcopy(complex1.structure[(insertion_index_1+1):])	
	else:
		# It's not clear to me whether this  condition will ever be executed, 
		# because find_external_strand_break should never return < 0 -CG 5/31
		d1 = complex1.strands[:]
		d4 = []
		s1 = copy.deepcopy(complex1.structure)
		s4 = []
	
	# If this condition is true, we will need to cut up complex2	
	if insertion_index_2 >= 0:
		d2 = complex2.strands[(insertion_index_2+1):]
		d3 = complex2.strands[0:(insertion_index_2+1)]
		s2 = copy.deepcopy(complex2.structure[(insertion_index_2+1):])
		s3 = copy.deepcopy(complex2.structure[0:(insertion_index_2+1)])
	else:
		# Likewise with this one; find_external_strand_break should always
		# produce an index >= 0
		d2 = []
		d3 = complex2.strands[:]
		s2 = []
		s3 = copy.deepcopy(complex2.structure)
	
	# We now need to update the structures based on the way things were shifted
	# We calculate the offsets for strand locations in each chunk
	s1_strand_offset = 0
	s2_strand_offset = len(d1) - len(d3)
	s3_strand_offset = len(d1) + len(d2)
	s4_strand_offset = len(d2) + len(d3)
	
	
	# We then iterate through the structure and update as needed
	if insertion_index_1 >= 0:
		for (strand_index, strand_list) in enumerate(s1):
			for (dom_index, pair) in enumerate(strand_list):
				# If the domain is paired, check if it is paired to something 
				# in s4
				if (pair and (pair[0] > insertion_index_1)):
					new_pair = (pair[0] + s4_strand_offset, pair[1])
					s1[strand_index][dom_index] = new_pair
	
	for (strand_index, strand_list) in enumerate(s2):
		for (dom_index, pair) in enumerate(strand_list):
			# If the domain is paired, check to see if it is paired to something
			# in s2 or s3
			if (pair):
				# Check if paired to s2
				if ((insertion_index_2 >= 0) and (pair[0] > insertion_index_2)):
					new_pair = (pair[0] + s2_strand_offset, pair[1])
					s2[strand_index][dom_index] = new_pair
				# Otherwise paired to s3
				else:
					new_pair = (pair[0] + s3_strand_offset, pair[1])
					s2[strand_index][dom_index] = new_pair
	
		
	for (strand_index, strand_list) in enumerate(s3):
		for (dom_index, pair) in enumerate(strand_list):
			# If the domain is paired, check to see if it is paired to something
			# in s2 or s3
			if (pair):
				# Check if paired to s2
				if ((insertion_index_2 >= 0) and (pair[0] > insertion_index_2)):
					new_pair = (pair[0] + s2_strand_offset, pair[1])
					s3[strand_index][dom_index] = new_pair
				# Otherwise paired to s3
				else:
					new_pair = (pair[0] + s3_strand_offset, pair[1])
					s3[strand_index][dom_index] = new_pair
	
	for (strand_index, strand_list) in enumerate(s4):
		for (dom_index, pair) in enumerate(strand_list):
			# Check if the domain is paired to something in s4
			if (pair and (pair[0] > insertion_index_1)):
				new_pair = (pair[0] + s4_strand_offset, pair[1])
				s4[strand_index][dom_index] = new_pair
	
	new_strands = d1 + d2 + d3 + d4
	new_structure = s1 + s2 + s3 + s4
	
	# Finally, we need to update given reaction indices to new indices
	# and then add the new pairing to the structure
	# if ((insertion_index_1 > 0) and (location1[0] > insertion_index_1)):
	if (location1[0] > insertion_index_1):
		location1 = (location1[0] + s4_strand_offset, location1[1])
	
	
	# NOTE: This case is unexercised by all tests written by Karthik, regardless
	# whether the first clause is included. Since this tests whether the 
	# location2 is within d2 or d3, I think the first clause is unnecessary.

	# if ((insertion_index_2 > 0) and (location2[0] > insertion_index_2)):
	if (location2[0] > insertion_index_2):
		location2 = (location2[0] + s2_strand_offset, location2[1])
	else:
		location2 = (location2[0] + s3_strand_offset, location2[1])
	
	# # do single bind11
	# new_structure[location1[0]][location1[1]] = location2
	# new_structure[location2[0]][location2[1]] = location1
		
	# remember strands participating in binding
	strand1 = complex1.get_strand(loc1[0])
	strand2 = complex2.get_strand(loc2[0])

	# make new complex
	new_complex = Complex(get_auto_name(), new_strands, new_structure)

	# strands may be re-ordered in new complex, so we need to back
	# out where the new strands ended up
	# new_strands = new_complex.strands
	# for index, strand in enumerate(new_strands):
	# 	if strand == strand1: location1 = (index, location1[1])
	# 	if strand == strand2: location2 = (index, location2[1])
	location1 = new_complex.rotate_location(location1)
	location2 = new_complex.rotate_location(location2)

	return new_complex, location1, location2



def do_single_open(reactant, loc):
	new_struct = copy.deepcopy(reactant.structure)
	loc1 = loc
	loc2 = new_struct[loc1[0]][loc1[1]]
	assert new_struct[loc2[0]][loc2[1]] == loc1
	new_struct[loc1[0]][loc1[1]] = None
	new_struct[loc2[0]][loc2[1]] = None
	out = Complex(get_auto_name(), reactant.strands[:], new_struct)
	return out


def open(reactant):
	"""
	Returns a list of reaction product sets that can be produced by the
	'open' reaction, in which a short helix dissociates. Each product
	set are the results of one particular dissociation; each strand in the
	reactant occurs exactly once in one of the complexes in the product set.
	
	A dissociation can happen to any helix under the threshold length		
	"""

	# remember the larger release cutoff; don't enumerate any reactions
	# for helices longer than this
	MAX_RELEASE_CUTOFF = max(RELEASE_CUTOFF_1_1, RELEASE_CUTOFF_1_N)

	reactions = []
	
	structure = reactant.structure
	strands = reactant.strands
	
	# for no-max-helix mode:
	if not UNZIP:
		# We iterate through all the domains
		for (strand_index, strand) in enumerate(reactant.strands):
			for (domain_index, domain) in enumerate(strand.domains):

				# The bound domain must be... bound
				if (structure[strand_index][domain_index] is None):
					continue

				bound_domain = strand.domains[domain_index]
				bound_loc = (strand_index, domain_index)

				release_reactant = do_single_open(reactant, bound_loc)
				product_set = find_releases(release_reactant)
				reactions.append((product_set, len(bound_domain)))

				# # search both directions around the loop for a bound domain that
				# # is complementary (and therefore can be bound to)
				# def criteria(dom1, struct1, loc1, dom2, struct2, loc2):
				# 	return struct1 == loc2 and struct2 == loc1 and dom1.can_pair(dom2)

				# bound_doms = (find_on_loop(reactant, bound_loc, -1, criteria) + 
				# 	find_on_loop(reactant, bound_loc, +1, criteria))	

				# for (displacing, bound, before, after) in bound_doms:
				# 	displacing_loc = list(displacing.locs)[0]
				# 	bound_loc = list(bound.locs)[0]
				# 	release_reactant = do_single_open(reaction, displacing_loc)
				# 	product_set = find_releases(release_reactant)
				# 	reaction = ReactionPathway('open', [reactant], sorted(product_set))

	# for max-helix mode:
	else:
		# We loop through all stands, domains
		for (strand_index, strand) in enumerate(strands):
			for (domain_index, domain) in enumerate(strand.domains):
				# If the domain is unpaired, skip it
				if (structure[strand_index][domain_index] == None):
					continue
				
				# A: Strand/domain position on "top" strand - CG 5/21 
				helix_startA = [strand_index, domain_index]
				
				# B: Strand/domain position on "bottom" strand - CG 5/21 
				helix_startB = list(structure[strand_index][domain_index])
			
				# If the domain is bound to an earlier domain, then we have
				# already considered it, so skip it
				if (helix_startB < helix_startA):
					continue
				
				helix_endA = helix_startA[:]
				helix_endB = helix_startB[:]
							
				helix_length = domain.length
				
				# Now iterate through the whole helix to find the other end
				# of this one
				# (The helix ends at the first strand break from either direction)
				while True:
					
					# Strands run in opposite directions, so A must be incremented 
					# and B decremented in order that both pointers move "right" 
					# along the helix- CG 5/21
					helix_endA[1] += 1
					helix_endB[1] -= 1
					
					# If one of the strands has broken, the helix has ended
					if (helix_endA[1] >= strands[helix_endA[0]].length):
						break
					elif (helix_endB[1] < 0):
						break
					
					
					# If these domains aren't bound to each other, the helix
					# has ended
					if (tuple(helix_endA) != structure[helix_endB[0]][helix_endB[1]]):
						break
					
					# Add the current domain to the current helix
					helix_length += strands[helix_endA[0]].domains[helix_endA[1]]\
										   .length
					
				# We must also iterate in the other direction
				while True:
					
					# Now we move 
					helix_startA[1] -= 1
					helix_startB[1] += 1

					# If one of the strands has broken, the helix has ended
					if (helix_startA[1] < 0):
						break
					elif (helix_startB[1] >= strands[helix_startB[0]].length):
						break
					
					# If these domains aren't bound to each other, the helix
					# has ended
					if (tuple(helix_startA) != structure[helix_startB[0]][helix_startB[1]]):
						break
					
					# Add the current domain to the current helix
					helix_length += strands[helix_startA[0]].domains[helix_startA[1]]\
										   .length
					
				
				# Move start location to the first domain in the helix
				helix_startA[1] += 1
				helix_startB[1] -= 1
				
				# If the helix is short enough, we have a reaction	
				if (helix_length <= MAX_RELEASE_CUTOFF):


					release_reactant = Complex(get_auto_name(), reactant.strands[:], 
						copy.deepcopy(reactant.structure))
					
					# Delete all the pairs in the released helix
					for dom in range(helix_startA[1], helix_endA[1]):
						bound_loc = reactant.structure[helix_startA[0]][dom]
						release_reactant.structure[helix_startA[0]][dom] = None
						release_reactant.structure[bound_loc[0]][bound_loc[1]] = None
			
			
					product_set = find_releases(release_reactant)
					
					
					reactions.append((product_set, helix_length))
		
	output = []
	for product_set,length in reactions:
		reaction = ReactionPathway('open', [reactant], sorted(product_set))
		
		# discard reactions where the release cutoff is greater than the threshold
		if len(reaction.products) == 1 and length > RELEASE_CUTOFF_1_1:
			continue
		elif len(reaction.products) > 1 and length > RELEASE_CUTOFF_1_N:
			continue

		reaction._const = opening_rate(length)
		output.append(reaction)
	
	output = sorted(list(set(output)))

	return output
	
	
def find_releases(reactant):
	"""
	Determines if reactant actually represents multiple unattached complexes.
	If so, returns a list of these complexes. Otherwise, returns [reactant].
	"""
	
	structure = reactant.structure
	strands = reactant.strands
	
	output_list = []
	
	# Iterate through the strands and determine if a split is necessary
	# Note that we don't iterate to the last strand because if the last
	# strand were free from the rest, we would catch that as the rest
	# of the strands being free from the last strand. Hence, we also
	# terminate loops if we ever get to the last strand
	for (strand_index, strand) in enumerate(strands[:-1]):
		
		domain_index = len(strand.domains)
		
		inner_index = (strand_index, domain_index - 1)
		
		# We now iterate through lower domains and see if we can find
		# a release point (Iterate while:
		#	0 > inner_index strand > index of last strand, and
		# 		inner_index strand < strand_index, or
		#		inner_index strand = strand_index and inner_index domain < domain_index)
		# while (inner_index[0] >= 0) and (inner_index[0] < (len(strands) - 1)) and \
		# 	  ( inner_index < (strand_index, domain_index)):
		while ( 0 <= inner_index[0] < (len(strands) - 1) ) and \
			  ( inner_index < (strand_index, domain_index) ):


			# If we have run off of the end of a strand,
			# then we have found a release point  
			if (inner_index[1] == -1):
				split_start = inner_index
				split_end = (strand_index, domain_index)
				split_list = split_complex(reactant, split_start, split_end)
				# We check the two resulting complexes to see if they can
				# be split further
				for complex in split_list:
					output_list += (find_releases(complex))
				return output_list					
			
			# Otherwise decide where to go next
			curr_structure = structure[inner_index[0]][inner_index[1]]
						
			# If this domain is unpaired, move to the next
			if (curr_structure == None):
				inner_index = (inner_index[0], inner_index[1] - 1)
			
			# Check if the structure points to a higher domain 
			# elif ((curr_structure[0] > inner_index[0]) or \
			#      ((curr_structure[0] == inner_index[0]) and \
			# 	  (curr_structure[1] > inner_index[1]))):
			elif (curr_structure > inner_index):

				# If the structure points to a domain above the start, this section
				# is connected to something higher, so abort this loop
				#if (curr_structure[0] > strand_index) and (curr_structure[1] > 0):
				if (curr_structure[0] > strand_index) and (curr_structure[1] >= 0):
					break
				  				  
				# Otherwise it points to a domain between this one and the start
				# -- we've already been there so move to the next
				inner_index = (inner_index[0], inner_index[1] - 1)			
				
			# Otherwise, follow the structure
			else:
				inner_index = curr_structure


		# Alright; the following statement precludes the while loop below from 
		# ever running, since the second clause will always fail. 
		# Then again, without ever running this while loop, all other tests pass.
		# OK, I'm fairly convinced this isn't needed, since the above loop 
		# effectively works by iterating (forwards) across each strand, then looking 
		# backwards -CG 5/29
		
		
		##### TODO: FIND OUT IF THIS IS NEEDED...
#		inner_index = (strand_index, 0)
#				
#		# If we didn't find a release point in the lower domains,
#		# we now try iterating through the higher domains
#		while (inner_index[0] < len(strands) - 1) and \
#			  ((inner_index[0] > strand_index) or \
#			  ((inner_index[0] == strand_index) and \
#			  (inner_index[1] > domain_index))):
#			# If we have run off of the end of a strand,
#			# then we have found a release point  
#			if (inner_index[1] == len(strands[inner_index[0]].domains)):
#				split_start = (strand_index, domain_index)
#				split_end = inner_index
#				split_list = split_complex(reactant, split_start, split_end)
#				# We check the two resulting complexes to see if they can
#				# be split further
#				for complex in split_list:
#					output_list.extend(find_releases(complex))
#				return output_list					
#							
#			# Otherwise decide where to go next
#			curr_structure = structure[inner_index[0]][inner_index[1]]	
#			
#			# If this domain is unpaired, move to the next
#			if (curr_structure == None):
#				inner_index = (inner_index[0], inner_index[1] + 1)
#
#			# Check if the structure points to a lower domain			
#			elif (curr_structure[0] < inner_index[0]) or \
#				 ((curr_structure[0] == inner_index[0]) and \
#				  (curr_structure[1] < inner_index[1])):
#				  
#				# If the structure points to a domain before the start,
#				# this section is connected to something lower, so abort this
#				# loop
#				if (curr_structure[0] < strand_index):
#					break
#				
#				# Otherwise it points to a domain between this one and the start
#				# -- we've already been there to move to the next
#				inner_index = (inner_index[0], inner_index[1] + 1)			
#			
#			
#			# Otherwise, follow the structure
#			else:
#				inner_index = curr_structure			
				
	# If we still haven't found any splits, then this complex cannot be split
	
	return [reactant]
	
def split_complex(reactant, split_start, split_end):
	"""
	Splits a disconnected complex between split_start and split_end (which
	should be at ends of strands), and returns the two resulting complexes.
	"""
	
	strands = reactant.strands
	structure = reactant.structure
	
	# We first take the parts apart
	strands1 = strands[0:split_start[0]]
	structure1 = structure[0:split_start[0]]
	
	strands2 = strands[split_start[0]:split_end[0]+1]
	structure2 = structure[split_start[0]:split_end[0]+1]
	
	strands3 = strands[split_end[0]+1:]
	structure3 = structure[split_end[0]+1:]
	
	out1_strands = strands1 + strands3
	out2_strands = strands2
		
		
	# These offsets are the changes to the strand numbers needed
	structure2_offset = len(strands1)
	structure3_offset = len(strands2)

	# Apply offset to structure2
	new_structure1 = []
			
	for strand in structure1:
		strand_out = []
		for struct_element in strand:
			if struct_element == None:
				strand_out.append(None)
			elif (struct_element[0] > split_end[0]):
				strand_out.append((struct_element[0] - structure3_offset, struct_element[1]))
			else:
				strand_out.append(struct_element)
		new_structure1.append(strand_out)
	
			
	# Apply offset to structure2
	new_structure2 = []
		
	for strand in structure2:
		strand_out = []
		for struct_element in strand:
			if struct_element == None:
				strand_out.append(None)
			else:
				strand_out.append((struct_element[0] - structure2_offset, struct_element[1]))
		new_structure2.append(strand_out)
		
	# Apply offset to structure3
	new_structure3 = []
	
	
	for strand in structure3:
		strand_out = []
		for struct_element in strand:
			if struct_element == None:
				strand_out.append(None)
			elif (struct_element[0] < split_start[0]):
				strand_out.append(struct_element)
			else:
				strand_out.append((struct_element[0] - structure3_offset, struct_element[1]))
		new_structure3.append(strand_out)
		
	out1_structure = new_structure1 + new_structure3
	out2_structure = new_structure2
	
	
	out1 = Complex(get_auto_name(), out1_strands, out1_structure)
	out2 = Complex(get_auto_name(), out2_strands, out2_structure)
	
	return [out1, out2]


def domains_adjacent(loc1, loc2):
	"""
	Returns true if the two locations represent adjacent domains;
	each location should be a (strand_index, domain_index) pair
	"""
	return (loc1[0] == loc2[0]) and (abs(loc1[0] - loc2[0]) == 1)

def branch_3way(reactant):
	"""
	Returns a list of reaction pathways that can be created through one 
	iteration of a 3 way branch migration reaction (more than one molecule may 
	be produced by a reaction because branch migration can liberate strands and 
	complexes).
	"""
	
	reactions = []
	structure = reactant.structure
	
	# We iterate through all the domains
	for (strand_index, strand) in enumerate(reactant.strands):
		for (domain_index, domain) in enumerate(strand.domains):

			# The displacing domain must be free			
			if (structure[strand_index][domain_index] != None):
				continue

			displacing_domain = strand.domains[domain_index]
			displacing_loc = (strand_index, domain_index)

			# search both directions around the loop for a bound domain that
			# is complementary (and therefore can be bound to)
			def criteria(dom1, struct1, loc1, dom2, struct2, loc2):
				return struct1 is None and struct2 is not None and dom1.can_pair(dom2)

			bound_doms = (find_on_loop(reactant, displacing_loc, -1, criteria) + 
				find_on_loop(reactant, displacing_loc, +1, criteria))


			# # search both directions around the loop for a bound domain that
			# # has the same sequence (and therefore can be displaced)
			# bound_doms = find_on_loop(reactant, displacing_loc, -1, \
			# 	lambda dom1, struct1, loc1, dom2, struct2, loc2: struct2 is not None and dom1 == dom2) + \
			# find_on_loop(reactant, displacing_loc, +1, \
			# 	lambda dom1, struct1, loc1, dom2, struct2, loc2: struct2 is not None and dom1 == dom2)

			# bound_doms = find_on_loop(reactant, displacing_loc, -1, \
			# 	lambda dom, struct, loc: struct is not None and dom == displacing_domain) + \
			# find_on_loop(reactant, displacing_loc, +1, \
			# 	lambda dom, struct, loc: struct is not None and dom == displacing_domain)


			# build products
			# [ (Loop([triple(start_loc)]), Loop([triple(bound_loc]), Loop(loop[:i]), Loop(loop[i+1:])) ]
			for (displacing, bound, before, after) in bound_doms:

				if UNZIP and LEGACY_UNZIP:
					displacing_loc = list(displacing.locs)[0]
					bound_loc = list(bound.locs)[0]
					reaction = ReactionPathway('branch_3way', [reactant], do_3way_migration_legacy(
						reactant, 
						displacing_loc, 
						bound_loc)
					)
				else:
					reaction = ReactionPathway('branch_3way', [reactant], do_3way_migration(
						reactant, displacing.locs, 
						bound.locs)
					)

				# length of invading domain
				length = len(displacing)

				# calculate reaction constant
				(after, before) = (before, after)
				reaction._const = branch_3way_remote_rate(length, before, after)
				# reaction._const = branch_3way_remote_rate(length, after, before)

				# skip remote toehold reactions if directed
				if REJECT_REMOTE:
					if not (not after.is_open and after.stems==1 and after.bases==0):
						# print "Rejecting... " + reaction.kernel_string()
						# import pdb; pdb.set_trace()
						continue

				reactions.append(reaction)

			# for (bound_loc, before, after) in bound_doms:
			# 	reaction = ReactionPathway('branch_3way', [reactant], do_3way_migration(\
			# 		reactant, displacing_loc, structure[bound_loc[0]][bound_loc[1]])
			# 	)

			# 	# length of invading domain
			# 	length = len(displacing_domain)

			# 	# calculate reaction constant
			# 	reaction._const = branch_3way_remote_rate(length, before, after)
			# 	# reaction._const = branch_3way_rate(length)

			# 	reactions.append(reaction)
				
	output = reactions

	# Remove any duplicate reactions
	output = sorted(list(set(output)))

	return output

def do_single_3way_migration(reactant, displacing_loc, new_bound_loc):
	"""
	displacing_loc will be bound to new_bound_loc; whatever new_bound_loc
	was bound to will be unbound.
	"""
	struct = copy.deepcopy(reactant.structure)
	displaced_loc = struct[new_bound_loc[0]][new_bound_loc[1]]

	assert struct[displacing_loc[0]][displacing_loc[1]] is None
	assert struct[new_bound_loc[0]][new_bound_loc[1]] is not None
	assert struct[displaced_loc[0]][displaced_loc[1]] is not None
	assert struct[displaced_loc[0]][displaced_loc[1]] == new_bound_loc

	struct[displacing_loc[0]][displacing_loc[1]] = new_bound_loc
	struct[new_bound_loc[0]][new_bound_loc[1]] = displacing_loc
	struct[displaced_loc[0]][displaced_loc[1]] = None
	
	out_reactant = Complex(get_auto_name(),reactant.strands[:], struct)

	return out_reactant


def do_3way_migration(reactant, displacing_locs, bound_locs):
	"""
	Each location in displacing_locs will end up bound to the corresponding
	location in bound_locs. The stuff bound to bound_locs will end up un-bound
	"""
	if isinstance(displacing_locs, tuple): displacing_locs = [displacing_locs]
	if isinstance(bound_locs, tuple):	bound_locs = [bound_locs]

	product = reactant
	for displacing_loc, bound_loc in zip(displacing_locs, bound_locs):
		product = do_single_3way_migration(product, displacing_loc, bound_loc)

	return find_releases(product)

def do_3way_migration_legacy(reactant, displacing_loc, new_bound_loc):
	"""
	Returns the product set which is the result of a 3-way branch migration
	reaction where the domain at displacing_loc displaces the domain bound to
	the domain at new_bound_loc.
	"""

	# out_reactant = copy.deepcopy(reactant)

	out_reactant_structure = copy.deepcopy(reactant.structure)

	out_reactant_structure[displacing_loc[0]][displacing_loc[1]] = new_bound_loc
	displaced_loc = out_reactant_structure[new_bound_loc[0]][new_bound_loc[1]]
	out_reactant_structure[new_bound_loc[0]][new_bound_loc[1]] = displacing_loc
	out_reactant_structure[displaced_loc[0]][displaced_loc[1]] = None
	
	out_reactant = Complex(get_auto_name(),reactant.strands[:], out_reactant_structure)
	
	global UNZIP
	
	# Check to see if an adjacent displacement is possible
	if (UNZIP):
		dstrand = displacing_loc[0]
		ddomain = displacing_loc[1]
		bstrand = new_bound_loc[0]
		bdomain = new_bound_loc[1]
		if (ddomain+1 < len(out_reactant.strands[dstrand].domains)) and \
			(out_reactant.structure[dstrand][ddomain+1] == None) and \
			(bdomain-1 >= 0) and \
			(out_reactant.structure[bstrand][bdomain-1] != None) and \
			(out_reactant.strands[bstrand].domains[bdomain-1].can_pair(out_reactant.strands[dstrand].domains[ddomain+1])):
			
			return do_3way_migration_legacy(out_reactant, (dstrand, ddomain+1), (bstrand, bdomain-1))

		elif (ddomain-1 >= 0) and \
			(out_reactant.structure[dstrand][ddomain-1] == None) and \
			(bdomain+1 < len(out_reactant.strands[bstrand].domains)) and \
			(out_reactant.structure[bstrand][bdomain+1] != None) and \
			(out_reactant.strands[bstrand].domains[bdomain+1].can_pair(out_reactant.strands[dstrand].domains[ddomain-1])):
			
			return do_3way_migration_legacy(out_reactant, (dstrand, ddomain-1), (bstrand, bdomain+1))
		else:
			return find_releases(out_reactant)
	else:
		return find_releases(out_reactant)


# def find_on_loop(reactant, start_loc, direction, filter):
# 	"""
# 	Finds the next domain within `reactant` that's on the same inner loop as 
# 	`start_loc` and matches the passed `filter` function. Looks in either the 
# 	5'->3' (+1) or 3'->5' (-1) `direction`.  

# 	Filter should accept the following arguments and return True or False:
# 		-	dom (utils.Domain) : the domain at `loc`
# 		-	struct (tuple or None): a (strand index, domain index) pair 
# 			indicating what `dom` is bound to, or None if `dom` is unpaired.
# 		-	loc (tuple) : a (strand index, domain index) pair
# 		-       Note that while every single-stranded domain is tested,
# 		        only the "first" domain of a stem helix (in the direction of 
# 			search) will be passed to the filter.


# 	Returns an array of tuples: `(loc, before, after)`, where:
# 		-	`loc` is a (strand index, domain index) pair indicating the 
# 			position of the matched domain
# 		-	`before` is a list of (domain, struct, loc) triples giving the 
# 			domains after `start_loc` but before the matched domain on the loop
# 			(or None instead of triple where there is a break in the loop)
# 		-	`after` is a list of (domain, struct, loc) triples giving the 
# 			domains after the matched domain but before `start_loc` on the loop
# 			(or None instead of triple where there is a break in the loop)

# 	Where a loop involves stems, only one of the complementary domains will be 
# 	listed in the array of tuples, specifically, the "first" one in the search 
# 	direction. Thus, a multiloop with n unpaired domains and m stems will 
# 	result, for closed loops, in `len(before+after) == n+m-2`, as the match 
# 	location and `start_loc` are omitted.

# 	`before` and `after` are converted to Loop objects (see utils.py) prior 
# 	to being returned, so that the number of bases and number of stems and 
# 	open/closed status is readily accessible.

# 	Note 1: `before` and `after` refer to the partial loops between `start_loc` 
# 	and each of the results, _in the `direction`_ of the search. For example:
  
#            A
# 	      ____
# 	     /    \ 
# 	 x  |     |  x*
# 	    |
# 	     \____> 3'

# 	        B

# 	If `start_loc` pointed to `x` and `direction` is +1, then `before` would 
# 	be `A` and `after` would be `B`. If instead `direction` is -1, then 
# 	`before` is `B` and `after` is `A`. 

# 	Note 2: If the domain passed to `start_loc` is a duplex, the results may 
# 	be unexpected:

# 	       ___  x  ___
# 	5' ___/   \___/   \___  
# 	3' ___  A  ___  B  ___)
# 	      \___/   \___/
# 	            x*

# 	Notice that the duplex x() participates in two internal loops (A and B). 
# 	By convention, the internal loop considered is the _internal loop which
# 	encloses this domain_. That means if you pass domain x and +1, you'll get
# 	loop A, whereas if you pass x and -1, you'll get loop B. This is in an
# 	attempt to be consistent with the case where you pass an unpaired domain
# 	(and therefore the internal loop searched is the one which encloses the 
# 	unpaired domain).
# 	"""	
# 	results = []
# 	loop = []

# 	def triple(loc):
# 		return (reactant.get_domain(loc),reactant.get_structure(loc),loc)

# 	# We now follow the external loop from the starting pair
# 	# searching for a bound domain to displace
# 	bound_loc = start_loc

# 	# Avoid getting stuck inside an internal loop enclosed by this domain,
# 	# if the starting domain is a duplex.
# 	# 
# 	#   1      2
# 	#  ___________
# 	#  ____  _____
# 	#   1*  /  2*
# 	#       
# 	#  If we start at domain 1, going in the - direction, then 
# 	#  immediately continue to the next domain, we'll go to 1* 
# 	if reactant.structure[bound_loc[0]][bound_loc[1]] is not None:
# 		bound_loc = reactant.structure[bound_loc[0]][bound_loc[1]]

# 	# Follow the external loop to the end
# 	while True:
# 		# move to the next domain in the indicated direction
# 		# (+1 = 5' -> 3', -1 = 3' -> 5')
# 		bound_loc = (bound_loc[0], bound_loc[1] + direction)

# 		# if we've reached the end of the strand (5')
# 		if (bound_loc[1] == -1):

# 			# Continue to next strand
# 			bound_loc = (wrap(bound_loc[0]-1,len(reactant.strands)),)
# 			bound_loc = (bound_loc[0], len(reactant.strands[bound_loc[0]]))
# 			loop.append( None )  #EW
# 			continue

# 		# if we've reached the end of the strand (3')
# 		elif (bound_loc[1] == len(reactant.strands[bound_loc[0]])):

# 			# Continue to next strand
# 			bound_loc = (wrap(bound_loc[0]+1, len(reactant.strands)), -1)
# 			loop.append( None ) #EW
# 			continue

# 		if bound_loc == start_loc:
# 			# We've returned to the original location of the 
# 			# displacing domain
# 			break

		
# 		# try to match the filter function
# 		elif (filter(reactant.get_domain(bound_loc), \
# 			reactant.get_structure(bound_loc), bound_loc)):

# 			# append the location
# 			results.append( (bound_loc, len(loop)) )

# 		# store unpaired domains and "first" domain of each stem
# 		loop.append(triple(bound_loc)) #EW 

# 		# if the domain at bound_loc is unbound
# 		if (reactant.structure[bound_loc[0]][bound_loc[1]] is None):
# 			# look to the next domain
# 			continue

# 		# so it's bound to something: follow the structure to stay on the same loop
# 		bound_loc = reactant.structure[bound_loc[0]][bound_loc[1]]

# 	return list( (bound_loc, Loop(loop[:i]), Loop(loop[i+1:]) ) for (bound_loc, i) in results )  #EW

def find_on_loop(reactant, start_loc, direction, filter):
	"""
	Finds the next domain within `reactant` that's on the same inner loop as 
	`start_loc` and matches the passed `filter` function. Looks in either the 
	5'->3' (+1) or 3'->5' (-1) `direction`.  

	Filter should accept the following arguments and return True or False:
		-	dom (utils.Domain) : the domain at `loc`
		-	struct (tuple or None): a (strand index, domain index) pair 
			indicating what `dom` is bound to, or None if `dom` is unpaired.
		-	loc (tuple) : a (strand index, domain index) pair
		-       Note that while every single-stranded domain is tested,
		        only the "first" domain of a stem helix (in the direction of 
			search) will be passed to the filter.


	Returns an array of tuples: `(loc, before, after)`, where:
		-	`loc` is a (strand index, domain index) pair indicating the 
			position of the matched domain
		-	`before` is a list of (domain, struct, loc) triples giving the 
			domains after `start_loc` but before the matched domain on the loop
			(or None instead of triple where there is a break in the loop)
		-	`after` is a list of (domain, struct, loc) triples giving the 
			domains after the matched domain but before `start_loc` on the loop
			(or None instead of triple where there is a break in the loop)

	Where a loop involves stems, only one of the complementary domains will be 
	listed in the array of tuples, specifically, the "first" one in the search 
	direction. Thus, a multiloop with n unpaired domains and m stems will 
	result, for closed loops, in `len(before+after) == n+m-2`, as the match 
	location and `start_loc` are omitted.

	`before` and `after` are converted to Loop objects (see utils.py) prior 
	to being returned, so that the number of bases and number of stems and 
	open/closed status is readily accessible.

	Note 1: `before` and `after` refer to the partial loops between `start_loc` 
	and each of the results, _in the `direction`_ of the search. For example:
  
           A
	      ____
	     /    \ 
	 x  |     |  x*
	    |
	     \____> 3'

	        B

	If `start_loc` pointed to `x` and `direction` is +1, then `before` would 
	be `A` and `after` would be `B`. If instead `direction` is -1, then 
	`before` is `B` and `after` is `A`. 

	Note 2: If the domain passed to `start_loc` is a duplex, the results may 
	be unexpected:

	       ___  x  ___
	5' ___/   \___/   \___  
	3' ___  A  ___  B  ___)
	      \___/   \___/
	            x*

	Notice that the duplex x() participates in two internal loops (A and B). 
	By convention, the internal loop considered is the _internal loop which
	encloses this domain_. That means if you pass domain x and +1, you'll get
	loop A, whereas if you pass x and -1, you'll get loop B. This is in an
	attempt to be consistent with the case where you pass an unpaired domain
	(and therefore the internal loop searched is the one which encloses the 
	unpaired domain).
	"""	
	results = []
	loop = []

	def triple(loc):
		return (reactant.get_domain(loc),reactant.get_structure(loc),loc)

	# We now follow the external loop from the starting pair
	# searching for a bound domain to displace
	bound_loc = start_loc

	# Avoid getting stuck inside an internal loop enclosed by this domain,
	# if the starting domain is a duplex.
	# 
	#   1      2
	#  ___________
	#  ____  _____
	#   1*  /  2*
	#       
	#  If we start at domain 1, going in the - direction, then 
	#  immediately continue to the next domain, we'll go to 1* 
	if reactant.structure[bound_loc[0]][bound_loc[1]] is not None:
		bound_loc = reactant.structure[bound_loc[0]][bound_loc[1]]

	# Follow the external loop to the end
	while True:
		# move to the next domain in the indicated direction
		# (+1 = 5' -> 3', -1 = 3' -> 5')
		bound_loc = (bound_loc[0], bound_loc[1] + direction)

		# if we've reached the end of the strand (5')
		if (bound_loc[1] == -1):

			# Continue to next strand
			bound_loc = (wrap(bound_loc[0]-1,len(reactant.strands)),)
			bound_loc = (bound_loc[0], len(reactant.strands[bound_loc[0]]))
			loop.append( None )  #EW
			continue

		# if we've reached the end of the strand (3')
		elif (bound_loc[1] == len(reactant.strands[bound_loc[0]])):

			# Continue to next strand
			bound_loc = (wrap(bound_loc[0]+1, len(reactant.strands)), -1)
			loop.append( None ) #EW
			continue

		if bound_loc == start_loc:
			# We've returned to the original location of the 
			# displacing domain
			break

		
		# try to match the filter function
		elif (filter(
			reactant.get_domain(start_loc),
			reactant.get_structure(start_loc),
			start_loc,
			reactant.get_domain(bound_loc), 
			reactant.get_structure(bound_loc), 
			bound_loc)):

			# append the location
			results.append( (bound_loc, len(loop)) )

		# store unpaired domains and "first" domain of each stem
		loop.append(triple(bound_loc)) #EW 

		# if the domain at bound_loc is unbound
		if (reactant.structure[bound_loc[0]][bound_loc[1]] is None):
			# look to the next domain
			continue

		# so it's bound to something: follow the structure to stay on the same loop
		bound_loc = reactant.structure[bound_loc[0]][bound_loc[1]]


	if UNZIP and not LEGACY_UNZIP:
		zipped_results = []
		for (bound_loc, i) in results:
			zipped_results.append( zipper(reactant, start_loc, bound_loc, loop[:i], loop[i+1:], direction, filter) )
		return zipped_results
	else:
		return [(Loop([triple(start_loc)]), Loop([triple(bound_loc)]), Loop(loop[:i]), Loop(loop[i+1:])) for (bound_loc, i) in results]
		# return list( (bound_loc, Loop(loop[:i]), Loop(loop[i+1:]) ) for (bound_loc, i) in results )  #EW
			

def zipper(reactant, start_loc, bound_loc, before, after, direction, filter):
	"""
	Takes a result from `find_on_loop` and "zips" it inwards (in the given 
	`direction`); that is, given some `start_loc` and some `bound_loc`, tries 
	to find as many adjacent domains as possible such that the `filter` 
	function still returns True.


	For example, if `start_loc` was b1 and `bound_loc` was b1*, and the filter
	function specified that the domain at `start_loc` must be complementary to 
	the domain at `bound_loc`, then the function would return [b1,b2] as 
	start_locs and [b1*, b2*] as bound_locs


	        b1* b2*
	        ______
	     __/      \__>
	    <__        __
	       \______/
	        b1  b2


	return start_locs, bound_locs, before, after

	"""
	def triple(loc):
		return (reactant.get_domain(loc),reactant.get_structure(loc),loc)

	def move_towards_middle(middle, direction):

		dstrand, ddomain = start_loc
		bstrand, bdomain = bound_loc

		while True:
			# move domain pointers "inwards" towards each other
			ddomain += direction
			bdomain -= direction

			# if ddomain is still on dstrand 
			if ((ddomain < len(reactant.strands[dstrand].domains) and ddomain >= 0) and 

				# and bdomain is still on bstrand 
				(bdomain < len(reactant.strands[bstrand].domains) and bdomain >= 0) and 

				# and ddomain hasn't passed bound_loc 
				# (cmp((dstrand, ddomain), bound_loc) == -direction ) and 
				(cmp((dstrand, ddomain), bound_loc) == cmp(start_loc, bound_loc) ) and 

				# and bdomain hasn't passed start_loc
				# (cmp((bstrand, bdomain), start_loc) == +direction ) and 
				(cmp((bstrand, bdomain), start_loc) == cmp(bound_loc, start_loc) ) and 

				# and filter condition still applies
				filter( 
					reactant.get_domain((dstrand, ddomain)),
					reactant.get_structure((dstrand, ddomain)),
					(dstrand, ddomain),
					reactant.get_domain((bstrand, bdomain)), 
					reactant.get_structure((bstrand, bdomain)), 
					(bstrand, bdomain))):

				# add new positions to list
				if direction == 1:
					start_locs.append( triple((dstrand, ddomain)) )
					bound_locs.append( triple((bstrand, bdomain)) )
				elif direction == -1:
					start_locs[:0] = [ triple((dstrand, ddomain)) ]
					bound_locs[:0] = [ triple((bstrand, bdomain)) ]

				# remove zipped positions from `middle` loop
				displacing_index = (direction-1)/2
				bound_index =   (-direction-1)/2

				try:
					if middle[displacing_index] is not None and middle[displacing_index][2] == (dstrand, ddomain):
						del middle[displacing_index]
				except IndexError: pass

				try:
					if middle[bound_index] is not None and middle[bound_index][2] == (bstrand, bdomain):
						del middle[bound_index]
				except IndexError: pass

			else: break

	start_locs = [triple(start_loc)]
	bound_locs = [triple(bound_loc)]

	for (d, middle) in [ (direction, before), (-direction, after) ]:
		move_towards_middle(middle, d)

	start_locs = Loop(start_locs)
	bound_locs = Loop(bound_locs)
	before = Loop(before)
	after = Loop(after)
	
	return start_locs, bound_locs, before, after		

# def zip(reactant, loc1, before, loc2, after, filter):
# 	dstrand = displacing_loc[0]
# 	ddomain = displacing_loc[1]
# 	bstrand = new_bound_loc[0]
# 	bdomain = new_bound_loc[1]
# 	if (ddomain+1 < len(reactant.strands[dstrand].domains)) and \
# 		(reactant.structure[dstrand][ddomain+1] == None) and \
# 		(bdomain-1 >= 0) and \
# 		(reactant.structure[bstrand][bdomain-1] != None) and \
# 		(reactant.strands[bstrand].domains[bdomain-1].can_pair(reactant.strands[dstrand].domains[ddomain+1])):
		
# 		return do_3way_migration(reactant, (dstrand, ddomain+1), (bstrand, bdomain-1))

# 	elif (ddomain-1 >= 0) and \
# 		(reactant.structure[dstrand][ddomain-1] == None) and \
# 		(bdomain+1 < len(reactant.strands[bstrand].domains)) and \
# 		(reactant.structure[bstrand][bdomain+1] != None) and \
# 		(reactant.strands[bstrand].domains[bdomain+1].can_pair(reactant.strands[dstrand].domains[ddomain-1])):
		
# 		return do_3way_migration(reactant, (dstrand, ddomain-1), (bstrand, bdomain+1))
# 	else:
# 		return


def branch_4way(reactant):
	"""
	Returns a list of complex sets that can be created through one iteration of
	a 4 way branch migration reaction (each set consists of the molecules that
	result from the iteration; more than one molecule may result because branch
	migration can liberate strands and complexes).	
	"""
	
	structure = reactant.structure
	reactions = []
	
	# We loop through all domains
	for (strand_index, strand) in enumerate(reactant.strands):
		for (domain_index, domain) in enumerate(strand.domains):
			
			# Unbound domains can't participate in branch migration
			if (structure[strand_index][domain_index] == None):
				continue

			displacing_domain = strand.domains[domain_index]
			displacing_loc = (strand_index, domain_index)

			# search both directions around the loop for a bound domain that
			# has the same sequence (and therefore can be displaced)
			# 
			#   z  ___  z* (displacing)
			#  ___/   \___>
			# <___     ___
			#     \___/
			#   z*      z
			# 
			bound_doms = find_on_loop(reactant, displacing_loc, +1, \
				# lambda dom, struct, loc: struct is not None and dom == displacing_domain)
				lambda dom1, struct1, loc1, dom2, struct2, loc2: struct2 is not None and dom1 == dom2)

			# bound_doms = find_on_loop(reactant, structure[strand_index][domain_index], -1, \
			# 	lambda dom, struct, loc: struct is not None and dom.can_pair(displacing_domain)) + \
			# find_on_loop(reactant, displacing_loc, +1, \
			# 	lambda dom, struct, loc: struct is not None and dom == displacing_domain))

			# build products
			for (displacing, displaced, before, after) in bound_doms:
				reaction = ReactionPathway('branch_4way', [reactant], do_4way_migration(
					reactant, 
					displacing.locs, 
					(structure[displacing_loc[0]][displacing_loc[1]] for displacing_loc in displacing.locs),
					(structure[bound_loc[0]][bound_loc[1]] for bound_loc in displaced.locs),
					displaced.locs)
				)

				# length of invading domain
				length = len(displacing)

				# calculate reaction constant
				reaction._const = branch_4way_remote_rate(length, before, after)

				# skip remote toehold reactions if directed
				if REJECT_REMOTE:
					if not (not after.is_open and after.stems==1 and after.bases==0 and 
						not before.is_open and before.stems==1 and before.bases==0):
						continue


				reactions.append(reaction)


			# for (bound_loc, before, after) in bound_doms:
			# 	reaction = ReactionPathway('branch_4way', [reactant], do_4way_migration(\
			# 		reactant, 
			# 		displacing_loc, structure[displacing_loc[0]][displacing_loc[1]],
			# 		structure[bound_loc[0]][bound_loc[1]],bound_loc)
			# 	)

			# 	# length of invading domain
			# 	length = len(displacing_domain)

			# 	# calculate reaction constant
			# 	reaction._const = branch_4way_remote_rate(length, before, after)

			# 	reactions.append(reaction)


	output = reactions
	
	# remove any duplicate reactions	
	output = sorted(list(set(output)))

	return output

def do_single_4way_migration(reactant, loc1, loc2, loc3, loc4):
	"""
	Performs a 4 way branch migration on a copy of reactant, with loc1 as the
	displacing domain, loc2 as the domain displaced from loc1, loc3 as the
	template domain, and loc4 as the domain displaced from loc3. Returns the
	set of complexes produced by this reaction (may be one or more complexes).

	loc1:loc2, loc3:loc4 -> loc1:loc3, loc2:loc4
	"""
	
	new_struct = copy.deepcopy(reactant.structure)
	new_struct[loc1[0]][loc1[1]] = loc3
	new_struct[loc3[0]][loc3[1]] = loc1
	new_struct[loc2[0]][loc2[1]] = loc4
	new_struct[loc4[0]][loc4[1]] = loc2
	
	out = Complex(get_auto_name(), reactant.strands[:], new_struct)
	
	return out

def do_4way_migration(reactant, loc1s, loc2s, loc3s, loc4s):
	if isinstance(loc1s, tuple): loc1s = [loc1s]
	if isinstance(loc2s, tuple): loc2s = [loc2s]
	if isinstance(loc3s, tuple): loc3s = [loc3s]
	if isinstance(loc4s, tuple): loc4s = [loc4s]

	product = reactant
	# for i in xrange(len(loc1s)):
	# 	product = do_single_4way_migration(product, loc1s[i], loc2s[i], loc3s[i], loc4s[i])
	for loc1, loc2, loc3, loc4 in zip(loc1s, loc2s, loc3s, loc4s):
		product = do_single_4way_migration(product, loc1, loc2, loc3, loc4)
	return find_releases(product)

# def do_4way_migration(reactant, loc1, loc2, loc3, loc4):
# 	"""
# 	Performs a 4 way branch migration on a copy of reactant, with loc1 as the
# 	displacing domain, loc2 as the domain displaced from loc1, loc3 as the
# 	template domain, and loc4 as the domain displaced from loc3. Returns the
# 	set of complexes produced by this reaction (may be one or more complexes).

# 	loc1:loc2, loc3:loc4 -> loc1:loc3, loc2:loc4


# 	"""
	
# 	new_struct = copy.deepcopy(reactant.structure)
# 	new_struct[loc1[0]][loc1[1]] = loc3
# 	new_struct[loc3[0]][loc3[1]] = loc1
# 	new_struct[loc2[0]][loc2[1]] = loc4
# 	new_struct[loc4[0]][loc4[1]] = loc2
	
# 	out = Complex(get_auto_name(), reactant.strands[:], new_struct)
	
# 	return find_releases(out)

