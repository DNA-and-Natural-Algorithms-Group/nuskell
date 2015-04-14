
import copy
import re

##################################
# Next stuff to do
# Figure out how properties of the domains will be stored
# Probably as a dictionary of objects that hold all the properties
# or something like that.
#
# Implement handling of lengths of domains. How will this be specified?
# Look at DNA compiler pages on wiki....use the standard!!!
# 
# Need to implement the maximal displacement and maximal binding reactions
# Without these, there will continue to be a state-space explosion and
# dividing the domains will continue to yield different results when it 
# crosses the release boundary
#
# Need to have a way to specify the depth/release length and all that 
# sort of thing in the file. At least make it an instance variable of 
# the Enumerator object. Ugh.
#
# This is getting ugly....I really should just rewrite the whole thing.
# Notes of stuff I hate about how I did this
# 1) Should have kept track of strands in the complexes. It would make
#     all the comparisons a lot easier. Also would have made the naming
#     of complexes more clear. 
# 2) Should have made the complexes start or end with a strand break. Would
#     Have made everything in searching for reactions a lot easier and more
#     intuitive.
# 3) Alternatively, should have thought up a better way to store the structure
#     No ideas on this one, storing them as indexing variables seems like the
#     ideal way to go to me, but the idea is there.
# 4) Should have made the definition of Complex and Enumerator more distinct.
#     The functionality of each is overlapping, could be segregated a bit
#     better. It is a little better now, but still could be cleaner.
# 5) 

class Complex:
  def __init__(self, name, domains, structure):
    """name = string that specifies a unique name for this complex
      domains = a list of integers specifying the domains (in order) 
        in the complex. -1 is used to indicate a strand break
      structure = a list of integers indicating the pairing of domains 
        in the complex. -1 is used to indicate unpaired. Strand breaks are
        currently shown as -2"""
    self.name = name
    self.domains = domains
    self.structure = structure
    self.update_available_domains()
    self.pretty_domains = ""
    self.dot_paren = ""

    # Used for updating finding strongly connected components in the fast
    # self-interaction subgraph
    self.outward_edges = []

  def __str__(self):
   if self.pretty_domains == "":
      return self.name + "\n" + str(self.domains) + "\n" + str(self.structure)
   else:
     return self.name + ":\n" + self.pretty_domains + "\n" + self.dot_paren

  def __repr__(self):
   if self.pretty_domains == "":
      return "Complex(" + self.name + ", " + str(self.domains) + ", " \
		    + str(self.structure) + ")"
   else:
     return "Complex(" + '"' + self.name + '"'+ ", " + \
       '"' + self.pretty_domains + '"'+ ", " + \
       '"' + self.dot_paren + '"' + ")"


  def update_available_domains(self):
    # find all the domains on external loops.
    # There has to be a more efficient way to do this....can't think of it
    # right now without spending some time to think about it. This'll work for
    # now.

    # One quick optimization, can skip the descent into the loops.
    #   Just need to set inner_index = struc[inner_index] if it is the start of
    #   a loop...
    doms = self.domains
    struc =self.structure + [-2]
    flag = False
    self.available_domains = []

    nStrands = 1
    for domain in doms:
      if domain == 0:
        nStrands += 1

    self.nStrands = nStrands

    for index in range(0,len(doms)):
      available = False
      if struc[index] == -1:
        depth = 0
        inner_index = index - 1
        exception_index = 0
        while inner_index >= 0 and inner_index < index \
            and not available:
          exception_index += 1
          if exception_index > len(struc):
            raise Exception("Update 1 Failed")
          if struc[inner_index] == -1:
            inner_index -= 1
          elif struc[inner_index] == -2:
            available = True
          else:
            inner_index = struc[inner_index] - 1
        inner_index = index + 1
        exception_index = 0
        while inner_index < len(struc) and inner_index > index and \
            not available:
          exception_index += 1
          if exception_index > len(struc):
            raise Exception("Update 2 Failed")
          if struc[inner_index] == -1:
            inner_index += 1
          elif struc[inner_index] == -2:
            available = True
          else:
            inner_index = struc[inner_index] + 1
        

      if available:
        self.available_domains.append( (doms[index],index)) 
    self.available_domains.sort()

  def rotate_strands(self):
    index = 0
    struc = self.structure
    doms = self.domains
    while index < len(struc) and struc[index] != -2 :
      index += 1
    
    break_index = index
    if break_index == len(struc): # 1 strand
      return self
    s1 = struc[:index]
    s2 = struc[index+1:]
    res_struc = s2 + [-2] + s1
    d1 = doms[:index]
    d2 = doms[index+1:]
    res_doms = d2 + [0] + d1
    for index in range(0,len(res_struc)):
      if res_struc[index] < 0:
        pass
      elif res_struc[index] < break_index:
        res_struc[index] += len(s2) + 1
      else:
        res_struc[index] -= len(s1) + 1
    res_comp = Complex(self.name + str("-rot"), res_doms, res_struc)
    return res_comp

  def equivalent(self, other):
    flag = False
    if other.structure == self.structure \
      and other.domains == self.domains:
      flag = True

    return flag

    
  def make_canonical(self, domain_list):
    max_rotation = self
    current_rotation = self.rotate_strands()
    index = 0
    self.make_pretty(domain_list)
    current_rotation.make_pretty(domain_list)
    while (self.pretty_domains != current_rotation.pretty_domains or \
            self.dot_paren != current_rotation.dot_paren) :
      if str(current_rotation.pretty_domains) > \
              str(max_rotation.pretty_domains) \
          or (str(current_rotation.pretty_domains) == \
              str(max_rotation.pretty_domains) and \
            str(current_rotation.dot_paren) > str(max_rotation.dot_paren)) :
        max_rotation = current_rotation
      current_rotation = current_rotation.rotate_strands()
      current_rotation.make_pretty(domain_list)

    self.structure = max_rotation.structure
    self.domains = max_rotation.domains

  def num_strands(self):
    return self.nStrands

  def make_pretty(self, dom_list):
    d_list = []
    for domain in self.domains:
      index = abs(domain)
      dom_name = dom_list[index].name
      if domain < 0:
        dom_name += "*"
      d_list.append(dom_name)

    self.pretty_domains = ' '.join(d_list)

    dot_paren_list = []
    for index in range(0,len(self.structure)):
      if self.structure[index] == -1:
        dot_paren_list.append('.')
      elif self.structure[index] == -2:
        dot_paren_list.append('+')
      elif self.structure[index] < index:
        dot_paren_list.append(')')
      elif self.structure[index] > index:
        dot_paren_list.append('(')
      else:
        raise Exception( "Error in structure, bad domain pair")
    self.dot_paren = ''.join(dot_paren_list)


identifier_pattern = """[A-Za-z0-9_]*"""

name_pattern = identifier_pattern + "[\t ]*:"
domain_pattern = "(" + identifier_pattern + "\*?[\t ]*)*"
structure_pattern = "[\(\)\.]*"

identifier_re = re.compile(identifier_pattern)
name_re = re.compile(name_pattern)
domain_re = re.compile(domain_pattern)
structure_re = re.compile(structure_pattern)

sequence_pattern = """sequence .+ : \-?[0-9]+"""
sequence_re = re.compile(sequence_pattern)


UNPAIRED = -1

##################################
# Default domain length
default_length = 6
default_release_cutoff = 8

class Domain:
  def __init__(self,name,length):
    self.name = name
    self.length = length

class Enumerator:
  def __init__(self):
    self.domain_dict = {'+': 0}
    self.domain_list = [Domain('+',0)]
    self.E = []
    self.T = []
    self.S = []
    self.N = []
    self.F = []
    self.domain_set = set(self.domain_dict.keys())
    self.name_int = 0
    self.MAX_REACTION_COUNT = 1000
    self.MAX_COMPLEX_COUNT = 200
    self.ACTIVATE_ZIP_REACTIONS = False
    self.initial_complexes = []

  

  def internal_hybridizations(self,complex):
    """ This function finds and returns a list of all the interactions which
    can occur within the complex that will not result in a pseudoknotted 
    structure. Note that this throws away some of the structures because of the
    unpseudoknotted requirement. The format of the returned values is:
    [(index,pair),(index,pair)...], indicating all changes from complex in
    the reaction
    domain_id is the unique identifier for that domain type, the index is the
    location in the complex """

    # First find internal hybridizations 
    match_list = []
    outer_index = 0
    inner_index = 0
    struc = complex.structure
    doms = complex.domains
    reactions = []
    for outer_index in range(0,len(struc)):
      if struc[outer_index] != UNPAIRED:
        # Paired or strand break, can't make a new pair
        continue
      inner_index = outer_index + 1
      while inner_index > outer_index and inner_index < len(struc):
        if struc[inner_index] == -2:
          inner_index += 1
        elif struc[inner_index] != UNPAIRED:
          inner_index = struc[inner_index] + 1
        elif doms[inner_index] == -doms[outer_index]:
          # Unpaired by previous statements, will pair, and is unpseudoknotted
          # because of the way that we are iterating through
          reactions.append((complex,[(outer_index,inner_index), \
            (inner_index,outer_index)],"Bind 1-1"))

          inner_index += 1
        else:
          inner_index += 1
    if self.ACTIVATE_ZIP_REACTIONS:
      
      zip_reacs = []
      for index in range(0,len(reactions)):
        for inner_index in range(index + 1,len(reactions)):
          reac1 = reactions[index]
          reac2 = reactions[inner_index]
          pair_list1 = reac1[1]
          pair_list2 = reac2[1]

          concat_flag = False

          for pair1 in pair_list1:
            for pair2 in pair_list2:
              if pair1[0] + 1 == pair2[0] and pair1[1] - 1 == pair2[1]:
                concat_flag = True
              if pair1[0] - 1 == pair2[0] and pair1[1] + 1 == pair2[1]:
                concat_flag = True

          if concat_flag:
            new_list = reac1[1] + reac2[1]
            reactions[inner_index] = (complex,new_list,"Bind1-1")
            reactions[index] = None
            break

      for index in range(len(reactions)-1,-1,-1):
        if reactions[index] == None:
          del reactions[index]

    return reactions

    
  def pairwise_interactions(self, complex, other):
    """Finds all the pairwise interactions which can occur between the
      complex and other. Returns a list of lists of domain pairs created
      [[(index1,index2),(index2,index1),( ...) ], [(index...) ...] ...]"""
    complex_index = 0
    other_index = len(other.available_domains) - 1

    sdoms = complex.available_domains
    odoms = other.available_domains

    interaction_list = []
    while complex_index < len(complex.available_domains) and other_index >= 0:
      # is this the start of a matching set?
      if sdoms[complex_index][0] == - odoms[other_index][0]:
        temp_complex_index = complex_index
        while temp_complex_index < len(sdoms) and \
            sdoms[temp_complex_index][0] == sdoms[complex_index][0]:
          temp_other_index = other_index
          while temp_other_index >= 0 and \
              odoms[temp_other_index][0] == odoms[other_index][0]:
            interaction_list.append([(sdoms[temp_complex_index][1], \
                odoms[temp_other_index][1])])
            temp_other_index -= 1
          temp_complex_index += 1
        complex_index = temp_complex_index
        other_index = temp_other_index
        
      # if sdoms > odoms, then increment sdoms 
      elif sdoms[complex_index][0] < -odoms[other_index][0]:
        complex_index += 1
        
      #if sdoms < odoms, then decrement odoms 
      elif sdoms[complex_index][0] > -odoms[other_index][0]:
        other_index -= 1

    if self.ACTIVATE_ZIP_REACTIONS:
      # look for reactions which can be concatenated.
      for index in range(0,len(interaction_list)):
        for inner_index in range(index + 1,len(interaction_list)):
          reac1 = interaction_list[index]
          reac2 = interaction_list[inner_index]
          concat_flag = False
          for pair1 in reac1:
            for pair2 in reac2:
              if pair1[0] + 1 == pair2[0] and pair1[1] - 1 == pair2[1]:
                concat_flag = True
              if pair1[0] - 1 == pair2[0] and pair1[1] + 1 == pair2[1]:
                concat_flag = True
          if concat_flag:
            interaction_list[index] = None
            interaction_list[inner_index] = reac1 + reac2
            break
      for index in range(len(interaction_list) - 1,-1,-1):
        if interaction_list[index] == None:
          del interaction_list[index]
    return interaction_list

  def displacements(self, complex):
    """Returns a list of tuples of 2 pairs that are created by displacement 
        reactions. e.g.
        [(((-3,1),(3,3)),((-3,5),(None,-1)))]
        This indicates that there was a single displacement reaction. Note that
        the original pair ((-3,5),(3,3)) goes to a different pair and an 
        unpaired strand"""
    dis_list = []
    dom = complex.domains
    
    for ind_i in range(0,len(complex.structure)):
      if complex.structure[ind_i] < 0:
        continue

      # displace forward
      ind1 = ind_i + 1 #Displacing domain
      if not ind1 < len(complex.structure):
        continue
      ind2 = complex.structure[ind1] # Displaced strand from initiator
      ind3 = complex.structure[ind_i] - 1 # Template strand (pair to initiator)
      if ind3 < 0:
        continue
      ind4 = complex.structure[ind3]     # Displaced from template

      if ind2 == -2 or ind4 == -2:
        continue
      if ind2 == -1 and ind4 == -1:   # This would just be a hybridization
        continue

      if ind2 == -1:
        dom2 = None
      else:
        dom2 = dom[ind2]
      
      if ind4 == -1:
        dom4 = None
      else:
        dom4 = dom[ind4]

      dom1 = dom[ind1]
      dom3 = dom[ind3]
      
      # Can they create a pair that doesn't already exist?
      if dom1 == -dom3 and ind2 != ind3:
        pair_tuples =  [(ind1,ind3),(ind3,ind1),(ind2,ind4),(ind4,ind2)]

        changes = []
        fourway = True
        # some of these checks aren't necessary any more, should be removed
        # left them in bc I'm not confident if they will be necessary when
        # this is changed
        for tuple in pair_tuples:
          if tuple[0] != -1:
            changes.append(tuple)
          else:
            fourway = False

        if fourway:
          reaction_type = "Four-way Displacement"
        else:
          reaction_type = "Three-way Displacement"

        reaction = (complex,changes,reaction_type)

        dis_list.append(reaction)
    return dis_list


  def releases(self,complex):
    """Format is [(ind1,pair1),(ind2,pair2)]. Both of pair1&2 are set to 
        unpaired to fulfill the reaction"""
    release_list = []

    struc = complex.structure
    dom = complex.domains


    for ind_i in range(0,len(struc)):
      flag = True

      helix_startA = ind_i 
      helix_startB = struc[ind_i]
      if helix_startB < helix_startA: 
        # strand break or unpaired, no release
        # or this has already been counted
        continue
      if helix_startA - 1 >= 0 and helix_startB + 1 < len(struc) and \
          helix_startA - 1 == struc[helix_startB + 1]:
        # Then this isn't the start of the helix, it has already been considered
        # ignore it
        continue

      helix_endA = helix_startA 
      helix_endB = helix_startB

      stop_flag = False

      helix_length = self.domain_list[abs(dom[helix_endA])].length

      while not stop_flag:
        helix_endA += 1
        helix_endB -= 1
        if helix_endA != struc[helix_endB]:
          stop_flag = True
        else:
          helix_length += self.domain_list[abs(dom[helix_endA])].length
          
      if helix_length < default_release_cutoff:
        change_list = []
        for i in range(helix_startA,helix_endA):
          change_list.append( (i,UNPAIRED))
        for j in range(helix_startB,helix_endB,-1):
          change_list.append( (j,UNPAIRED))
        release_list.append( \
          (complex,change_list,"Open"))
    return release_list

    
  def process_str_complex(self, name, str_domains, dot_paren_structure):
    cur_domain_set = set([s.rstrip("*") for s in str_domains])
    self.domain_set = set(self.domain_dict.keys())
    unassigned_domains = cur_domain_set - self.domain_set
    for domain in unassigned_domains:
      self.domain_dict[domain] = len(self.domain_list)
      self.domain_dict[domain + "*"] = -len(self.domain_list)
      self.domain_list.append(Domain(domain,default_length))
    self.domain_set = set(self.domain_dict.keys())
    int_domains = [self.domain_dict[s] for s in str_domains]

    int_struc = []

    struc_stack = []
    for index in range(0,len(dot_paren_structure)):
      int_struc.append(-10)
      char = dot_paren_structure[index]
      if char == '(':
        struc_stack.append(index) 
      elif char == ')':
        match_index = struc_stack.pop()
        if int_domains[index] != -int_domains[match_index] \
            or int_domains[index] == 0:
          raise Exception(str(int_domains[index]) + \
            " cannot pair with " + str(int_domains[match_index]))
        int_struc[match_index] = index
        int_struc[index] = match_index
      elif char == '.':
        int_struc[index] = -1
      elif char == '+':
        int_struc[index] = -2
      else:
        raise Exception("Invalid character in the structure")
    self.initial_complexes.append(Complex(name, int_domains, int_struc))

  def self_hybrid_complexes(self, reactions):
    complex = reactions[0]
    reacs = reactions[1]

    compList = []

    structure = complex.structure
    doms = complex.domains
    for reac in reacs:
      struc = copy.copy(structure)
      struc[reac[0][1]] = reac[1][1]
      struc[reac[1][1]] = reac[0][1]
      comp = Complex(str(self.name_int),doms,struc)
      self.name_int += 1
      compList.append(comp)

    return compList

  def find_insertion_index(self,break_index, structure):
    stop_flag = False
    index = break_index + 1
    insertion_index = -1
    found_flag = False


    while not stop_flag:
      if index >= len(structure):
        stop_flag = True
      elif structure[index] == -2:
        insertion_index = index
        found_flag = True
        stop_flag = True
      elif structure[index] == -1:
        index += 1
      elif structure[index] > index:
        index = structure[index] + 1
      else:
        stop_flag = True

    stop_flag = False
    index = break_index - 1
    while not found_flag and not stop_flag:
      if index < 0:
        stop_flag = True
        insertion_index = -1
        found_flag = True
      elif structure[index] == -2:
        insertion_index = index
        found_flag = True
        stop_flag = True
      elif structure[index] == -1:
        index -= 1
      elif structure[index] < index:
        index = structure[index] - 1
      else:
        stop_flag = True
    
    if not found_flag:
      raise Exception("Invalid call to find_insertion_index()")
      
    return insertion_index


  def cross_hybrid_complexes(self, reactions):
    comp1 = reactions[0]
    comp2 = reactions[1]
    reacs = reactions[2]
    comp_list = []
    # Find the insertion point (the insertion domain has to be on an exterior
    # loop, means that there is one and only one place where the insertion
    # can actually occur)
    for reac in reacs:
      copy_comp1 = copy.copy(comp1)
      copy_comp2 = copy.copy(comp2)
      comp1_pair_index = reac[0][0]
      comp2_pair_index = reac[0][1]
      comp1_dom = copy_comp1.domains
      comp2_dom = copy_comp2.domains
      comp1_structure = copy_comp1.structure
      comp2_structure = copy_comp2.structure
      # Find the insertion location in complex 1

      insertion_index = \
        self.find_insertion_index(comp1_pair_index, comp1_structure)
      insertion_index2 = \
        self.find_insertion_index(comp2_pair_index, comp2_structure)

      # The domains remain unchanged, except for the addition of a strand break
      # which was previously at the end of the second complex (and so unmarked)
      d2 = comp2_dom[insertion_index2:]
      d3 = [0] + comp2_dom[0:insertion_index2]
      if insertion_index > 0:
        d1 = comp1_dom[0:insertion_index]
        d4 = comp1_dom[insertion_index:]
        s1 = copy.copy(comp1_structure[0:insertion_index])
        s4 = copy.copy(comp1_structure[insertion_index:])
      else:
        d1 = comp1_dom
        d4 = []
        s1 = copy.copy(comp1_structure)
        s4 = []
      if insertion_index2 > 0:
        d2 = [0] + comp2_dom[insertion_index2 + 1:]
        d3 = [0] + comp2_dom[0:insertion_index2]
        s2 = copy.copy([-2] + comp2_structure[insertion_index2 + 1:])
        s3 = copy.copy([-2] + comp2_structure[0:insertion_index2])
      else:
        d2 = []
        d3 = [0] + comp2_dom
        s2 = []
        s3 = copy.copy([-2] + comp2_structure)


      s1_offset = 0
      s2_offset = len(s1) - len(s3) + 1 
      s3_offset = len(s1) + len(s2) + 1 # adding 1 bc of extra '+' sign in front
      s4_offset = len(s2) + len(s3)

      # e.g. .(..+.) b a b c + c a* pairing with .(..+.) c a b* c + c a*
      # pairing of the b's
      # insertion_index = 3
      # insertion_index2 = 4
      # orig struc1 = [-1,5,-1,-1,-2,-1,0]
      # orig struc2 = [-1,6,-1,-1,-2,-1,1]
  
      # s1 = [-1,6,-1,-1,-2,-1,1]
      # s2 = [-2,-1,1]
      # s3 = [-2,-1,6,-1,-1]
      # s4 = []

      # --> .(..+.)+.(+.)..
      # --> ((..+.)+.(+.)).
      # --> 13,6,-1,-1,-2,-1,1,-2,-1,12,-2,-1,9,0,-1


      # s1 = [-1,6,-1,-1]
      # s2 = [-2,-1,1] --> []
      # s3 = [-2,-1,6,-1,-1] --> []
      # s4 = [-2,-1,1] --> []
      # --> (..+.(+.)..+.)
      # --> -1,14,-1,-1,-2,-1,9,-2,-1,6,-1,-1,-2,-1,1
      #     b  a   b  c  +  c a* +  c a  b* c  +  c a*  
      # --> -1,14,10,-1,-2,-1,9,-2,-1,6, 2,-1,-2,-1,1

      # Remap the structure integers to take into account the fact that 
      # more domains were added

      for index in range(0,len(s1)):
        if s1[index] < 0:
          pass
        elif insertion_index < 0:
          pass
        elif s1[index] >= insertion_index:
          s1[index] += s4_offset

      for index in range(0,len(s2)):
        if s2[index] < 0:
          pass
        elif s2[index] >= insertion_index2 and insertion_index2 > 0:
          s2[index] += s2_offset
        elif s2[index] >= 0:
          s2[index] += s3_offset

      for index in range(0,len(s3)):
        if s3[index] < 0:
          pass
        elif s3[index] >= insertion_index2 and insertion_index2 > 0:
          s3[index] += s2_offset
        elif s3[index] >= 0:
          s3[index] += s3_offset

      for index in range(0,len(s4)):
        if s4[index] < 0:
          pass
        elif s4[index] >= insertion_index:
          s4[index] += s4_offset


      new_structure = s1+s2+s3+s4
      # Remap the reaction to the new indices
      for pair in reac:
        comp1_pair_index = pair[0]
        comp2_pair_index = pair[1]
        if comp1_pair_index >= insertion_index and insertion_index >= 0:
          comp1_pair_index += s4_offset
      
        if comp2_pair_index >= insertion_index2 and insertion_index2 >= 0:
          comp2_pair_index += s2_offset
        elif comp2_pair_index >= 0:
          comp2_pair_index += s3_offset
        else:
          print "Error, invalid reaction present in ch Reactions"

        new_structure[comp1_pair_index] = comp2_pair_index
        new_structure[comp2_pair_index] = comp1_pair_index
      
      new_domains = d1+d2+d3+d4


      # Create the complex
      new_complex = Complex(str(self.name_int), new_domains, new_structure)
      self.name_int += 1
      comp_list.append(([comp1,comp2],[new_complex],"Bind 2-1"))
    return comp_list

  def complex_split(self,comp, split_start, split_end):
    struc = comp.structure
    struce =[-2] + comp.structure 
    dom = comp.domains

    if split_start == -1:
      str1 = []
      dom1 = []
    else:
      str1 = struc[0:split_start]
      dom1 = dom[0:split_start] 

    str2 = struc[split_start + 1:split_end]
    dom2 = dom[split_start + 1:split_end]

    str3 = struc[split_end + 1:]
    dom3 = dom[split_end + 1:]
    if struce[split_start + 1] != -2 or struce[split_end + 1] != -2:
      raise Exception("Invalid split call, split_start or split end not at break")

    domains1 = dom1 + [0] + dom3 
    domains2 = dom2 

    if split_start == -1:
      str2_offset = 0
      str3_offset = -len(str2) - 1
    else:
      str2_offset = - len(str1) - 1
      str3_offset = - len(str2) - 1


    for index in range(0,len(str1)):
      if str1[index] >=split_end:
        str1[index] += str3_offset
      elif str1[index] >= split_start:
        raise Exception("Invalid split call, the strands are connected 1\n")
    for index in range(0,len(str2)):
      if str2[index] < 0:
        continue
      elif str2[index] < split_end and str2[index] > split_start:
        str2[index] += str2_offset
      else:
        raise Exception("Invalid split call, the strands are connected 2\n")
    for index in range(0,len(str3)):
      if str3[index] >= split_end:
        str3[index] += str3_offset
      elif str3[index] > split_start:
        raise Exception("Invalid split call, the strands are connected 3\n")

    structure1 = str1 + [-2] + str3
    structure2 = str2
    

    if structure1[0] == -2:
      del structure1[0]
      del domains1[0]
    if structure2[0] == -2:
      del structure2[0]
      del domains2[0]
    if structure1[-1] == -2:
      del structure1[-1]
      del domains1[-1]
    if structure2[-1] == -2:
      del structure2[-1]
      del domains2[-1]

    comp1 = Complex(str(self.name_int),domains1,structure1)
    self.name_int += 1
    comp2 = Complex(str(self.name_int),domains2,structure2)
    self.name_int += 1
    return [comp1,comp2]


  def find_detached_complexes(self, comp):
    """ Returns 2 complexes if the complex can be broken up. 
        Returns the complex in a list [complex] otherwise """
    struc = [-2] + comp.structure
    comp_list = []
    flag = False
    for outer_index in range(0,len(struc)):
      if struc[outer_index] == -2: # Strand break
        inner_index = outer_index - 1
        while inner_index >= 0 and inner_index < outer_index:
          if struc[inner_index] == -1:
            inner_index -= 1
          elif struc[inner_index] == -2:
            # Split the complex at outer_index and inner_index
            split_start = min(inner_index,outer_index) - 1 # bc of added -2
            split_end = max(inner_index, outer_index) - 1
            comp_list = self.complex_split(comp,split_start, split_end)
            flag = True
            break
          elif struc[inner_index] > inner_index and \
              struc[inner_index] < outer_index:
            inner_index -= 1
          else:
            inner_index = struc[inner_index]
        inner_index = outer_index + 1
        while inner_index < len(struc) and inner_index > outer_index and\
              (not flag):
          if struc[inner_index] == -1:
            inner_index += 1
          elif struc[inner_index] == -2:
            split_start = min(inner_index, outer_index) - 1
            split_end = max(inner_index, outer_index) - 1
            comp_list = self.complex_split(comp,split_start, split_end)
            flag = True
            break
          elif struc[inner_index] < inner_index and \
              struc[inner_index] > outer_index:
            inner_index += 1
          else:
            inner_index = struc[inner_index] + 2
            
          
      if flag:
        break

    if len(comp_list) == 0:
      comp_list = [comp]
    return comp_list


  def displacement_complexes(self, reactions):
    complex = reactions[0]
    reacs = reactions[1]
    
    results = []
    for reac in reacs:
      cur_comp = copy.deepcopy(complex)
      struc = cur_comp.structure
      cur_comp.name = str(self.name_int)
      self.name_int += 1
      four_way = False
      for index in [0,1]:
        if reac[index][0][1] < 0 and reac[index][1][1] >= 0:
          struc[reac[index][1][1]] = -1
          four_way = False
        if reac[index][1][1] < 0 and reac[index][0][1] >= 0:
          struc[reac[index][0][1]] = -1
          four_way = False
        if reac[index][1][1] >= 0 and reac[index][0][1] >= 0:
          struc[reac[index][0][1]] = reac[index][1][1]
          struc[reac[index][1][1]] = reac[index][0][1]
          four_way = True
      
      complex_list = self.find_detached_complexes(cur_comp)
      if complex_list == []:
        cur_comp.name = str(self.name_int)
        self.name_int += 1
        cur_comp.update_available_domains()
        complex_list = [cur_comp]
        if four_way:
          results.append(([complex],complex_list,"SF"))
        else:
          results.append(([complex],complex_list,"SD"))
      else:
        if four_way:
          results.append(([complex],complex_list,"CF"))
        else:
          results.append(([complex],complex_list,"CD"))
      
    return results      


  def create_fast_complexes(self, index_reactions):
    # on the input the format is 
    # (reactant complex,[(index,pair),...],"TYPE_INDICATOR")
    # output format ([reactant complex],[output complex,..],"TYPE_INDICATOR")

    reactions = []
    
    for ireaction in index_reactions:
      input_complex = ireaction[0]
      reaction_changes = ireaction[1]
      reaction_type = ireaction[2]

      output_complex = copy.deepcopy(input_complex)
      for change in reaction_changes:
        changed_index = change[0]
        changed_value = change[1]
        output_complex.structure[changed_index] = changed_value
        output_complex.name = str(self.name_int)
        self.name_int += 1

      
      output_list = self.find_detached_complexes(output_complex)

      for output in output_list:
        output.make_canonical(self.domain_list)
        output.update_available_domains()

      reactions.append(([input_complex],output_list,reaction_type))
    return reactions
      
  def parse_input_file(self, infilename):
    infile = open(infilename, "r")
    line = infile.readline()
    complex_list = []
    while line != "":
      
      if sequence_re.match(line):
        # Parse domain specification
        split_line = line.split()
        domain_name = split_line[1]
        if not domain_name.isalnum():
          print "domain name " + domain_name + " is not alphanumeric, ignoring"
        domain_length_str = split_line[-1]
        try:
          domain_length = int(domain_length_str)
        except ValueError:
          print "Invalid sequence length specification"
          print "Ignoring constraint"
          domain_length = default_length
        if domain_name in self.domain_dict:
          domain_int = self.domain_dict[domain_name]
          self.domain_list[domain_int].length = domain_length
        else:
          self.domain_dict[domain_name] = len(self.domain_list)
          self.domain_dict[domain_name + "*"] = -len(self.domain_list)
          self.domain_list.append(Domain(domain_name,domain_length))

      if name_re.match(line):
        name = identifier_re.match(line).group(0)
        line2 = infile.readline()
        if not domain_re.match(line2):
          raise Exception("Invalid domain specification")
        line3 = infile.readline()
        if not structure_re.match(line3):
          raise Exception("Invalid structure specification")
        str_domains = line2.split()
        str_structure = ''.join(line3.strip().split())
        self.process_str_complex(name,str_domains,str_structure)
        

      line = infile.readline()
    
    return len(self.initial_complexes)

  def make_slow_reactions(self, source):
    # This function finds the next set of bimolecular reactions that should be
    # considered and any other slow reactions. In the future it may be
    # beneficial to have a more general rate threshold which determines
    # whether the reaction is actually fast or slow.
    cross_interaction_list = []

    cross_interactions = self.pairwise_interactions(source,source)
    if len(cross_interactions) > 0:
      cross_interaction_list.append((source,source,cross_interactions))

    for complex in self.E:
      cross_interactions = self.pairwise_interactions(source, complex)
      if len(cross_interactions) > 0:
        cross_interaction_list.append((source,complex,cross_interactions))

    reactions = []
    for reac in cross_interaction_list:
      reactions.extend(self.cross_hybrid_complexes(reac))

    return reactions

      

  def make_fast_reactions(self, source):
    # This function currently makes the reactions which can occur from the 
    # source complex through only self-interactions, which are assumed to
    # be fast reactions. source is a complex.
    
    self_interaction_list = []
    
    self_interaction_list.extend(self.internal_hybridizations(source))

    self_interaction_list.extend(self.displacements(source))

    self_interaction_list.extend(self.releases(source))

    fast_reactions = self.create_fast_complexes(self_interaction_list)

    return fast_reactions
    

  def prune_complex(self, complex):
    return complex.num_strands() > 4

    # Note: this may include complexes which are already present. These will
    #     be needed to generate any edges leading into these complexes, but 
    #     will need to be removed before adding them all to the list of
    #     structures.

  def get_products(self, reactions):
  # get_products: get the products and check them against currently enumerated
  #               complexes/existing but not enumerated complexes
  #               will change the reference from the newly enumerated complex
  #               to the duplicate in the more-enumerated set.
  #   if product is in E,S,T: adjust pointer, move on, all interactions are
  #               enumerated already. can be simplified by noting that I
  #               already explored all the strongly connected components for
  #               these complexes
  #   if product is in N, then it is in the current neighborhood, so change the
  #               pointer and move on (not returned)
  #   if product is in F, then it is also in the current neighborhood, change
  #               the pointer and move on (not returned)
  #   if the product is in B, then we need to enumerate self-interactions now,
  #               remove it from B, change the pointer, and return the complex.
  #               (need to change the pointer so that the reactions are all
  #               accurate)

    res_products = []
    for reaction in reactions:
      products = reaction[1]
      
      for index in range(0,len(products)):
        products[index].make_canonical(self.domain_list)
        products[index].update_available_domains()
        enumerated = False
        for complex in self.E + self.S + self.T + self.N + self.F:
          if products[index].equivalent(complex):
            enumerated = True
            products[index] = complex
            break
        for complex in self.B:
          if products[index].equivalent(complex):
            products[index] = complex
            break
        for complex in res_products:
          if products[index].equivalent(complex):
            products[index] = complex

        if not enumerated:
          res_products.append(products[index])
    return res_products

  def tarjans(self, node):
    
    node.index = self.tarjan_index
    node.lowlink = self.tarjan_index
    self.tarjan_index += 1
    self.tarjan_stack.append(node)
    node.onStack = True
    
    for next in node.outward_edges:
      if next.index == -1:
        self.tarjans(next)
        node.lowlink = min(node.lowlink,next.lowlink)
      elif next.onStack:
        node.lowlink = min(node.lowlink, next.lowlink)
    if node.lowlink == node.index:
      stop_flag = False
      scc = []
      while stop_flag == False:
        nextNode = self.tarjan_stack.pop()
        nextNode.onStack = False
        scc.append(nextNode)
        if nextNode == node:
          stop_flag = True

      self.SCC_stack.append(scc)


      
    

  def strong_components(self, N, N_reacs):
    # construct the graph from the reactions, use the complexes as the nodes
    # of a graph, this could be useful to save until later so that we know
    # all of the strongly connected components downstream of a given node.
    # Isn't necessary to save the graph, but will for now.
    # This algorithm is fairly standard, this one was copied largely from
    # wikipedia's article on Tarjan's algorithm

    ################
    # Initialize the graph
    self.tarjan_index = 0
    self.tarjan_stack = []
    self.SCC_stack = []

    for node in N:
      node.outward_edges = []
      node.full_outward_edges = []
      node.index = -1
      node.lowlink = -1

    for reac in N_reacs:
      # assuming all fast reactions are unimolecular sources
      
      for product in reac[1]:
        prod_in_N = False
        for complex in N:
          if product == complex:
            prod_in_N = True
            break
        if prod_in_N:
          reac[0][0].outward_edges.append(product)

      reac[0][0].full_outward_edges.extend(reac[1])

    for node in N:
      if node.index == -1:
        self.tarjans(node)

    end_states = []
    for scc in self.SCC_stack:
      scc_products = []
      are_end_states = True
     
      for node in scc:
        for product in node.full_outward_edges:
          scc_products.append(product)
      for product in scc_products:
        prod_in_scc = False
        for complex in scc:
          if product.equivalent(complex):
            prod_in_scc = True
            break
        if not prod_in_scc:
          are_end_states = False
          break
      if are_end_states:
        end_states.extend(scc)


    return end_states
        
  
  def find_reactions(self):
    #####################################################################
    # Enumerates all states as described in documentation available from
    # Brian Wolfe or another person who takes over the project.
    # Returns E, the end-state complexes, T, the transition state complexes, and
    # the reactions which move between these complexes.
    # Will find all the reactions of an initialized Enumerator object
    #
    # E = enumerated end states. only cross-reactions with other end states
    #     have yet to be considered for these complexes. These states will
    #     remain here until the end of the program.
    # S = end states which have not yet had cross-reactions with E considered
    #     yet (but have all self-interactions enumerated already). After
    #     cross reactions with E are considered (and self bimolecular reactions)
    #     then this complex will be added to E.
    # T = transition states. Enumerated states which are not end states. 
    #     all self-reactions have been considered. No bimolecular reactions are
    #     possible. These states will remain here until program end.
    # N = current neighborhood. These are the states in the current fast-
    #     neighborhood that have had all self-interactions enumerated, but
    #     the neighborhood hasn't been characterized yet (these states will
    #     be partitioned into S (if they are strongly connected components) or
    #     T (if they aren't) after the neighborhood is completely enumerated.
    # F = fast states which have been discovered, but self interactions haven't
    #     been enumerated yet (in the current neighborhood). These will be 
    #     moved to N as soon as they are enumerated.
    # B = products of bimolecular reactions that have had no reactions 
    #     enumerated yet. These will be moved to F when their neighborhood is
    #     about to be considered.

    for complex in self.initial_complexes:
      complex.make_canonical(self.domain_list)
      complex.update_available_domains()

    self.E = [] # no complexes fully enumerated
    self.S = []
    self.T = []
    self.N = []
    self.F = []
    self.B = self.initial_complexes
    self.reactions = []

    while len(self.B) > 0:
      N_reacs = []
      fast_source = self.B.pop()
      self.F = [fast_source]
      while len(self.F) > 0:
        fast_element = self.F.pop()
        fast_reactions = self.make_fast_reactions(fast_element)
        fast_products = self.get_products(fast_reactions)
        self.F.extend(fast_products)
        N_reacs.extend(fast_reactions)
        self.N.append(fast_element)

      strong_components = self.strong_components(self.N,N_reacs)
      for complex in strong_components:
        # If it is in the strong_components then it is an end-state, so add
        # it to S
        self.N.remove(complex)
        self.S.append(complex)
        # Otherwise, add it to the transition states T.
      self.T.extend(self.N)
      self.reactions.extend(N_reacs)
      self.N = []



    while len(self.S) > 0:
      # print "E"
      # print "-----"
      # print self.E
      # print "-----------------------------"
      # print "S"
      # print "----"
      # print self.S
      # print "T"
      # print "----"
      # print self.T
      element = self.S.pop()
      slow_reacs = self.make_slow_reactions(element) # bimolecular reaction enum
      self.E.append(element)
      if len(self.B) != 0:
        print "len(B) = " + str(len(B))
        print "Exiting..."
        raise Exception("bug in outer loop")
      self.B = self.get_products(slow_reacs) 
                                        # add resulting products to B (or note
                                        # that they are already enumerated and
                                        # change references in slow_reacs and
                                        # don't return them
      self.reactions.extend(slow_reacs)
      while len(self.B) > 0:
        if len(self.E) + len(self.T) + len(self.S) > self.MAX_COMPLEX_COUNT:
          raise Exception("MAX_COMPLEX_COUNT exceeded")
        if len(self.reactions) > self.MAX_REACTION_COUNT:
          raise Exception("MAX_REACTION_COUNT exceeded")
        N_reacs = []
        fast_source = self.B.pop()
        self.F = [fast_source]
        while len(self.F) > 0:
          fast_element = self.F.pop() # current element to enumerate
          # make all the fast reactions
          fast_reactions = self.make_fast_reactions(fast_element)
          # get the products from the reactions (and resolve any duplicates 
          # that were already enumerated in a past reaction)
          fast_products = self.get_products(fast_reactions)
          # add these back to the source set so the whole neighborhood gets
          # the fast reactions enumerated.
          self.F.extend(fast_products)
          # save the reactions so that we can construct the graph for finding
          # strongly connected components
          N_reacs.extend(fast_reactions)
          # fast elements were
          self.N.append(fast_element)
        # Now we have the neighborhood from fast_source, just need to find the
        # strongly-connected components. This function uses the graph specified
        # by the reactions and uses Tarjan's algorithm for identifying the
        # strongly connected components. It returns a list of nodes
        # in component that have no fast outgoing edges.
        # print "\n\n"
        # print "N"
        # print self.N
        # print "Reactions"
        # print N_reacs
        # print "--------\n\n"
        strong_components = self.strong_components(self.N,N_reacs)
        for complex in strong_components:
          # If it is in the strong_components then it is an end-state, so add
          # it to S
          self.N.remove(complex)
          self.S.append(complex)
          # Otherwise, add it to the transition states T.
        self.T.extend(self.N)
        self.reactions.extend(N_reacs)
        self.N = []


    ############################################
    # Translate the results into dot-parens and the given domain names
    # Complex objects
    
    for complex in self.E + self.T:
      complex.make_pretty(self.domain_list)

    return (self.E,self.T, self.reactions)
