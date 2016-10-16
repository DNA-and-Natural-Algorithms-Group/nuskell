
class iupac_translator(object):
  # A class to handle constraints in IUPAC notation.
  T = 'T'
  def __init__(self, molecule='DNA'):
    assert molecule == 'DNA' or molecule == 'RNA'
    self.T = 'T' if self.molecule == 'DNA' else 'U'

  def iupac(self, nuc):
    T = self.T if self.T else iupac_translator.T
    iupac_dict = {
      'A' : 'A',
      'C' : 'C',
      'G' : 'G',
       T  :  T,
      'R' : 'AG',   # purine
      'Y' : 'C'+T,  # pyrimidine
      'S' : 'CG',   # strong
      'M' : 'AC',
      'W' : 'A'+T,  # weak
      'K' : 'G'+T,
      'V' : 'ACG',  # not T
      'H' : 'AC'+T, # not G
      'D' : 'AG'+T, # not C
      'B' : 'CG'+T, # not A
      'N' : 'ACG'+T }
    return iupac_dict[nuc]

  def iupac_union(self, nucs):
    # Return the maximal common constraint
    u = 'N'
    for n in nucs :
      u = self.bin_iupac(self.iupac_bin(u) & self.iupac_bin(n))
    return u

  def iupac_neighbor(self, nuc):
    T = self.T if self.T else iupac_translator.T
    neighbor_dict = {  # ACGT => ACGT
      'A' :  T,        # 1000 => 0001 
      'C' : 'G',       # 0100 => 0010 
      'G' : 'Y',       # 0010 => 0101 
       T  : 'R',       # 0001 => 1010 
      'T' : 'R',       # 0001 => 1010 
      'R' : 'Y',       # 1010 => 0101 
      'Y' : 'R',       # 0101 => 1010 
      'S' : 'B',       # 0110 => 0111 
      'M' : 'K',       # 1100 => 0011 
      'W' : 'D',       # 1001 => 1011 
      'K' : 'N',       # 0011 => 1111 
      'V' : 'B',       # 1110 => 0111 
      'H' : 'D',       # 1101 => 1011 
      'D' : 'N',       # 1011 => 1111 
      'B' : 'N',       # 0111 => 1111 
      'N' : 'N'}       # 1111 => 1111 
    return neighbor_dict[nuc]

  def iupac_bin(self, nuc):
    T = self.T if self.T else iupac_translator.T
    iupac_bin_dict = { # ACGT
      'A' : 8,      # 1000,
      'C' : 4,      # 0100,
      'G' : 2,      # 0010,
       T  : 1,      # 0001,
      'R' : 10,     # 1010,  # purine
      'Y' : 5,      # 0101,  # pyrimidine
      'S' : 6,      # 0110, 
      'M' : 12,     # 1100, 
      'W' : 9,      # 1001, 
      'K' : 3,      # 0011, 
      'V' : 14,     # 1110,  # not T
      'H' : 13,     # 1101,  # not G
      'D' : 11,     # 1011,  # not C
      'B' : 7,      # 0111,  # not A
      'N' : 15}     # 1111,
    return iupac_bin_dict[nuc]
 
  def bin_iupac(self, nuc):
    T = self.T if self.T else iupac_translator.T
    bin_iupac_dict = [ # ACGT
      '',           # 0000  0 
       T,           # 0001  1 
      'G',          # 0010  2
      'K',          # 0011  3
      'C',          # 0100  4
      'Y',          # 0101  5
      'S',          # 0110  6
      'B',          # 0111  7
      'A',          # 1000  8
      'W',          # 1001  9 
      'R',          # 1010 10
      'D',          # 1011 11
      'M',          # 1100 12
      'H',          # 1101 13
      'V',          # 1110 14
      'N']          # 1111 15
    return bin_iupac_dict[nuc]

class Domain(iupac_translator):
  """The Nuskell Domain class.

  SeqDomain = SeqDomain | [Domains]

  A domain is a sequence of cnsequtive nucleotides, or a sequence of
  consequtive domains. A domain can be present in multiple DNA complexes, and
  its secondary structure can change dependent on the context. However, a
  single domain must always correspond to a single structure. Every domain has
  a unique descriptor (domain-ID) that allows for accessing/modifying every
  occurence in the System.
  
  NOTE: For domains on the nucleotide level, or first-level abstractions
  commonly used in strand displacement circuits, this is must (?) be dot-parens
  notation. For domains of domains, one may invent a new character mapping, e.g. 
  #"d1 = d1a d1b d1c" "h = ( . )"

  Every 'Domain' can have exactly one 'ComplementDomain', but it can be
  complementary to multiple Domains.

  Global:
    id_counter: stores the last automatically assigend domain-ID in a System.

  """

  id_counter = 0
 
  def __init__(self, name='', constraints=[], subdomains=[], complements=set()):
    """Initialization of the Domain Object.

    Arguments:
      name (str)             -- Name of this domain. If not given,
                                 an automatic one is generated.
      constraints (str)      -- Sequence constraints on this domain.
                                 Specifiers may be any of
                                 A,T,C,G,R,Y,W,S,M,K,B,V,D,H,N.
                                 Either constraints or subdomains
                                 must be specified.
      subdomains ([Domains]) -- List of component Domain objects. Either
                                 constraints or subdomains must be specified.
    """
    # Assign id
    self.id = Domain.id_counter
    Domain.id_counter += 1
 
    # Assign name
    if name :
      self._name = name
    else :
      self._name = 'domain_{0}'.format(self.id)
    
    # Assign constraints or subdomain list
    if constraints and not subdomains :
      assert all(isinstance(c, str) for c in constraints)
      self._constraints = constraints
      self._subdomains = None
    elif subdomains and not constraints :
      assert all(isinstance(d, Domain) for d in subdomains)
      self._subdomains = subdomains
      self._constraints = None
    else:
      raise ValueError("Must pass one of 'constraints' or 'subdomains' \
          keyword argument.")

    self._ComplementDomain = None

    # Initialize an empty set of complementary Domains. This is an experimental
    # feature that might come in handy for RNA. It will raise Warnings when
    # used at this point
    self._complements = set()
    if complements :
      raise Warning("Specifying complements this is an experimental feature.")
      assert all(isinstance(d, Domain) for d in complements)
      self.add_complements(complements)

  @property
  def length(self):
    return len(self.sequence)

  @property
  def base_length(self):
    return len(self.base_sequence)
    
  @property
  def sequence(self):
    if self._constraints :
      # This should also allow to overwrite the constraints
      return self._constraints
    else :
      return self._subdomains

  @property
  def base_sequence(self):
    """Breaks down the domain into non-composite domains."""
    if self._constraints :
      return self._constraints
    return list(''.join(map(lambda x: ''.join(x.base_sequence), self._subdomains)))

  def update_constraints(self, con):
    """ Unify new and old constraint """
    if len(con) != len(self._constraints) :
      raise ValueError("Length of new constraint != old constraint")
    for i in range(len(con)):
      con1 = con[i]
      con2 = self._constraints[i]
      con[i] = self.iupac_union([con1, con2])
      if not con[i] :
        raise ValueError("Constraints cannot be satisfied")
    return con

  def get_ComplementDomain(self, constraints=[], subdomains=[]):
    #TODO: This ignores the 
    if not self._ComplementDomain :
      if constraints :
        if len(constraints) != len(self._constraints) :
          raise ValueError("Length of constraint != complementary constraint")
        for i in range(len(constraints)):
          con1 = constraints[i]
          con2 = self.iupac_neighbor(self._constraints[i])
          constraints[i] = self.iupac_union([con1, con2])
          if not constraints[i] :
            raise ValueError("Constraints cannot be satisfied")
      self._ComplementDomain = ComplementDomain(self, 
          constraints=constraints, subdomains=subdomains)
    return self._ComplementDomain

  # def add_complements(self, other):
  #   # TODO: Check if constraints or subdomains are really complementary!
  #   raise Warning("add_complements: Experimental Feature.")
  #   for c in other :
  #     other.add_complements(set(self))
  #     self.complements.add(other)

  # def all_complements(self):
  #   # return all self.complements
  #   raise NotImplementedError

  # Built-in functions
  def __eq__(self, other):
    # Needs to have the same ID
    if isinstance(other, Domain): 
      if type(other) == ComplementDomain :
        # NOTE: ComplementDomain is an instance of Domain, removing this line
        # breaks a test case when the class is inherited elsewhere.
        return False
      return self.id == other.id
    else :
      return False
  def __ne__(self, other):
    return not (self == other)

  def __str__(self):
    return self._name

  def __invert__(self):
    """ Return a complement of this Domain. 

    This function will not automatically make a ComplementDomain, as it is
    unclear how constraints should be handled in this case. However,
    self._ComplementDomain can be initilized with get_ComplementDomain()
    """
    if self._complements:
      raise Warning("Multiple complements")
      return self._complements
    else :
      return self._ComplementDomain

class ComplementDomain(Domain):
  """
  Represents a complemented domain. Note that this is always
  defined in terms of an original domain and does not have the same
  data members, instead providing an interface to the complementary
  members.

  """
  def __init__(self, CompDomain, constraints=[], subdomains=[]):
    """Create the ComplementDomain for a given domain. A ComplementDomain can
    only be initialized with respect to a given Domain. There can only exist
    one ComplementDomain for every Domain. 
    """
    self.id   = CompDomain.id 
    self._name = CompDomain._name + '*'
    self._ComplementDomain = CompDomain

    # Assign constraints or subdomain list
    if constraints and not subdomains :
      self._constraints = constraints
    elif subdomains and not constraints :
      self._subdomains = self.update_subdomains(subdomains)
    else:
      raise ValueError("Must pass one of 'constraints' or 'subdomains' \
          keyword argument.")

  def __invert__(self):
    return self._ComplementDomain

  def __eq__(self, other):
    """ Returns True iff their complements are equal."""
    if isinstance(other, ComplementDomain): 
      return self._ComplementDomain.__eq__(other._ComplementDomain)
    else :
      return False
  def __ne__(self, other):
    """ Returns True iff they are not equal."""
    return not self.__eq__(other)

class NusDomain(Domain):
  def __init__(self, domaintag='', **kwargs):
    """ 
    Args: 
      domaintag: You can specify different types of domains for 
        your sequence designer: 'toehold', 'branch-migration', etc...
    """
    super(NusDomain, self).__init__(**kwargs)

    if domaintag == 'toehold':
      self._name = 't{}'.format(self.id)
    elif domaintag == 'wildcard' :
      self._name = '?{}'.format(self.id)
    else :
      self._name = 'd{}'.format(self.id)

class Complex(object):
  # Replace the Structure Object
  pass

class Reaction(object):
  pass

class TestTube(object):
  pass

class Solution(TestTube):
  # Replace the Nuskell Solution Object
  pass

