"""The nuskell objects are to a large ident identical with "DNAObjects" from
Joseph Schaeffer and Joseph Berleant. 

Currently, the following Objects are implemented:
  *IUPAC_translator: handling of nucleic acid constaints
  *Domain: A constrained subsequence of a molecule
  *ComplementDomain: A domain complementary to a Domain Object.
  *Complex: A sequence/structure pair.
  #*TestTube: A set of complexes
"""

def find(l, key):
  for i in range(len(l)):
    if l[i] == key:
      return i
  return None

def flatten(l) :
  if l == []:
    return l
  if isinstance(l[0], list):
    return flatten(l[0]) + flatten(l[1:])
  return l[:1] + flatten(l[1:])

class IUPAC_translator(object):
  # A class to handle constraints in IUPAC notation.
  T = 'T'
  def __init__(self, molecule='DNA'):
    assert molecule == 'DNA' or molecule == 'RNA'
    self.T = 'T' if self.molecule == 'DNA' else 'U'

  def iupac(self, nuc):
    T = self.T if self.T else IUPAC_translator.T
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
    T = self.T if self.T else IUPAC_translator.T
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
    T = self.T if self.T else IUPAC_translator.T
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
    T = self.T if self.T else IUPAC_translator.T
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

class Domain(IUPAC_translator):
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
 
  def __init__(self, constraints=[], subdomains=[], name='', prefix='d'):
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

    # Assign name
    # -----------
    # NOTE: The Domain function does not allow to specify names or prefixes to
    # end with digits, as proximal digits always are an automatically assigned
    # ID of a complex. A mixture of both naming modes is forbidden, as it must
    # be forbidden to initialize two domains with the same name. The
    # information of two domains with the same ID being different would be lost
    # when the domains are written into a file. 
    self.id = Domain.id_counter
    Domain.id_counter += 1

    if name :
      if name[-1]=='*' :
        raise ValueError('Invalid name for Domain Object!')
      if name[-1].isdigit() :
        raise ValueError('Domain name must not end with a digit!', name)
      self._name = name
    else :
      if prefix == '' :
        raise ValueError('Domain prefix must not be empty!')
      if prefix[-1].isdigit():
        raise ValueError('Domain prefix must not end with a digit!')
      self._name = prefix + str(self.id)

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
    #if complements :
    #  raise Warning("Specifying complements this is an experimental feature.")
    #  assert all(isinstance(d, Domain) for d in complements)
    #  self.add_complements(complements)

  @property
  def name(self):
    return self._name

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
    return list(''.join(
      map(lambda x: ''.join(x.base_sequence), self._subdomains)))

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
    else :
      self._ComplementDomain.update_constraints(constraints)
    return self._ComplementDomain

  @property
  def is_ComplementDomain(self):
    return self._name[-1:] == '*'

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
    # NOTE: this might actually be equivalent to Python's "is"
    if type(self) == type(other):
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

#class NusDomain(Domain):
#  def __init__(self, domaintag='', **kwargs):
#    """ 
#    Args: 
#      domaintag: You can specify different types of domains for 
#        your sequence designer: 'toehold', 'branch-migration', etc...
#    """
#    super(NusDomain, self).__init__(**kwargs)
#
#    if domaintag == 'toehold':
#      self._name = 't{}'.format(self.id)
#    elif domaintag == 'wildcard' :
#      self._name = '?{}'.format(self.id)
#    else :
#      self._name = 'd{}'.format(self.id)

#class Structure(object):
#  """The Structure of a DNA complex.
#
#  The structure is specified using a list of domains and a corresponding
#  dot-parens notation. 
#
#  Args:
#    domains (List[str]): names of domains
#    dotparens (List[str]): A list of dot-parens characters ".(.+)."
#    attr (Dict[name]=Domain): a mapping between names and Domain objects
#  """
#  def __init__(self, domains, dotparens, attr):
#    self.domains = domains
#    self.dotparens = dotparens
#    self.attributes = attr

class Complex(object):
  """A complex is a sequence & structure pair. 
  
  This implementation requires both sequence and structure to be of the same
  length, either on the domain-level or on the nucleotide-level.   
  """
  # TODO: is there something like a complex of complexes??
  id_counter = 0

  def __init__(self, sequence=[], structure=[], name='', prefix='cplx'):
    """
    Initialization:
    
    Keyword Arguments:
    name [type=str]                -- Name of the complex. An automatic name
                                      "complex_###" is assigned if this is not
                                      given.
    strands [type=list of strands] -- List of strands in this complex
    structure [type=str OR list]   -- Structure for this complex, in dot-paren
                                      notation or strand-list notation
                                      as used by Casey's enumerator. Note that
                                      only the latter representation can be
                                      used for pseudoknotted structures.
    """
    # Assign id
    self.id = Complex.id_counter
    Complex.id_counter += 1

 
    # Assign name
    if name :
      if name[-1]=='*' :
        raise ValueError('Invalid name for Complex Object!')
      if name[-1].isdigit() :
        raise ValueError('Complex name must not end with a digit!', name)
      self._name = name
    else :
      if prefix == '' :
        raise ValueError('Complex prefix must not be empty!')
      if prefix[-1].isdigit():
        raise ValueError('Complex prefix must not end with a digit!')
      self._name = prefix + str(self.id)

    #TODO: this might change!
    if sequence == [] :
      raise ValueError('Complex: requires Sequence and Structure Argument')

    if len(sequence) != len(structure) :
      raise ValueError("Complex: sequence and structure must have same length")
    self._sequence = sequence
    self._structure = structure

  @property
  def sequence(self):
    return self._sequence

  @property
  def lol_sequence(self):
    """ Returns sequence as a list of lists, rather than one flat list with the
    '+' separator.
    
    Example: 
      ['d1', 'd2', '+', 'd3'] ==> [['d1', 'd2'], ['d3']]
    """
    indices = [-1] + [i for i, x in enumerate(self._sequence) if x == "+"]
    indices.append(len(self._sequence))
    return [self._sequence[indices[i-1]+1: indices[i]] 
        for i in range(1, len(indices))]

  @property
  def nucleotide_sequence(self):
    """Return the complex sequence in form of a flat list of nucleotides. """
    lol = self.lol_sequence
    def my_base_sequence(seq) :
      if all(isinstance(d, Domain) for d in seq):
        return map(lambda x: x.base_sequence, seq)
      else :
        return seq
    return flatten(map(lambda x: my_base_sequence(x) + ['+'], lol))[:-1]

  @property
  def structure(self):
    return self._structure

  @property
  def lol_structure(self):
    # NOTE: It is not clear if this makes sense
    raise NotImplementedError

  @property
  def nucleotide_structure(self):
    # NOTE: It is not clear if this makes sense
    raise NotImplementedError

  @property
  def rotate_once(self):
    """Rotate the strands within the complex and return the updated object. """
    if "+" in self._sequence :
      p = find(self._sequence, "+")
      self._sequence = self._sequence[p + 1:] + ["+"] + self._sequence[:p]

      stack = []
      for i in range(p):
        if self._structure[i] == "(": stack.append(i)
        elif self._structure[i] == ")": stack.pop()
      for i in stack:
        self._structure[i] = ")"

      stack = []
      for i in reversed(range(p + 1, len(self._structure))):
        if self._structure[i] == ")": stack.append(i)
        elif self._structure[i] == "(": stack.pop()
      for i in stack:
        self._structure[i] = "("
      self._structure = self._structure[p + 1:] + ["+"] + self._structure[:p]
    return self

  @property
  def rotate(self):
    """Generator function yielding every rotation of the complex. """
    for i in range(len(self.lol_sequence)):
      yield self.rotate_once

  def __str__(self):
    return self._name

  def __eq__(self, other): 
    """ Two complexes are equal if they have the same coarse-graining in terms
    of domains and the same secondary structure. 
    NOTE: This function might change the strand-ordering of the complex!
    """
    if type(self) != type(other):
      return False
    if len(self.sequence) != len(other.sequence):
      return False
    if len(self.nucleotide_sequence) != len(other.nucleotide_sequence):
      return False
    if self._sequence == other.sequence and self.structure == other.structure:
      return True
    else :
      for r in self.rotate :
        if r.sequence == other.sequence and r.structure == other.structure:
          return True
      return False
  def __ne__(self, other):
    """ Returns True if the complexes are not equal. """
    return not self.__eq__(other)
 
class Reaction(object):
  pass

class TestTube(object):
  """ A graph representation of a test tube containing nucleic acids.

  The nodes of the graph are complexes, the edges define a neighborhood
  relation. 

  Node attributes:
    - free energy
    - concentration (mol/L, INF)

  Edge attributes:
    - transition rate
  
  """

  def __init__(self, domains=[], complexes=[]):
    assert all(isinstance(d, Domain) for d in domains)
    assert all(isinstance(c, Complex) for c in complexes)
    self._domains = dict()
    self._complexes = dict()

  @property
  def domains(self):
    return self._domains

  @property
  def complexes(self):
    return self._complexes

  #def load_pilfile(self, pilfile):
  #  from nuskell.parser import parse_pil_file
  #  [complexes] = parse_pil_file(pilfile)

  def add_complex(self, cplx):
    if cplx.name in self._complexes :
      return False
    else :
      self._complexes[cplx.name] = cplx
      return True

  def add_complex_by_name(self, name, sequence, structure):
    if name in self._complexes :
      return False
    else :
      self._complexes[name] = Complex(
          sequence=sequence, structure=structure, name=name)
      return True

  def add_domain(self, domain):
    # This method does not care about complementarity
    if domain.name in self._domains :
      return False
    else :
      self._domains[domain.name] = domain
      return True
 
  def add_domain_by_name(self, name, sequence):
    """Adds a domain with a particular Name to the TestTube. 
    
    If a domain with the same name has been added before, there will be no
    changes. Note that the Domain function does not allow to specify names that
    end with digits, as such names are reserved for auto-generated domains.
    """
    if name in self._domains and self._domains[name] :
      # Temporary placeholders for complement domains (None)
      return False
    elif name in self._domains :
      if name+'*' not in self._domains:
        raise ValueError("Complement Domain missing in the TestTube")
      self._domains[name] = self._domains[name+'*'].get_ComplementDomain(
          constraints=sequence)
    else :
      # Initialize new Domain
      if name[-1] == '*':
        tmp = Domain(constraints=list('N' * len(sequence)), name=name[:-1])
        self._domains[name[:-1]] = None
        self._domains[name] = tmp.get_ComplementDomain(constraints=sequence)
      else :
        self._domains[name] = Domain(constraints=sequence, name=name)
    return True
 
  def concentrations(self):
    pass
  
class Solution(TestTube):
  # Replace the Nuskell Solution Object
  pass

