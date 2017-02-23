# -*- coding: utf-8 -*-
#
# Copyright (c) 2016 Caltech. All rights reserved.
# Written by Stefan Badelt (badelt@caltech.edu)
#
# nuskell.objects: shared between different components of nuskell
#
# nuskell.objects are to some extent identical to "DNAObjects" coded by 
# Joseph Schaeffer and Joseph Berleant. 
#
# The following Objects are implemented:
#  *IUPAC_translator: handling of nucleic acid constraints
#  *Domain: A constrained subsequence of a molecule
#  *ComplementDomain: A domain complementary to a Domain Object.
#  *Complex: A sequence/structure pair.
#  *TestTube: A set of Complex objects and interface to
#     - read/write PIL files
#     - enumerate species

import networkx as nx
from nuskell.parser import parse_pil_file

def pair_table(ss, chars=['.']):
  """Return a secondary struture in form of pair table:

  :param ss: secondary structure in dot-bracket format
  :param base: choose between a pair-table with base 0 or 1
  :param chars: a list of characters that are ignored, default: ['.']

  :exmaple:
     base=0: ((..)). => [5,4,-1,-1,1,0,-1]
      i.e. start counting from 0, unpaired = -1
     base=1: ((..)). => [7,6,5,0,0,2,1,0]
      i.e. start counting from 1, unpaired = 0, pt[0]=len(ss)

  :return: a pair-table
  :return-type: list
  """
  stack=[];

  pt=[-1] * len(ss);

  for i, char in enumerate(ss):
    if (char == '('):
      stack.append(i);
    elif (char == ')'):
      try :
        j=stack.pop();
      except IndexError, e :
        raise RuntimeError("Too many closing brackets in secondary structure")
      pt[i]=j
      pt[j]=i
    elif (char == '+') :
      pt[i] = '+'
    elif (char not in set(chars)):
      raise ValueError("unexpected character in sequence: '" + char + "'")

  if stack != [] :
    raise RuntimeError("Too many opening brackets in secondary structure")
  return pt

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
  """The Nuskell Domain object.

  A domain is a sequence of consecutive nucleotides and/or domains. A domain
  can be present in multiple complexes, and its secondary structure can change
  dependent on the context. However, a single domain must always correspond to
  a single structure. Every domain has a unique descriptor (domain-ID) that
  allows for accessing/modifying every occurrence in the System.
  
  Every Domain can have exactly one ComplementDomain, but it can be
  complementary to multiple Domains. Hence, the concept of a ComplementDomain
  may or may not be used.

  Global Variables:
    id_counter: the next automatically assigned domain-ID in a system.
  """

  id_counter = 0
 
  def __init__(self, sequence=[], name='', prefix='d', complements=None):
    """Initialization of the Domain object.

    Arguments:
      sequence <list>         = A list of sequence constraints on this domain.
      name <optional: str>    = Name of this domain. If not specified, an
                                automatic name is generated.
      prefix <optional: str>  = A prefix for autmated naming of Domains. 
                                Default: 'd' for 'domain'

    Note: The Domain function does not allow names or prefixes ending with
    digits or '*'. Proximal digits are reseved to name domains automatically
    using their unique domain-IDs. A mixture of both naming modes is forbidden.
    Beware, that it is still possible to initialize two different domains with
    the same name, and that the information that these domain are different is
    lost when writing them to a file.
    """

    # Assign name
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

    # Assign domain sequence 
    if not sequence or type(sequence) != list :
      raise ValueError("Must pass a sequence of constraints or domain objects.")
    else :
      assert all(isinstance(s, (str, Domain)) for s in sequence)
      self._sequence = sequence

    self._ComplementDomain = None

    # Initialize an empty set of complementary Domains. This is an experimental
    # feature and an alternative to defining *one* ComplementDomain.
    if complements :
      raise NotImplementedError(
          "Specifying complements this is an experimental feature.")
      assert all(isinstance(d, Domain) for d in complements)
      #self.add_complements(complements)
    else :
      self._complements = set()

  @property
  def name(self):
    return self._name

  @property
  def length(self):
    return len(self.sequence)

  @property
  def sequence(self):
    return self._sequence

  @property
  def base_length(self):
    """Return the length of the nucleotide sequence."""
    return len(self.base_sequence)
    
  @property
  def base_sequence(self):
    """Breaks down the domain into non-composite domains."""
    if all(isinstance(s, str) for s in self._sequence):
      return self._sequence
    else :
      if all(isinstance(s, Domain) for s in self._sequence):
        return list(''.join(
          map(lambda x: ''.join(x.base_sequence), self._sequence)))
      else :
        # Implement mixed composite / non-composite case here
        raise NotImplementedError

  def _merge_constraints(self, con, con2):
    """Return a new list of unified constraints. """
    return map(self.iupac_union, zip(con,con2))

  def update_constraints(self, con):
    """Unify new and old constraint. """
    if not con or type(con) != list :
      raise ValueError("Must pass a constraints as a list.")

    # Implement this when needed!
    if not all(isinstance(s,str) for s in self.sequence + con) :
      raise NotImplementedError('Cannot update constraints on composite domains.')

    if len(self._sequence) != len(con):
      raise ValueError("Constraints differ in length.")

    new = self._merge_constraints(self.sequence, con)

    if '' in new :
      raise ValueError("Constraints cannot be satisfied")
    else :
      self._sequence = new

  def get_ComplementDomain(self, compseq):
    """This function initializes or updates a ComplementDomain. 

    Note: To return the complement, use the __invert__ operator.
    """
    # Implement when needed!!!
    if not all(isinstance(s,str) for s in self.sequence + compseq) :
      raise NotImplementedError('Cannot initialize composite ComplementDomain.')

    if self._ComplementDomain :
      self._ComplementDomain.update_constraints(compseq)
    else :
      if len(compseq) != len(self._sequence) :
        raise ValueError("Length of constraint != complementary constraint")
      complement = map(self.iupac_neighbor, self._sequence)
      new = self._merge_constraints(complement, compseq)
      if '' in new :
        raise ValueError("Constraints cannot be satisfied")
      self._ComplementDomain = ComplementDomain(self, new)
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
      raise NotImplementedError("Multiple complements not implemented")
      return self._complements
    else :
      return self._ComplementDomain

class ComplementDomain(Domain):
  """ Represents a ComplementDomain, assuming that every Domain has exactly one
  ComplementDomain. Note that this object is always defined in terms of an
  original domain.
  """
  def __init__(self, CompDomain, sequence=[]):
    """Create the ComplementDomain for a given domain. A ComplementDomain can
    only be initialized with respect to a given Domain. There can only exist
    one ComplementDomain for every Domain. 
    """
    self.id   = CompDomain.id 
    self._name = CompDomain._name + '*'
    self._ComplementDomain = CompDomain

    # Assign domain sequence 
    if not sequence or type(sequence) != list :
      raise ValueError("Must pass a sequence of constraints or domain objects.")
    else :
      assert all(isinstance(s, (str, Domain)) for s in sequence)
      self._sequence = sequence

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

class Complex(object):
  """A complex is a sequence & structure pair. 
  
  Sequence and structure can be specified on the domain or on the nucleotide
  level, but they have to be of the same length. Although not yet implemented,
  one may define special characters for secondary structures that are more
  diverse than a regual dot-parens string, e.g. 'h' = '(', '.', ')'

  Global Variables:
    id_counter: the next automatically assigned complex-ID in a system.
  """

  id_counter = 0

  def __init__(self, sequence=[], structure=[], name='', prefix='cplx'):
    """Initialization of the Complex object.

    Arguments:
      sequence <list>         = A list of sequence constraints on this domain.
      structure <list>        = A list of dot-parens characters corresponding
                                to the specified sequence.
      name <optional: str>    = Name of this domain. If not specified, an
                                automatic name is generated.
      prefix <optional: str>  = A prefix for autmated naming of Domains. 
                                Default: 'cplx' for 'complex'
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

    if sequence == [] :
      raise ValueError('Complex() requires Sequence and Structure Argument')

    if len(sequence) != len(structure) :
      raise ValueError("Complex() sequence and structure must have same length")
    self._sequence = sequence
    self._structure = structure

  @property
  def name(self):
    return self._name

  @name.setter
  def name(self, name):
    if name[-1]=='*' :
      raise ValueError('Invalid name for Complex Object!')
    if name[-1].isdigit() :
      raise ValueError('Complex name must not end with a digit!', name)
    self._name = name

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
      elif not all(isinstance(d, str) for d in seq) :
        raise NotImplementedError("mixed sequences are not supported")
      else :
        return seq
    return flatten(map(lambda x: my_base_sequence(x) + ['+'], lol))[:-1]

  @property
  def structure(self):
    return self._structure

  @property
  def lol_structure(self):
    indices = [-1] + [i for i, x in enumerate(self._structure) if x == "+"]
    indices.append(len(self._structure))
    return [self._structure[indices[i-1]+1: indices[i]] 
        for i in range(1, len(indices))]

  @property
  def nucleotide_structure(self):
    """Return the complex structure on nucleotide-level. """
    lol = zip(self.lol_sequence, self.lol_structure)
    def my_base_sequence(seq, sst) :
      if all(isinstance(d, Domain) for d in seq):
        tups = zip(seq, sst)
        return map(lambda (x,y): y * x.base_length, tups)
      elif not all(isinstance(d, str) for d in seq) :
        raise NotImplementedError("mixed sequences are not supported")
      else :
        return sst
    return flatten(map(lambda (x,y): my_base_sequence(x,y) + ['+'], lol))[:-1]

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

class TestTube(object):
  """A reaction network of nucleic acid complexes.

  **Description:**
    TestTube() objects are Nuskell's standard interface to enumerate and
    simulate nucleic acid systems.  Domain-level reaction networks are
    enumerated using the Python package ``Peppercorn'' and simulated using the
    Python package ``crnsimulator''. TestTube() can be exported to standard data
    structures, e.g. to verify the equivalence between two TestTube() objects.

    Single or multiple Complex() and/or Reaction() objects can be accessed,
    added and removed from the system.  TestTube() provides (optional) assert
    statements checking if Complex() and Domain() instances have been
    duplicated, but they might be time-consuming for large networks.

    Built-in functions can process domain-level networks to remove strands with
    wildcard-domains (also called history-domains). 
    There are a few options to replace a wildcard-strand:
      * replace with any single-stranded matching strand in the TestTube()
      * replace with all single-stranded matching strands in the TestTube()
      * replace with any matching strand in the TestTube()
      * replace with all matching strands in the TestTube()
      * remove the wildcard-domain.
    It is recommended to adjust experimental protocols to the method chosen here.

  **Structure:**
    TestTube() is based on networkx.MultiDiGraph(), with two types of nodes: (a)
    nuskell.Complex() and (b) nuskell.Reaction().  The graph is bipartite, edges
    are directed (from reactants to prodcuts), they connect reactants to a
    reaction node and a reaction node to products.
    
    TestTube() provides an additional *concentration* attribute and *constant*
    atribute to Complex() nodes, as well as a rate attribute to Reaction()
    nodes. These attributes are accessed when writing ODE systems and
    (optionally) updated after simulations. This means a TestTube() **can** be
    initialized without using Complex() and Reaction() objects, but by using a
    consistent naming scheme.

  **Developers:**
    It is recommended to store all other node attributes (e.g. free energies,
    etc.) directly in the nuskell.Complex() and nuskell.Reaction() objects.
    TestTube() does not provide an I/O interface for file formats. There is a
    separate TestTubeIO() object explicitly to parse and write compatible file
    formats (*.pil, *.dom, *.dna, etc.).

  **TODO:**
    The interface allows boolean operations for TestTube(): AND, OR, NOT, XOR
    The interface allows set operations for TestTube(): UNION, ...
  """

  # A global TestTube() variable to make (time-consuming) sanity-checks when
  # separate instances are subject to boolean or arithmetic opertions.
  sanitychecks = True
  warnings = True

  def __init__(self, complexes=None, reactions=None):
    """Initialization of TestTube object. 

    Arguments:
      complexes <optional: dic()> =  A dictionary of complex names that stores
        a tuple of [0]: the respective Complex() object or None, [1]: the
        concentration, and [2]: boolean value to specify if the concentration
        stays constant during a simulation: True = constant; False = initial
        For example: complexes['A'] = (None, 100, True) is a new
        species called 'A', which has no Complex() object initialized, and
        remains constant at 100 [nM] concentration.

        Note: The constant attribute defaults to False if not specified, in
        which case the concentrations specify the initial state and might be
        updated after a simulation.
    
    """
    self._RG = nx.MultiDiGraph() 

    if complexes :
      for name, data in complexes.items():
        if len(data) == 2 or len(data) == 3:
          assert isinstance(data[0], (Complex, None))
          assert isinstance(data[1], (float, None))
          if len(data) == 3:
            assert isinstance(data[2], bool)
        else :
          raise RuntimeError('wrong initialization of arguments for TestTube()')

        cplx = data[0]
        conc = data[1]
        const = data[2] if len(data) == 3 else False
        if cplx :
          self._RG.add_node(cplx, concentration=conc, constant=const) 
        else :
          self._RG.add_node(name, concentration=conc, constant=const) 

    if reactions :
      # Also reactions have to be uniquely adressable
      for name, react in reactions.items():
        if isinstance(react, Reaction) :
          self._RG.add_node(Reaction, rate=Reaction.rate)
          for r in Reaction.reactants :
            assert self._RG.has_node(r)
            self._RG.add_edge_from(r, Reaction)
          for p in Reaction.products:
            assert self._RG.has_node(p)
            self._RG.add_edge_from(Reaction, p)
        elif isinstance(react, list):
          assert len(react) == 3
          self._RG.add_node(name, rate=react[2])
          for r in react[0]:
            assert self._RG.has_node(r)
            self._RG.add_edge_from(r, name)
          for p in react[1]:
            assert self._RG.has_node(p)
            self._RG.add_edge_from(name, p)

    self._domains = None
    self._strands = None

  @property
  def ReactionGraph(self):
    return self._RG
  
  @ReactionGraph.setter
  def ReactionGraph(self, RG):
    self._RG = RG
    self._domains = None
    self._strands = None

  @property
  def complexes(self):
    # sorted = sorted(CG.nodes(data=True), key=lambda x: x[1]['concentration'], reverse=False)
    return self._RG.nodes()

  def complex_concentration(self, cplx):
    return self._RG.node[cplx]['concentration'], self._RG.node[cplx]['constant']

  def add_complex(self, cplx, (conc, const), sanitycheck=True):
    """Add a complex to the TestTube. 

    A new complex resets .domains and .strands
    """ 
    if not isinstance(cplx, Complex) :
      #TODO: do the formal stuff later
      raise NotImplementedError

    if self._RG.has_node(cplx):
      # Same Name, same object!
      assert self._RG[cplx] is cplx, Warning("Conflicting Complex Objects")
      assert self._RG[cplx]['concentration'] == conc, Warning("Conflicting complex concentrations")
      assert self._RG[cplx]['const'] == const, Warning("Conflicting complex concentrations")
    else :
      # NOTE: This might become inefficient at some point, but it has been
      # introduced to overcome issues with some translation schemes that
      # produce the same fuel strand multiple times.
      if sanitycheck and (cplx.sequence, cplx.structure) in map(
          lambda x: (x.sequence, x.structure), self._RG.nodes()):
        print 'WARNING: One complex, one name! Skipping complex:', cplx.name, \
            map(str, cplx.sequence), cplx.structure
      else :
        self._RG.add_node(cplx, concentration=conc, constant=const) 
        self._domains = None
        self._strands = None

  def rm_complex(self, cplx):
    """Remove a Complex from the TestTube. 

    Removing a complex resets .domains and .strands
    """
    if self._RG.has_node(cplx):
      if G.in_edges(cplx) or G.out_edges(cplx):
        raise RuntimeError("Cannot remove a complex engaged in reactions.")
      self._RG.remove_node(cplx)
      self._domains = None
      self._strands = None

  @property
  def strands(self):
    """Return a dictionary of strands present in the TestTube. 

    A strand is a nucleic-acid molecule connected by a single covalent
    backbone. Strands are named automatically, and their names may change
    whenever a new Complex is added to the TestTube.

    Returns: 
      strands[strand_1] = [Domain(X), Domain(Y), Domain(Z)]
    """
    if not self._strands :
      count = 0
      self._strands = dict()
      self._strand_names = dict()
      for cplx in self._RG.nodes_iter():
        for s in cplx.lol_sequence :
          strand = tuple(map(str,s))
          if strand not in self._strand_names :
            name = 'strand_{}'.format(count); count += 1
            self._strand_names[strand] = name
            self._strands[name] = s
    return self._strands

  @property
  def domains(self):
    """Return a dictionary of Domain Objects present in the TestTube. 
    
    Returns:
      domains[Domain.name] = Domain
    """
    if not self._domains :
      self._domains = dict()
      for cplx in self._RG.nodes_iter():
        for d in cplx.sequence :
          if d == '+' : continue
          if d.name in self._domains :
            assert self._domains[d.name] is d
          else :
            self._domains[d.name] = d
    return self._domains
 
  def enumerate_reactions(self, args):
    # call peppercorn to enumerate the reaction network starting from all
    # complexes in the testtube.
    # initialized a 'networkx' MultiDiGraph()
    raise NotImplementedError

  def add_irrev_reaction(self, reaction):
    # reaction = irreversible: [[react],[prod],rate]
    raise NotImplementedError

  @property
  def reactions(self):
    raise NotImplementedError

  def compose(self, other):
    assert isinstance(other, TestTube)
    combined = TestTube()
    return nx.compose(self.ReactionGraph, other.ReactionGraph)


  def __add__(self, other):
    assert isinstance(other, TestTube)
    combined = TestTube()

    # global TestTube() variable
    if not TestTube.sanitychecks :
      if TestTube.warnings :
        raise Warning('TestTube() - sanity checks turned off!')
      combined.ReactionGraph = nx.compose(self.ReactionGraph, other.ReactionGraph)

    elif False and len(other.complexes) > len(self.complexes) :
      combined.ReactionGraph = other.ReactionGraph
      map(lambda c : combined.add_complex(c, self.complex_concentration(c), 
        sanitycheck=True), self.complexes)

    else :
      combined.ReactionGraph = self.ReactionGraph
      map(lambda c : combined.add_complex(c, other.complex_concentration(c), 
        sanitycheck=True), other.complexes)
      #map(lambda (x, y): combined.add_complex(x, sanitycheck=True, **y), other.complexes)
    return combined

  def __radd__(self, other):
    # Reverse add is used for: sum([Testtube1, Testtube2, ...])
    if other == 0:
      return self
    else:
      return self.__add__(other)

class TestTubeIO(object):
  """A wrapper class to handle I/O of TestTube objects."""

  def __init__(self, ttube):
    assert isinstance(ttube, TestTube)
    self._testtube = ttube

  @property
  def testtube(self):
    return self._testtube

  def write_domfile(self, dom):
    domains = self._testtube.domains

    # Print Sequences
    for k, v in sorted(domains.items(), key=lambda x : x[1].id):
      if v.name[-1]=='*' : continue
      dom.write("sequence {:s} = {:s} : {:d}\n".format(
          v.name, ''.join(v.sequence), v.length))

    # Print Complexes
    for cplx in sorted(self._testtube.complexes):
      dom.write("{:s} = {:s} : {:s}\n".format(cplx.name, 
        ' '.join(map(str,cplx.sequence)), ' '.join(cplx.structure)))

  def load_domfile(self, domfile):
    raise NotImplementedError

  def write_pilfile(self, pil):
    """Write the contents of TestTube() into a pilfile. 

    This function supports a particular version of pilfiles, used by nuskell
    and peppercorn.

    sequence d1 = NNN : 3
    sequence d2 = NNNN : 4
    sequence hist = N : 1
    strand s1 = hist d1 d2 : 8
    strand s2 = d2* d1* : 7
    structure c1 = s1 : .(((((((+)))))))
    """
    domains = self._testtube.domains

    # Print Sequences
    for k, v in sorted(domains.items(), key=lambda x : x[1].id):
      if v.name[-1]=='*' : continue
      pil.write("sequence {:s} = {:s} : {:d}\n".format(
          v.name, ''.join(v.sequence), v.length))

    # Print Strands
    strands = self._testtube.strands
    for k, v in sorted(strands.items(), key=lambda x: int(x[0].split('_')[1])):

      pil.write("strand {:s} = {:s} : {:d}\n".format(k, 
        ' '.join(map(str, v)), sum(map(lambda x : x.length, v))))

    # Print Structures
    for cplx in sorted(self._testtube.complexes):
      pil.write("structure {:s} = ".format(cplx.name))
      for s in cplx.lol_sequence :
        strand = tuple(map(str,s))
        name = self._testtube._strand_names[strand]
        pil.write("{:s} ".format(name))
      pil.write(": {:s}\n".format(''.join(cplx.nucleotide_structure)))

  def load_pilfile(self, pilfile):
    """Parses a pil file written in KERNEL notation! """
    ppil = parse_pil_file(pilfile)

    def resolve_loops(loop):
      """ Return a sequence, structure pair from kernel format with parenthesis. """
      sequen = []
      struct = []
      for dom in loop :
        if isinstance(dom, str):
          sequen.append(dom)
          if dom == '+' :
            struct.append('+')
          else :
            struct.append('.')
        elif isinstance(dom, list):
          struct[-1] = '('
          old = sequen[-1]
          se, ss = resolve_loops(dom)
          sequen.extend(se)
          struct.extend(ss)
          sequen.append(old + '*' if old[-1] != '*' else old[:-1])
          struct.append(')')
      return sequen, struct

    domains = []
    for line in ppil :
      if line[0] == 'domain':
        domains.append(Domain(name=line[1]+'_', sequence = list('N'* int(line[2]))))
      elif line[0] == 'complex':
        name = line[1]+'_'
        sequence, structure = resolve_loops(line[2])

        for e in range(len(sequence)):
          d = sequence[e]
          if d == '+': continue
          if d[-1] == '*' : 
            dname = d[:-1] + '_'
            dom = filter(lambda x: x.name == dname, domains)
            if len(dom) < 1 :
              raise RuntimeError('Missing domain specification', d)
            elif len(dom) > 1 :
              raise RuntimeError('Conflicting matches for domain specification', d)
            sequence[e] = dom[0].get_ComplementDomain(list('R'*dom[0].length))

          else :
            dname = d + '_'
            dom = filter(lambda x: x.name == dname, domains)
            if len(dom) < 1 :
              raise RuntimeError('Missing domain specification', d)
            elif len(dom) > 1 :
              raise RuntimeError('Conflicting matches for domain specification', d)
            sequence[e] = dom[0]

        self._testtube.add_complex(Complex(sequence = sequence, structure =
          structure, name=name), (float("inf"), True))
      else :
        raise NotImplementedError('Weird expression returned from pil_parser!')

    return self._testtube

  def write_dnafile(self, fh, formal=[], crn=None, ts=None):
    """ Write a TestTube Object into VisualDSD *.dna format.

    Note: This function assumes that toehold domains are named starting with a
    't', history domains start with a 'h' and anti-sense domains end with '*'.

    Args:
      fh <filehandle>: The function prints to this filehandle.
      formal <optional: list(str)>: A list of formal species.
      crn <optional: list(list)>: a nuskell-style parsed CRN expression
      ts <optional: str>: name of the translation scheme

    """
    fh.write("(* File autogenerated by nuskell. ")

    if ts:
      fh.write("\n - Translation Scheme: {}".format(ts))
    if crn:
      fh.write("\n - Input CRN: \n")
      for rxn in crn :
        if len(rxn) > 2 and rxn[2] == 'reversible' :
          fh.write("    {} <=> {}\n".format(
            ' + '.join(rxn[0]), ' + '.join(rxn[1])))
        else :
          fh.write("    {} -> {}\n".format(
            ' + '.join(rxn[0]), ' + '.join(rxn[1])))

    fh.write("*)\n\n".format(crn))

    fh.write("def Fuel = 20\n")
    fh.write("def Formal = 5\n\n")

    first = True
    for cplx in sorted(self._testtube.complexes):

      if first:
        fh.write ('( ')
        first=False
      else :
        fh.write ('| ')

      if cplx.name in formal :
        fh.write("Formal * ")
      else :
        fh.write("constant Fuel * ")

      name = cplx.name
      sequ = cplx.sequence
      stru = cplx.structure

      ptab = pair_table(stru)

      dnaexpr = [[]]
      pos = 0
      for e, d in enumerate(ptab) :
        if d == '+' :
          flag = 'top' if flag == 'bound' else flag
          expr = 'cut'

        elif d == -1 :
          toe = '^' if sequ[e].name[0] == 't' else ''
          if sequ[e].name[-1] == '*' :
            flag = 'bottom'
            expr = sequ[e].name[:-1] + toe + '*'
          elif sequ[e].name[0] == 'h':
            flag = 'top'
            expr = '_'
          else :
            flag = 'top'
            expr = sequ[e].name + toe


        elif d > e : # '('
          flag = 'bound'
          toe = '^' if sequ[e].name[0] == 't' else ''
          if sequ[e].name[-1] == '*' :
            expr = sequ[e].name[:-1] + toe + '*'
          elif sequ[e].name[0] == 'h':
            raise Exception('unexpected bound history domain.')
          else :
            expr = sequ[e].name + toe

          dnaexpr.append([])
          pos += 1

        elif d < e : # ')'
          flag = 'bottom'
          expr = None
          pos -= 1
          if pos < 0 :
            raise Exception('too many closing base-pairs')
          continue
        else :
          raise Exception('strange case:', e, d)

        if dnaexpr[pos] == [] :
          dnaexpr[pos] = [[flag, expr]]
        else :
          dnaexpr[pos].append([flag, expr])

      # decode dnaexpr
      dnaflat = []
      for d in dnaexpr:
        for dd in d :
          dnaflat.append(dd)

      # PRINT TO FILE
      close = None
      for e, d in enumerate(dnaflat) :
        if d[1] == 'cut': 
          fh.write(close)
          close = None
          if e == len(dnaflat)-1:
            continue
          if d[0] == 'bottom' :
            fh.write('::')
          else :
            fh.write(':')
          continue
        
        if d[0] == 'bottom' :
          if close is None :
            fh.write('{')
            close = '}'
          elif close == ']' or close == '>' :
            fh.write('{}{'.format(close))
            close = '}'

        if d[0] == 'bound' :
          if close is None :
            fh.write('[')
            close = ']'
          elif close == '}' or close == '>' :
            fh.write('{}['.format(close))
            close = ']'

        if d[0] == 'top' :
          if close is None :
            fh.write('<')
            close = '>'
          elif close == '}' or close == ']' :
            fh.write('{}<'.format(close))
            close = '>'
        fh.write(" {} ".format(d[1]))
      if close :
        fh.write("{} (* {} *)\n".format(close, name))
      else:
        fh.write(" (* {} *)\n".format(name))

    fh.write (")\n")

  def load_dnafile(self, dnafile):
    raise NotImplementedError

class Reaction(object):
  """ A formal reaction object. 

  Stores common attributes such as reactants, products, rates, rtype, etc.
  """

  def __init__(self, reactants, products, rtype = None, rate = None):
    """ """
    self._reactants = sorted(reactants)
    self._products = sorted(products)
    self._type = rtype
    self._rate = rate

