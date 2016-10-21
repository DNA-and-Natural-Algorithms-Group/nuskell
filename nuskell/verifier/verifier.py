
import os
import sys
import string # printRxn
from collections import Counter

from nuskell.parser import parse_crn_string, parse_dom_file
from nuskell.parser import split_reversible_reactions
import crn_bisimulation_equivalence
import crn_pathway_equivalence

def printRxn(rxn, inter={}):
  # First, replace every Instance of a variable by its interpretation
  # Second, print it
  react = map(lambda x: ' + '.join(
    sorted(inter[x].elements())) if x in inter else x, rxn[0])
  produ = map(lambda x: ' + '.join(
    sorted(inter[x].elements())) if x in inter else x, rxn[1])
  print ' + '.join(react), '->', ' + '.join(produ)

def patternMatch(x, y):
  """Matches two complexes if they are the same, ignoring history domains. 

  Note: The strand order of the second complex changes to the strand order of
  the first complex, if there is a rotation under which both complexes are
  patternMatched.

  Args: 
    x (Complex()) : A nuskell Complex() object.
    y (Complex()) : A nuskell Complex() object.

  Returns: True/False
  """
  if len(x.sequence) != len(y.sequence) :
    return False

  def pM_check(pMx, pMy):
    """Recursively parse the current sequences and structures. 

    Args: 
      pMx [seqX,strX]: A list of two lists (sequence, structrure)
      pMy [seqY,strY]: A list of two lists (sequence, structrure)

    Returns: True/False
    """
    if len(pMx[0]) == 0 :
      return True

    if (pMx[0][0] != '?' and pMy[0][0] != '?') and \
        (pMx[0][0] != pMy[0][0] or pMx[1][0] != pMy[1][0]):
          return False
    return pM_check([pMx[0][1:], pMx[1][1:]], [pMy[0][1:], pMy[1][1:]])

  pMx = [x.sequence, x.structure]
  pMy = [y.sequence, y.structure]
  if pM_check(pMx,pMy) :
    return True
  elif '+' in x.sequence and '+' in y.sequence :
    for yr in y.rotate :
      pMy = [yr.sequence, yr.structure]
      if pM_check(pMx,pMy) :
        return True
  return False

def removeFuels(crn, fuel):
    crn = [[filter(lambda s: s not in fuel, rxn[0]),
            filter(lambda s: s not in fuel, rxn[1])]
           for rxn in crn]
    return crn

def remove_duplicates(l):
    r = []
    if len(l) == 0: return []
    l.sort()
    while len(l) > 1:
        if l[0] != l[1]:
            r.append(l[0])
        l = l[1:]
    r.append(l[0])
    return r

def get_interpretation(input_fs, init_cplxs, enum_cplxs):
  """Infer a partial interpretation from species names.

  Some schemes, e.g. Soloveichik et al. (2010), use history domains.  During
  enumeration, these are treated as regular unique domains, but that leads to a
  larger network than necessary. This routine identifies species with wildcards
  to remove them from the enumerated network. This can make the CRN
  exponentially easier to verify.

  Args:
    input_fs:   A list of formal species
    init_cplxs: A TestTube object with all complexes initially present in a 
                DSD system (formal+constant)
    enum_cplxs: A TestTube object that contains enumerated complexes.

  Returns:
    enum_to_formal (dict): A mapping between enumerated species and their formal
                           interpretation, e.g. enum_to_formal['A_1']=['A']
    enum_rename (dict):    A list of names in the enumerated network that are to
                           be replaced with names of formal species 
                           e.g. enum_rename['i86'] = ['A_1']
  """
  enum_to_formal = dict() # Interpretation: dict['A_i'] = Counter('A':1)
  enum_rename = dict()    # Rename 
  remove_ihist = set()
  for nx, x in init_cplxs.complexes.items() :
    if '?' in x.sequence :
      cnt = 0
      for ny, y in enum_cplxs.complexes.items() :
        if nx == ny : # formal species!
          enum_rename[ny] = nx + "_i"
          if nx in input_fs :
            enum_to_formal[nx+"_i"]=Counter(nx)
          else :
            # SB: this is just because I want to observe a case, 
            # should be save to remove this else statement.
            raise ValueError('Unexpected constant species')
        else :
          if patternMatch(x, y) :
            cnt += 1
            enum_rename[ny] = nx + "_" + str(cnt)
            if nx in input_fs :
              enum_to_formal[nx+"_"+str(cnt)] = Counter(nx)
              remove_ihist.add(nx+"_i")
            else :
              # SB: this is just because I want to observe a case, 
              # should be save to remove this else statement.
              raise ValueError('Unexpected constant species')
    else :
      if nx in input_fs :
        enum_to_formal[nx]=Counter(nx)

  # removing initial history species, if replaceable
  map(lambda x: enum_to_formal.pop(x), remove_ihist)

  return enum_to_formal, enum_rename

def verify(input_crn, enum_crn, init_cplxs, enum_cplxs, 
    method = 'bisimulation', verbose = False):
  """Wrapper-function to speed up the verification of a CRN translation.

  This function serves two purposes: First, it finds a mapping between species
  in the formal CRN and species in the enumerated CRN based on species names.
  Second, it reduces the size of the enumerated CRN. First of all, fuel and
  waste species are ignored for verification and duplicate species are removed.
  Some translation schemes use history domains, which are only specified for
  products, but not for reactants. This function removes all species with
  unspecified history domains, given the same species with a specified history
  domain exists in the system. For schemes with history domains, the
  interpretation may map multiple species of the enumerated CRN to one formal
  species.

  Args: 
    input_crn (string): formal input CRN
    enum_crn (list of lists): peppercorn-enumeration of the initial complexes 
    init_cplxs (TestTube()): A nuskell TestTube object with all initial species.
    enum_cplxs (TestTube()): A nuskell TestTube object with the enumerated CRN
    method (string): equivalence notion
    verbose (bool): print status information during verification

  Returns:
    True: formal and implementation CRN are equivalent
    False: formal and implementation CRN are not equivalent
  """
  interactive = False

  (input_crn, input_fs, input_cs) = parse_crn_string(input_crn) 
  irrev_crn = split_reversible_reactions(input_crn)

  # TODO: The format of interpret may change for a new verification interface
  # NOTE: enum_rename is empty for schemes without history domains
  interpret, enum_rename = get_interpretation(input_fs, init_cplxs, enum_cplxs)

  if verbose:
    print "======================="
    print "== CRN preprocessing =="
    print "======================="
    print "Original enumerated CRN:"
    for r in enum_crn: printRxn(r)

  if enum_rename :
    for i in range(len(enum_crn)):
      [reactants, products] = enum_crn[i]
      def get_name(x):
        if x in enum_rename.keys():
          x = enum_rename[x]
        return x
      reactants = map(get_name, reactants)
      products = map(get_name, products)
      enum_crn[i] = [reactants, products]

  # Reduce the enumerated CRN, 
  cs = map(str, init_cplxs.complexes)
  cs = filter(lambda x: x not in input_fs, cs)

  if verbose:
    print "Renamed enumerated CRN:"
    for r in enum_crn: printRxn(r)
 
  enum_crn = removeFuels(enum_crn, cs)
  enum_crn = sorted(map(lambda x: [sorted(x[0]), sorted(x[1])], enum_crn))
  enum_crn = remove_duplicates(enum_crn) # TODO: WHAT FOR?

  if verbose:
    print "Processed enumerated CRN:"
    for r in enum_crn: printRxn(r)
 
  # Get rid of all reactions that start with an *initial* history domain but
  # then later can get replaced by a produce molecule.
  if enum_rename:
    [prev, total] = [set(), set(interpret.keys())]
    while prev != total:
      prev = set(list(total))
      for [r,p] in enum_crn:
        if set(r).intersection(total) == set(r):
          total = total.union(set(p))
    enum_crn = filter(
        lambda x: set(x[0]).intersection(total) == set(x[0]), enum_crn)

  if verbose:
    print "History-reduced enumerated CRN:"
    for r in enum_crn: printRxn(r)
    print "Interpreted enumerated CRN:"
    for r in enum_crn: printRxn(r, interpret)
    print "======================="

  if method == 'bisimulation':
    return crn_bisimulation_equivalence.test(
        (irrev_crn, input_fs), (enum_crn, input_fs), verbose)
  elif method == 'bisim-loopsearch':
    return crn_bisimulation_equivalence.test(
      (irrev_crn, input_fs), (enum_crn, input_fs), verbose, [[],[]], 'pspace')
  elif method == 'bisim-wholegraph':
    return crn_bisimulation_equivalence.test(
      (irrev_crn, input_fs), (enum_crn, input_fs), verbose, [[],[]], 'whole')

  elif method == 'pathway':
    # NOTE: Adaptation to pathway interface
    pinter = dict()
    for k,v in interpret.items() :
      v = sorted(v.elements())[0]
      pinter[k]=[v]
    return crn_pathway_equivalence.test((irrev_crn, input_fs), 
        (enum_crn, pinter.keys()), pinter, verbose, False, interactive)
  elif method == 'integrated':
    # TODO: integrated-hybrid -> first consider some species as formal for
    # pathway decomposition, then do bisimulation. This is necessary for
    # history domains in some schemes, but it can be used for more general
    # things. E.g. if you make reversible reactions formal, then it will get
    # accepted and bisimulation can do the rest. In any case, the current
    # implementation does not cover the full hybrid theory, only some special
    # cases and some kind of bisimulation.
    pinter = dict()
    for k,v in interpret.items() :
      v = sorted(v.elements())[0]
      pinter[k]=[v]
    return crn_pathway_equivalence.test((irrev_crn, input_fs), 
        (enum_crn, pinterpret), pinter, verbose, True, interactive)
  else:
    print "Verification method unknown."
    return False

