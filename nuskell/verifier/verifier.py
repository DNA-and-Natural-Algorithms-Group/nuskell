
import os
import sys
import string # printRxn

from nuskell.parser import parse_crn_string, parse_dom_file
from nuskell.parser import split_reversible_reactions
import crn_bisimulation_equivalence
# But it doesn't import the pspace algorithm yet...
import crn_pathway_equivalence

def printRxn(rxn, inter = {}):
  #NOTE: temporarily taken from crn_pathway_equivalence
  first = True
  for x in rxn[0]:
    if x[0] not in string.letters:
      if x in inter.keys() and inter[x]==[]:
        x = "w" + x
      else:
        x = "i" + x
    if not first:
      print "+",
    else:
      first = False
    print x,
  print "->",
  first = True
  for x in rxn[1]:
    if x[0] not in string.letters:
      if x in inter.keys() and inter[x]==[]:
        x = "w" + x
      else:
        x = "i" + x
    if not first:
      print "+",
    else:
      first = False
    print x,
  print

def pM2(original,target):
  if len(original[0]) != len(target[0]): 
    return False
  p = rotate(target)
  while True:
    flag = True
    for i in range(len(original[0])):
      if not (
          (original[0][i] == p[0][i] and original[1][i] == p[1][i]) or 
          (original[0][i] == "?" and p[1][i] == ".")):
        flag = False
    if flag:
      return True
    if p == target:
      return False
    p = rotate(p)
  raise SystemExit('wtf')

def patternMatch(x, y):
  # TODO: Checks if two complexes are the same, ignoring history domains
  if len(x[0]) != len(y[0]) :
    return False

  if len(x[0]) == 0 :
    return True

  if (x[0][0] != '?' and y[0][0] != '?') and \
      (x[0][0] != y[0][0] or x[1][0] != y[1][0]):
        return False

  return patternMatch([x[0][1:], x[1][1:]], [y[0][1:], y[1][1:]])

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
  enum_to_formal = {}
  enum_rename = {}
  remove_ihist = set()
  for nx, x in init_cplxs.complexes.items() :
    # x = ['A', ['?', ['t0'], ['d2']], ['.', '.', '.']]
    if '?' in x.sequence :
      cnt = 0
      for ny, y in enum_cplxs.complexes.items() :
        if nx == ny : # formal species!
          enum_rename[ny] = nx + "_i"
          if nx in input_fs :
            enum_to_formal[nx+"_i"]=[nx]
          else :
            # this is just because I want to observe a case, 
            # should be save to remove this else statement.
            raise ValueError('Unexpected constant species')
        else :
          # TODO: patternMatch is untested, also there might be cases where one
          # has to use the 'rotate' function, which was used previously.
          if '+' in x.sequence :
            # In this case we need to rotate one complex to see if there is a
            # pattern match
            raise NotImplementedError
          pMx = [x.sequence, x.structure]
          pMy = [y.sequence, y.structure]
          if patternMatch(pMx, pMy) :
            cnt += 1
            enum_rename[ny] = nx + "_" + str(cnt)
            if nx in input_fs :
              enum_to_formal[nx+"_"+str(cnt)] = [nx]
              remove_ihist.add(nx+"_i")
            else :
              # this is just because I want to observe a case, 
              # should be save to remove this else statement.
              raise ValueError('Unexpected constant species')
    else :
      if nx in input_fs :
        enum_to_formal[nx]=[nx]

  # removing initial history species, if replaceable
  map(lambda x: enum_to_formal.pop(x), remove_ihist)

  return enum_to_formal, enum_rename

def verify(input_crn, enum_crn, init_cplxs, enum_cplxs, 
    method = 'bisimulation', verbose = True):
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
    input_crn: formal input CRN
    enum_crn: peppercorn-enumeration of the initial complexes 
    init_cplxs: formal and complex species that are initially present
    enum_cplxs: resting states of the enumerated CRN
    method: equivalence notion
    verbose: print status information during verification

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

  enum_crn = removeFuels(enum_crn, cs)
  enum_crn = sorted(map(lambda x: [sorted(x[0]), sorted(x[1])], enum_crn))
  enum_crn = remove_duplicates(enum_crn) # TODO: WHAT FOR?

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

  if False : #Temporary hack to get implementation CRNs
    print "Enumerated: CRN"
    print enum_crn
    for r in enum_crn:
      printRxn(r)
    raise SystemExit
  
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
    return crn_pathway_equivalence.test(
        (irrev_crn, input_fs), 
        (enum_crn, interpret.keys()), interpret, verbose, False, interactive)
  elif method == 'integrated':
    # TODO: integrated-hybrid -> first consider some species as formal for
    # pathway decomposition, then do bisimulation. This is necessary for
    # history domains in some schemes, but it can be used for more general
    # things. E.g. if you make reversible reactions formal, then it will get
    # accepted and bisimulation can do the rest. In any case, the current
    # implementation does not cover the full hybrid theory, only some special
    # cases and some kind of bisimulation.
    return crn_pathway_equivalence.test(
        (irrev_crn, input_fs), 
        (enum_crn, interpret.keys()), interpret, verbose, True, interactive)
  else:
    print "Verification method unknown."
    return False

