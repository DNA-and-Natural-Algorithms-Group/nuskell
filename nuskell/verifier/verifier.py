
import os
import sys
from nuskell.parser import parse_crn_string, parse_dom_file
from nuskell.parser import split_reversible_reactions
import crn_bisimulation_equivalence
# But it doesn't import the pspace algorithm yet...
import crn_pathway_equivalence

def find(l, key):
  for i in range(len(l)):
    if l[i] == key:
      return i
  return None

def rotate(complex):
  def hardcopyList(l):
    if type(l) != list:
      return l
    return map(hardcopyList, l)
  complex = hardcopyList(complex)

  if "+" not in complex[0]:
    return complex        
  else:
    p = find(complex[0], "+")
    dom = complex[0][p + 1:] + ["+"] + complex[0][:p]

    # change parentheses appropriately
    dpr = complex[1]
    stack = []
    for i in range(p):
        if dpr[i] == "(": stack.append(i)
        elif dpr[i] == ")": stack.pop()
    for i in stack:
        dpr[i] = ")"
    stack = []
    for i in reversed(range(p + 1, len(dpr))):
        if dpr[i] == ")": stack.append(i)
        elif dpr[i] == "(": stack.pop()
    for i in stack:
        dpr[i] = "("
    
    dpr = dpr[p + 1:] + ["+"] + dpr[:p]
    return [dom, dpr]

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
  # TODO: Checks if two complexes are the same, assuming one has a history
  # domain. This function was here, as is. It needs to be checked to see if it
  # does what it should.
  if "+" in x[0]:
    if "+" not in y[0]:
      return False
    px = find(x[0], "+")
    py = find(y[0], "+")
    return patternMatch([x[0][:px], x[1][:px]],
                        [y[0][:py], y[1][:py]]) and \
           patternMatch([x[0][px + 1:], x[1][px + 1:]],
                        [y[0][py + 1:], y[1][py + 1:]])

  if len(x[0]) == 0:
    if len(y[0]) > 0:
      return False
    else:
      return True

  if x[0][0] != "?":
    if len(y[0]) == 0 or x[0][0] != y[0][0] or x[1][0] != y[1][0]:
      return False
    else:
      return patternMatch([x[0][1:], x[1][1:]],
                              [y[0][1:], y[1][1:]])
  else:
    for i in range(len(y) + 1):
      if patternMatch([x[0][1:], x[1][1:]],
                      [y[0][i:], y[1][i:]]):
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
    input_fs: A list of formal species
    init_cplxs: All complexes initially present in a DSD system (formal+constant)
    enum_cplxs: List of enumerated complexes.

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
  for x in init_cplxs :
    # x = ['A', ['?', ['t0'], ['d2']], ['.', '.', '.']]
    if '?' in x[1] :
      cnt = 0
      for y in enum_cplxs :
        if x[0] == y[0]: # formal species!
          enum_rename[y[0]] = x[0] + "_i"
          if x[0] in input_fs :
            enum_to_formal[x[0]+"_i"]=[x[0]]
          else :
            # this is just because I want to observe a case, 
            # should be save to remove this else statement.
            raise ValueError('Unexpected constant species')
        else :
          # TODO: patternMatch is untested, also there might be cases where one
          # has to use the 'rotate' function, which was used previously.
          if patternMatch(x[1:], y[1:]) :
            cnt += 1
            enum_rename[y[0]] = x[0] + "_" + str(cnt)
            if x[0] in input_fs :
              enum_to_formal[x[0]+"_"+str(cnt)] = [x[0]]
              remove_ihist.add(x[0]+"_i")
            else :
              # this is just because I want to observe a case, 
              # should be save to remove this else statement.
              raise ValueError('Unexpected constant species')
    else :
      if x[0] in input_fs :
        enum_to_formal[x[0]]=[x[0]]

  # removing initial history species, if replaceable
  map(lambda x: enum_to_formal.pop(x), remove_ihist)

  return enum_to_formal, enum_rename

def pre_process(enum_crn, input_fs, complexes, slow_cplxs):
  ##DEPRICATED: This is the original preprocessing routine.
  # Some things are implemented different than in the new version. Needs some
  # testing to see if the new version is actually doing the exact same thing.
  # There might be differences in patternMatching

  cs = map(lambda x: x[0], complexes)
  cs = filter(lambda x: x not in input_fs, cs)

  inter = {}
  dictionary = {}
  fsp = []
  rm = set()
  for x in complexes:
    if "?" not in x[1]:
      if x[0] in input_fs:
        inter[x[0]] = [x[0]]
        fsp.append(x[0])
      continue
    cnt = 0
    for y in slow_cplxs:
      if x[0] == y[0]:
        dictionary[y[0]] = x[0] + "_i"
        if x[0] in input_fs:
          inter[x[0]+"_i"] = [x[0]]
          fsp.append(x[0]+"_i")
        continue
      original = x[1:]
      target = y[1:]
      if len(original[0]) != len(target[0]): continue
      p = rotate(target)
      while True:
        flag = True
        for i in range(len(original[0])):
          if not (
              (original[0][i] == p[0][i] and original[1][i] == p[1][i]) or 
              (original[0][i] == "?" and p[1][i] == ".")):
            flag = False
        if flag:
          cnt += 1
          dictionary[y[0]] = x[0] + "_" + str(cnt)
          if x[0] in input_fs:
            inter[x[0]+"_"+str(cnt)] = [x[0]]
            fsp.append(x[0]+"_"+str(cnt))
            rm.add(x[0]+"_i")
          break
        if p == target:
          break
        p = rotate(p)

  for i in range(len(enum_crn)):
    [reactants, products] = enum_crn[i]
    def get_name(x):
      if x in dictionary.keys():
        x = dictionary[x]
      return x
    reactants = map(get_name, reactants)
    products = map(get_name, products)
    enum_crn[i] = [reactants, products]
  
  # preprocess fuel
  enum_crn = removeFuels(enum_crn, cs)
  enum_crn = sorted(map(lambda x: [sorted(x[0]), sorted(x[1])], enum_crn))
  enum_crn = remove_duplicates(enum_crn)

  # removing initial signals that are unnecessary
  fsp = set(fsp)
  for x in rm:
    if x in inter.keys(): del inter[x]
    if x in fsp: fsp.remove(x)
  norm = set(fsp)-rm

  flag = None
  while flag != norm:
    flag = set(list(norm))
    for [r,p] in enum_crn:
      if set(r).intersection(norm) == set(r):
        norm = norm.union(set(p))
  enum_crn = filter(lambda x: set(x[0]).intersection(norm) == set(x[0]), enum_crn)

  return inter, enum_crn, fsp

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
  cs = map(lambda x: x[0], init_cplxs)
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

