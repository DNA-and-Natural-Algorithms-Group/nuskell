
import os
import sys
from nuskell.parser import parse_crn_string, parse_dom_file
from nuskell.parser import split_reversible_reactions
import crn_bisimulation_equivalence
import crn_pathway_equivalence

try :
  # TODO: only used in do_enumerator_things()
  import peppercorn
except ImportError:
  sys.exit("""nuskell depends on the peppercorn package
  -- download at http://dna.caltech.edu/peppercorn """)

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

def patternMatch(x, y):
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

def pre_process(enum_crn, input_fs, complexes, slow_cplxs):
  """ this might modify enum_crn"""

  # Extract the constant species from all complexes
  cs = map(lambda x: x[0], complexes)
  cs = filter(lambda x: x not in input_fs, cs)

  inter = {}
  dictionary = {}
  fsp = []
  rm = set()
  #TODO why?
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

#TODO DEPRICATED
def enumerator_input(dom):
  import random, string
  efile = "".join(random.sample(string.letters + string.digits, 8)) + "._tmp"
  F = open(efile, "w")
  if len(dom) == 2:   # there is sequence information
    complexes = dom[1]
    for i in dom[0]:
      print >> F, "domain", i[0], ":", i[1]
  else:
    complexes = dom[0]
  print >> F, "domain dummy : 15"

  strand_n = 0
  out_strands = {}
  out_complexes = []

  for i in complexes:
    c = [i[0],"",""]
    s = ""
    i[1].append("+")
    for j in i[1]:
      if j == "?":
        s += " dummy"
      if type(j) == list:
        if len(j) == 2:
          j = j[0] + j[1]
        else:
          j = j[0]
        s += " " + j
      elif j == "+":
        if s not in out_strands.keys():
          strand_n += 1
          sname = "strand"+str(strand_n)
          out_strands[s] = sname
        else:
          sname = out_strands[s]
        s = ""
        c[1] += " "+sname
    i[1] = i[1][:-1]
    for j in i[2]:
      if j == "+": j = " + "
      c[2] += j
    out_complexes.append(c)

  for s in out_strands.keys():
    sname = out_strands[s]
    print >> F, "strand", sname, ":", s

  for [c0,c1,c2] in out_complexes:
    print >> F, "complex", c0, ":"
    print >> F, c1
    print >> F, c2
    print >> F
  F.close()
  return efile, complexes

def set_enumargs(enum, args):
  """Transfer options to Enumerator Object. 
  
  Set Nuskell-defaults here.
  """

  enum.MAX_COMPLEX_COUNT  = 10000
  enum.MAX_REACTION_COUNT = 50000
  enum.MAX_COMPLEX_SIZE   = 100
  #enum.k_fast = 2.0
  #enum.REJECT_REMOTE = True

  if args.k_slow is not None:
    enum.k_slow = args.k_slow
  if args.k_fast is not None:
    enum.k_fast = args.k_fast

  if args.MAX_REACTION_COUNT is not None:
    enum.MAX_REACTION_COUNT = args.MAX_REACTION_COUNT

  if args.MAX_COMPLEX_COUNT is not None:
    enum.MAX_COMPLEX_COUNT = args.MAX_COMPLEX_COUNT

  if args.MAX_COMPLEX_SIZE is not None:
    enum.MAX_COMPLEX_SIZE = args.MAX_COMPLEX_SIZE

  if args.RELEASE_CUTOFF is not None:
    enum.RELEASE_CUTOFF = args.RELEASE_CUTOFF

  if args.RELEASE_CUTOFF_1_1 is not None:
    enum.RELEASE_CUTOFF_1_1 = args.RELEASE_CUTOFF_1_1

  if args.RELEASE_CUTOFF_1_N is not None:
    enum.RELEASE_CUTOFF_1_N = args.RELEASE_CUTOFF_1_N

  if args.REJECT_REMOTE is not None:
    enum.REJECT_REMOTE = args.REJECT_REMOTE

  if args.UNZIP is not None:
    enum.UNZIP = args.UNZIP

  if args.LEGACY_UNZIP is not None:
    enum.LEGACY_UNZIP = args.LEGACY_UNZIP

  # TODO: what is this?
  # enum.DFS = not args.bfs

  # Modify enumeration events based on command line options.
  if args.ignore_branch_3way:
    if reactions.branch_3way in enum.FAST_REACTIONS:
      enum.FAST_REACTIONS.remove(reactions.branch_3way)

  if args.ignore_branch_4way:
    if reactions.branch_4way in enum.FAST_REACTIONS:
      enum.FAST_REACTIONS.remove(reactions.branch_4way)

  return

def verify(input_crn, enum_crn, complexes, slow_cplxs, 
    method = 'bisimulation', verbose = True):
  """Prepare the verification of a translation scheme.

  This function initializes the:

    *) formal CRN
    *) implementation CRN (formal and constant species)
    *) enumerated implementation CRN (formal, constant and enumerated species)
    *) interpretation CRN (a mapping from implementation to the formal CRN)

  ... and then calls the verification method chosen by the user in order to 
  compute whether an implementation CRN is equivalent to the formal CRN:

    *) pathway (Seung Woo Shin's notion of pathway equivalence)
    *) pathway-integrated (Seung Woo Shin's integrated hybrid notion of pathway
    equivalence)
    *) bisimulation (Robert Johnson's ...)
    *) bisimulation-xyz (Quin Dong's ...)
    ...

  Args: 
    input_crn: formal input CRN
    enum_crn: peppercorn-enumerated implementation CRN
    complexes: all complexes in the (non-enumerated) implementation CRN
    slow_cplxs: resting states in the enumerated implementation CRN
    method: chose equivalence notion
    verbose: print more information to STDOUT

  Returns:
    True: formal and implementation CRN are equivalent
    False: formal and implementation CRN are not equivalent
  
  """
  interactive = False

  # Parse the CRN
  (input_crn, input_fs, input_cs) = parse_crn_string(input_crn) 
  irrev_crn = split_reversible_reactions(input_crn)

  print 'f', input_fs
  print 'e', enum_crn

  # TODO: figure out what exactly pre_process does and then 
  # go on...

  # Extract the constant species from all complexes
  cs = map(lambda x: x[0], complexes)
  cs = filter(lambda x: x not in input_fs, cs)

  # Initialization of an interpretation dictionary, a formal species set, and
  # the enumerated species as renamed in the preprocessing code

  inter = {}  # { A_i : ['A'] }
  fsp = set() # (['A_i', 'B_1', 'X_1'])

  dictionary = {}
  rm = set()
  for x in complexes:
    # x example: ['A', ['?', ['d1'], ['d2'], ['d3']], ['.', '.', '.', '.']]
    if "?" not in x[1]:
      if x[0] in input_fs:
        inter[x[0]] = [x[0]]
        fsp.add(x[0])
      continue
    cnt = 0
    for y in slow_cplxs:
      if x[0] == y[0]:
        dictionary[y[0]] = x[0] + "_i"
        if x[0] in input_fs:
          inter[x[0]+"_i"] = [x[0]]
          fsp.add(x[0]+"_i")
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
            fsp.add(x[0]+"_"+str(cnt))
            rm.add(x[0]+"_i")
          break
        if p == target:
          break
        p = rotate(p)

  inter, enum_crn, fsp = pre_process(
      enum_crn, 
      input_fs, 
      complexes, # need these! (cs = compelexes - input_fs)
      slow_cplxs)

  print 'i', inter
  print 'e', enum_crn
  print 'f', fsp

  # fs = formal species; cs = fuels or constant species
  if method == 'bisimulation':
    return crn_bisimulation_equivalence.test(
        (irrev_crn, input_fs), (enum_crn, input_fs), verbose)
  elif method == 'pathway':
    return crn_pathway_equivalence.test(
        (irrev_crn, input_fs), 
        (enum_crn, fsp), inter, verbose, False, interactive)
  elif method == 'integrated':
    return crn_pathway_equivalence.test(
        (irrev_crn, input_fs), 
        (enum_crn, fsp), inter, verbose, True, interactive)
  else:
    print "Verification method unknown."
    return False

