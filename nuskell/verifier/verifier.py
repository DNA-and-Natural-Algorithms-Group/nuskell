
import os
import sys
import nuskell.parser.dom_parser as dom_parser
import nuskell.parser.crn_parser as crn_parser
import crn_bisimulation_equivalence
import crn_pathway_equivalence

try :
  # TODO: only used in do_enumerator_things()
  import enumerator
except ImportError:
  sys.exit("""nuskell depends on the enumertor package
  -- download at http://dna.caltech.edu/enumerator """)

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

def post_process(enum_crn, inp_fs, cs, complexes, slow_cplxs):
  """ this might modify enum_crn"""
  inter = {}
  dictionary = {}
  fsp = []
  rm = set()
  for x in complexes:
    if "?" not in x[1]:
      if x[0] in inp_fs:
        inter[x[0]] = [x[0]]
        fsp.append(x[0])
      continue
    cnt = 0
    for y in slow_cplxs:
      if x[0] == y[0]:
        dictionary[y[0]] = x[0] + "_i"
        if x[0] in inp_fs:
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
          if not ((original[0][i] == p[0][i] and original[1][i] == p[1][i]) or (original[0][i] == "?" and p[1][i] == ".")):
            flag = False
        if flag:
          cnt += 1
          dictionary[y[0]] = x[0] + "_" + str(cnt)
          if x[0] in inp_fs:
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

def do_enumerator_things(efile):
  import enumerator
  import enumerator.input as enum_in
  import enumerator.output as enum_out
  condense = True

  enumerator = enum_in.input_enum(efile)
  enumerator.MAX_COMPLEX_COUNT = 10000
  enumerator.MAX_REACTION_COUNT = 50000
  enumerator.MAX_COMPLEX_SIZE = 100
  #enumerator.k_fast = 2.0
  #enumerator.REJECT_REMOTE = True
  #TODO: check this and uncomment
  #enumarg(enumerator)
  enumerator.enumerate()
  enum_out.output_crn(enumerator,efile, output_condensed = condense)

  F = open(efile, "r")
  e_crn = []
  for line in F:
    line = line.split("->")
    line[0] = line[0].split("+")
    line[1] = line[1].split("+")
    line[0] = map(lambda x: x.strip(), line[0])
    line[1] = map(lambda x: x.strip(), line[1])
    line[0] = sorted(line[0])
    line[1] = sorted(line[1])
    e_crn.append(line)
  F.close()

  # convert the output from state enumerator to cmpdna format
  slow_cplxs = []
  for z in enumerator.resting_states:
    name = str(z)
    for y in z.complexes:
      y1 = []
      for x in y.strands:
        y1.append("+")
        for w in x.domains:
          w = str(w)
          if w[-1] == "*":
            y1.append([w[:-1], "*"])
          else:
            y1.append([w])
      if len(y1)>0: y1 = y1[1:]
      y = [name, y1, list(y.dot_paren_string())]
      slow_cplxs.append(y)

  return e_crn, slow_cplxs

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

def verify(input_crn, domfile, method = 'bisimulation', verbose = True):
  """ Initilize the verification of a translation scheme
  """

  interactive = False

  # Parse the CRN
  (inp_crn, inp_fs, inp_cs) = crn_parser.parse_string(input_crn) 
  t = inp_crn
  inp_crn = []
  for [x, r, p] in t:
    inp_crn.append([r,p])
    if x == "reversible":
      inp_crn.append([p,r])

  # Parse the DOM
  dom = dom_parser.parse(domfile)

  # Generate an input file for state enumerator
  efile, complexes = enumerator_input(dom)

  #print efile, complexes

  # constant species
  cs = map(lambda x: x[0], complexes)
  cs = filter(lambda x: x not in inp_fs, cs)
  #print "constant species:", cs

  # call state enumerator
  enum_crn, slow_cplxs = do_enumerator_things(efile)
  os.remove(efile)

  # postprocess
  inter, enum_crn, fsp = post_process(enum_crn, inp_fs, cs, complexes, slow_cplxs)
  
  # fs = formal species; cs = fuels or constant species
  if method == 'bisimulation':
      return crn_bisimulation_equivalence.test((inp_crn, inp_fs), (enum_crn, inp_fs), verbose)#, inter)
  elif method == 'pathway':
      return crn_pathway_equivalence.test((inp_crn, inp_fs), (enum_crn, fsp), inter, verbose, False, interactive)
  elif method == 'integrated':
      return crn_pathway_equivalence.test((inp_crn, inp_fs), (enum_crn, fsp), inter, verbose, True, interactive)


