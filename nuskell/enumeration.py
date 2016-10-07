import os
from nuskell.parser import parse_dom_file
import nuskell.include.peppercorn.enumerator as pepen
import nuskell.include.peppercorn.input as pepin
import nuskell.include.peppercorn.output as pepout
import nuskell.include.peppercorn.reactions as reactions

# Once peppercorn is its own package, we import it as depencency.
# try :
#   import peppercorn
# except ImportError:
#   sys.exit("""nuskell depends on the peppercorn package
#   -- download at http://dna.caltech.edu/peppercorn """)

def enumerate_crn_old(args, domfile):
  dom = parse_dom_file(domfile)
  efile, ctmp = enumerator_input(dom)
  enum = pepin.input_enum(efile)
  os.remove(efile)

  set_enumargs(enum, args)
  enum.enumerate()

  return get_enum_data(args, enum)

def enumerate_crn(args, pilfile):
  # A new interface that requires a pilfile rather than a domfile.
  enum = pepin.input_pil(pilfile)

  set_enumargs(enum, args)
  enum.enumerate()
  
  return get_enum_data(args, enum)

def get_enum_data(args, enum):
  # Write the output of peppercorn into the enumfile
  enumfile = args.output + '.enum'
  pepout.output_crn(enum, enumfile, output_condensed = True)

  # Post-process enumerator results to extract the condensed crn
  enum_crn = []
  with open(enumfile, 'r') as enu :
    for line in enu :
      react = line.split('->')
      react[0] = sorted([x.strip() for x in react[0].split('+')])
      react[1] = sorted([x.strip() for x in react[1].split('+')])
      enum_crn.append(react)

  # Extraction of complexes in the implementation CRN:
  init_cplxs = []
  for cx in enum.initial_complexes :
    cxs = []
    for sd in cx.strands:
      cxs.append('+')
      for ds in map(str, sd.domains):
        if ds == 'dummy' : 
          cxs.append('?')
        elif ds[-1] == '*':
          cxs.append([ds[:-1], '*'])
        else:
          cxs.append([ds])
    # remove the first '+' again
    if len(cxs)>0: cxs = cxs[1:]
    cx = [cx.name, cxs, list(cx.dot_paren_string())]
    init_cplxs.append(cx)

  # Extraction of complexes in the enumerated CRN:
  slow_cplxs = []
  for rs in enum.resting_states:
    name = str(rs)
    for cx in rs.complexes:
      cxs = []
      for sd in cx.strands:
        cxs.append('+')
        for ds in map(str, sd.domains):
          if ds[-1] == '*':
            cxs.append([ds[:-1], '*'])
          else:
            cxs.append([ds])
      # remove the first '+' again
      if len(cxs)>0: cxs = cxs[1:]
      cx = [name, cxs, list(cx.dot_paren_string())]
      slow_cplxs.append(cx)

  # Alternative way to parse the dom-file for implementation complexes
  #dom = parse_dom_file(domfile)
  #cplxs = dom[1] if len(dom)==2 else dom[0] # else: no sequence information
  return enum_crn, init_cplxs, slow_cplxs

def set_enumargs(enum, args):
  """Transfer options to Enumerator Object. 
  
  Do NOT set Nuskell-defaults here. Defaults are set with the argparse object
  of peppercorn: nuskell/include/peppercorn/enumerator.py 
  
  """

  #enum.MAX_REACTION_COUNT = 5000000
  #enum.MAX_COMPLEX_COUNT  = 100000
  #enum.MAX_COMPLEX_SIZE   = 1000
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

#DEPRICATED, used by enumerate_crn_old()
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

