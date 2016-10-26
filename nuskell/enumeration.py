import os
from nuskell.parser import parse_dom_file
import nuskell.include.peppercorn.enumerator as pepen
import nuskell.include.peppercorn.input as pepin
import nuskell.include.peppercorn.output as pepout
import nuskell.include.peppercorn.reactions as reactions
from nuskell.include.peppercorn.condense import condense_resting_states

from nuskell.objects import Domain, Complex, TestTube


# Once peppercorn is its own package, we import it as depencency.
# try :
#   import peppercorn
# except ImportError:
#   sys.exit("""nuskell depends on the peppercorn package
#   -- download at http://dna.caltech.edu/peppercorn """)

def peppercorn_enumerate(args, pilfile, solution, verbose = False):
  # Prepare
  enum = pepin.input_pil(pilfile)
  set_enumargs(enum, args)

  # Do it 
  enum.enumerate()

  # Get Data
  enum_crn = get_enum_crn(enum, solution, verbose)
  enum_solution = get_enum_solution(enum, solution, verbose)

  return enum_crn, enum_solution

def enum_cplx_rename(x) :
  # Translating enum output to be compatible with nuskell Objects.
  x = 'e'+x if x[0].isdigit() else x
  return x +'p' if x[-1].isdigit() else x

# NOTE: This function is obsolete, because domains are initialized using
# the original solution object
#def enum_domain_rename(x) :
#  # Translating enum output to be compatible with nuskell Objects.
#  # ... and do some weird translation of history domain names
#  x = 'e'+x if x[0].isdigit() else x
#
#  if x == 'dummy':
#    return 'dummy'
#  elif x[-1]=='*':
#    return x[:-1] + 'p' + '*'
#  elif x[-1].isdigit() :
#    return x + 'p'
#  else :
#    return x

def get_enum_crn(enum, solution, verbose):
  # Extract condensed reactions
  condense_options={}
  condensed = condense_resting_states(enum, **condense_options)
  reactions = condensed['reactions']

  enum_crn = []
  for r in sorted(reactions): #utils.natural_sort(reactions):
    react = map(enum_cplx_rename, map(str, r.reactants))
    prod  = map(enum_cplx_rename, map(str, r.products))
    if verbose :
      print ' + '.join(react), '->', ' + '.join(prod)
    enum_crn.append([react, prod])

  return enum_crn


def get_enum_solution(enum, solution, verbose):
  # Extraction of complexes in the enumerated CRN:
  enum_cplxs = TestTube()
  enum_cplxs += solution
  for rs in enum.resting_states:
    enum_cplxs.load_enum_cplxs(rs.complexes,
        enum_cplx_rename)

  return enum_cplxs

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

