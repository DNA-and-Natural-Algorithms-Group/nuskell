#/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from nuskell.compiler import translate
from nuskell.enumeration import TestTubePeppercornIO
from nuskell.verifier import preprocess, verify
from nuskell.objects import TestTube, TestTubeIO
from nuskell.parser import parse_crn_string
from nuskell.parser import split_reversible_reactions
from nuskell.include.peppercorn.enumerator import get_peppercorn_args

def main():
  """Compare different Tranlation schemes for different CRNs.

  A number of descriptors are gathered in the main loop, stored in a DataFrame
  and then plotted.
  """
  parser = argparse.ArgumentParser()
  parser.add_argument("--ts_dir", action = 'store', default = "examples/ts/",
      help="Specify directory that contains translation schemes" + \
        "files that has a *.ts ending will be interpreted as translation scheme")

  parser.add_argument("--verify", default='bisimulation', action = 'store',
      help="Specify name a verification method: \
          (bisimulation, pathway, integrated, bisim-loop-search,\
          bisim-depth-first, bisim-whole-graph)") 

  parser = get_peppercorn_args(parser)

  args = parser.parse_args()

  ts_list = args.ts_dir

  # read CRN from STDIN
  # input_crn = sys.stdin.readlines()
  # input_crn = "".join(input_crn)

  # A few common tests to find bugs in reaction schemes
  reactions = [
      'A -> ', 
      '-> B', 
      'A -> B', 
      'A <=> B', 
      'A -> A + A', 
      'A + A -> A', 
      'A + B -> B + B' ]

  networks = [
      'A + B -> B + B; B + C -> C + C; A + C -> A + A',
      'A <=> A + A; A + B -> B + B; B -> ; A + C -> ; C <=> C + C',
      'A <=> X; B + X <=> Y + A; C + Y + X <=> Z + B + A']


  plotdata = [] # Scheme, CRN, Cost, Speed
  for crn in reactions + networks:
    for scheme in os.listdir(args.ts_dir) :
      if scheme[-3:] != '.ts' : 
        #print "# Ignoring file:", args.ts_dir + scheme
        continue

      #if scheme[0] != 's' and scheme[0] != 'q' : continue

      (fcrn, fs, cs) = parse_crn_string(crn) 
      fcrn = split_reversible_reactions(fcrn)

      # Settings that might be subject to change later
      args.REJECT_REMOTE = False
      args.MAX_COMPLEX_COUNT = 10000

      print '\n# Scheme:', scheme, 'CRN:', crn
      solution, _ = translate(crn, args.ts_dir + scheme, args.verbose)
      
      # TRANSLATION
      # NOTE: The actual number of nucleotides is smaller for translations using history domains.
      print 'Number of Nucleotides:',
      print sum(map(len, 
        map(lambda x: x.nucleotide_sequence, 
          solution.complexes.values())))
          #filter(lambda x: x.name not in fs, solution.complexes.values()))))

      cost = sum(map(len, map(lambda x: x.nucleotide_sequence, solution.complexes.values())))

      # ENUMERATION
      if scheme == 'soloveichik2010_gen.ts':# and crn == 'A + A -> A':
        print "# WARNING: changing eqivalence notion --reject-remote"
        args.REJECT_REMOTE = True

      if scheme == 'srinivas2015_gen.ts':# and crn == 'A + A -> A':
        print "# WARNING: changing eqivalence notion --reject-remote"
        args.REJECT_REMOTE = True

      print "Size of enumerated network:",
      peppercorn = TestTubePeppercornIO(solution, args)
      peppercorn.enumerate()
      #print peppercorn.condense_reactions

      # Reduce the enumerated CRN and find an interpretation
      enum_crn = peppercorn.condense_reactions
      enum_solution = peppercorn.testtube
      icrn, interpret = preprocess(fcrn, enum_crn, fs,
          solution, enum_solution, verbose=(args.verbose>1))

      print len(icrn)
      speed = len(icrn)

      # VERIFICATION
      #print "\n# Verifying Implementation CRN:"
      if len(crn) < 20 :
        v, _ = verify(fcrn, icrn, fs, interpret=interpret, method='bisimulation',
            verbose=(args.verbose>1))
        if v:
          print "{}: CRNs are bisimulation equivalent.".format(v)
        else:
          print "{}: CRNs are not bisimulation equivalent.".format(v)
        vb = v

        if scheme == 'srinivas2015_gen.ts' and crn == 'A -> A + A':
          print "# WARNING: skipping pathway equivalence"
          v = None
        elif scheme == 'cardelli2011_2D_gen.ts' and crn == '-> B':
          print "# WARNING: skipping pathway equivalence"
          v = None
        elif scheme == 'cardelli2011_2D_gen.ts' and crn == 'A -> A + A':
          print "# WARNING: skipping pathway equivalence"
          v = None
        elif scheme == 'cardelli2011_2D_gen.ts' and crn == 'A + B -> B + B':
          print "# WARNING: skipping pathway equivalence"
          v = None
        else :
          v, _ = verify(fcrn, icrn, fs, interpret=interpret, method='pathway',
              verbose=(args.verbose>1))
          if v:
            print "{}: CRNs are pathway equivalent.".format(v)
          else:
            print "{}: CRNs are not pathway equivalent.".format(v)
        vp = v
      else :
        vb = None
        vp = None

      plotdata.append([scheme, crn, vb, vp, cost, speed])


  # Results:
  for p in plotdata: print p

  # Normalize data to reference scheme/crn
  normalize = 'soloveichik2010_gen.ts'
  if normalize :
    normdata = []
    norm_values = filter(lambda x: x[0] == normalize, plotdata)
    for n in norm_values: # A particular CRN
      current = filter(lambda x: x[1] == n[1], plotdata)
      for c in current:
        normdata.append(c[:4] + [float(x)/y for x,y in zip(n[4:], c[4:])])
      
    for p in normdata: print p
    plotdata = normdata


  df = pd.DataFrame(plotdata, columns=['Translation scheme', 'CRN', 'bisimulation equivalent', 
    'pathway equivalent', 'relative counts of nucleotides', 'relative size of DSD network'])
  sns.set(style="ticks", color_codes=True)
  g = sns.PairGrid(data=df, hue='Translation scheme', size=4, 
      hue_order=['srinivas2015_gen.ts', 'soloveichik2010_gen.ts', 'qian2011_gen.ts', 
        'cardelli2011_FJ_gen.ts', 'cardelli2011_NM_gen.ts', 'cardelli2013_2D_gen.ts'],
      y_vars=["relative counts of nucleotides"], 
      x_vars=["relative size of DSD network"])

  #for ax in g.axes.flat: 
  #  ax.set_ylim(0,2)
  #  ax.set_xlim(0,2)
 
  g = g.map(plt.scatter)
  g = g.add_legend(bbox_to_anchor=(1.5, 0.5))
  #g = g.add_legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)

  pfile = 'mpp_plot.pdf'
  plt.savefig(pfile, bbox_inches="tight")


if __name__ == '__main__':
  main()
