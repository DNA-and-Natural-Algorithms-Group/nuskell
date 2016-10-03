#/usr/bin/env python

import os
import sys
import string
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pyparsing import ParseException

from nuskell.compiler import translate
from nuskell.enumeration import enumerate_crn_old
from nuskell.include.peppercorn.enumerator import get_peppercorn_args

def test_scheme_directory(input_crn, ts_dir, args, normalize = ''):
  rawdata = [] # all data, but not in Pandas format
  norm_values = []
  if normalize : # then do the respective scheme first
    crn_string = input_crn.replace('\n','; ')
    print 'Scheme:', normalize, 'CRN:', crn_string
    circuit_objects = try_to_compile(input_crn, ts_dir + normalize, args)

    if circuit_objects == []:
      raise SystemExit('Scheme for normalization failed!')
    else :
      plot_vars = [normalize, crn_string]
      plot_values = num_circuit_properties(*circuit_objects)
      # normalize everything to 1
      norm_values = plot_values
      plot_values = [float(x)/y for x,y in zip(plot_values, norm_values)]
      rawdata.append(plot_vars + plot_values)

  for scheme in os.listdir(ts_dir) :
    #if len(rawdata) > 5 : break # testing
    if scheme[-3:] != '.ts' : continue
    if scheme == normalize : continue

    crn_string = input_crn.replace('\n','; ')
    print 'Scheme:', scheme, 'CRN:', crn_string
    circuit_objects = try_to_compile(input_crn, ts_dir + scheme, args)

    plot_vars = [scheme, crn_string]
    if circuit_objects == []:
      plot_values = [None, None, None, None]
    else :
      plot_values = num_circuit_properties(*circuit_objects)
      # normalize everything to 1
      if normalize :
        plot_values = [float(x)/y for x,y in zip(plot_values, norm_values)]
      
      # Only append if it worked, but you can change the indnet to include NaNs
      rawdata.append(plot_vars + plot_values)
  return rawdata

def num_circuit_properties(domains, strands, fs, cs, ecrn):
  init = [None, None, None, None, None]

  #print '# Number of Domains:', len(domains)
  init[0] = len(domains) 
  #print '# Lengths of Domains', [d.length for d in domains]
 
  #print '# Number of Strands:', len(strands)
  init[1] = len(strands)
  #print '# Lengths of Strands:', [s.length for s in strands]

  #print '# Number of formal species:', len(fs)
  #print '# Strands in each formal complex:', [len(clx.strands) for clx in fs]

  #print '# Strands in each constant complex:', [len(clx.strands) for clx in cs]
  init[2] = max([len(clx.strands) for clx in cs])

  clen = 0
  for clx in cs :
    for snd in clx.strands :
      clen += snd.length
  #print '# Number of nucleotides in all complexes:',clen
  init[3] = clen

  #print '# Size of enum crn', len(ecrn)
  init[4] = len(ecrn)

  print clen, len(ecrn)

  # print '# Number of constant species:', len(cs)
  #init[2] = len(cs)
  #init[3] = 500
  # size of enumerated network
  # total number of nucleotides
  return init

def try_to_compile(input_crn, scheme, args):
  from nuskell.interpreter.environment import RuntimeError

  try :
    args.output = 'test'
    pil = args.output+'.pil'
    dom = args.output+'.dom'
    domains, strands, fs, cs = translate(
        input_crn, scheme, pilfile=pil, domfile=dom)

    #if scheme[-12:] == 'lakin2011.ts' :
    #  args.ignore_branch_4way = True
    #else :
    #  args.ignore_branch_4way = False

    print input_crn, scheme

    if input_crn == 'O1 <=> I1; O2 + I1 <=> I2 + O1; O3 + I2 + I1 <=> I3 + O2 + O1' and scheme[-12:] == 'lakin2011.ts' :
      return []
    elif input_crn == 'A <=> A+A; A+B -> B+B; B -> ; A + C -> ; C <=> C + C'\
        and scheme[-12:] == 'lakin2011.ts':
      return []
    #elif input_crn == 'A <=> A+A; A+B -> B+B; B -> ; A + C -> ; C <=> C + C'\
    #    and scheme[-18:] == 'soloveichik2010.ts':
    #  return []
    elif input_crn == 'A <=> A+A; A+B -> B+B; B -> ; A + C -> ; C <=> C + C'\
        and scheme[-15:] == 'srinivas2015.ts':
      return []
    else :
      enum_crn, cplxs, slow_cplxs = enumerate_crn_old(args, dom)


    return [domains, strands, fs, cs, enum_crn]
  except ParseException as e:
    print 'cannot parse the translation scheme'
    return []
  except RuntimeError as e:
    print 'cannot translate the CRN using this scheme'
    return []
  except SystemExit as e:
    print 'Scheme exits with:', e
    return []


def main():
  """Compare different Tranlation schemes for different CRNs.

  A number of descriptors are gathered in the main loop, stored in a DataFrame
  and then plotted.
  """
  parser = argparse.ArgumentParser()
  parser.add_argument("--ts_dir", action = 'store', default = "examples/ts/",
      help="Specify directory that contains translation schemes" + \
        "files that has a *.ts ending will be interpreted as translation scheme")

  parser = get_peppercorn_args(parser)

  args = parser.parse_args()

  ts_list = args.ts_dir

  # read CRN from STDIN
  # input_crn = sys.stdin.readlines()
  # input_crn = "".join(input_crn)

  #crn_list = []
  crn_list = ['A + B -> B + B',
              'A + B -> B + B; B + C -> C + C; A + C -> A + A',
              'O1 <=> I1; O2 + I1 <=> I2 + O1; O3 + I2 + I1 <=> I3 + O2 + O1',
              #'O1 <=> I1; I1 + O2 <=> O1 + I2; I1 + I2 + O3 <=> O1 + O2 + I3',
              'A <=> A+A; A+B -> B+B; B -> ; A + C -> ; C <=> C + C'
             ]
 
  crn_list.extend(['A -> ', '-> B', 'A -> B', 'A <=> B', 'A -> B + C'])

  #crn_list.extend(['A <=> B \n A -> B + C \n A+B -> X+Y \n A+B ->B+B' ])

  cols = ['Scheme', 'CRN', 
      '# of Domains', '# of Strands', 
      'max(Strand in Complex)', 
      '# of Nucleotides in all Complexes', 
      'size of enumerated CRN']

  rawdata = []
  for crn in crn_list:
    rawdata.extend(
        test_scheme_directory(
          crn, args.ts_dir, args, normalize='qian2011.ts')) 
    print 

  print rawdata

  df = pd.DataFrame(rawdata, columns=cols)
  #print df

  sns.set(style="ticks", color_codes=True)

  g = sns.PairGrid(data=df, hue='Scheme', size=4,
      y_vars=["# of Nucleotides in all Complexes"], 
      x_vars=["size of enumerated CRN"])

  #g = sns.PairGrid(data=df, hue='Scheme', size=4,
  #    x_vars=["max(Strand in Complex)"], 
  #    y_vars=["# of Nucleotides in all Complexes"])

  #for ax in g.axes.flat: 
  #  ax.set_ylim(0,2)
  #  ax.set_xlim(0,2)
 
  g = g.map(plt.scatter)
  g = g.add_legend()

  pfile = 'compare.pdf'
  plt.savefig(pfile)

  # rawdata = []
  # for crn in crn_list:
  #   rawdata.extend(
  #       test_scheme_directory(crn, args.ts_dir, args))

  # print r for r in rawdata


if __name__ == '__main__':
  main()
