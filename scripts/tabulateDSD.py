#/usr/bin/env python

from os import listdir
import sys
import argparse
import nuskell
from pyparsing import ParseException
from nuskell.interpreter.environment import RuntimeError

def builtin_schemes():
  schemes = ({ 
    soloveichick2010 : 'share/ts/soloveichick.ts',
    lakin2011: 'share/ts/lakin.ts',
    })

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--ts_dir", action = 'store', default = "$PEPPERHOME", 
      help="Specify directory that contains translation schemes" + \
        "files that has a *.ts ending will be interpreted as translation scheme")
  args = parser.parse_args()

  ts_list = args.ts_dir

  # read CRN from STDIN
  input_crn = sys.stdin.readlines()
  input_crn = "".join(input_crn)

  for scheme in listdir(args.ts_dir) :
    if scheme[-3:] != '.ts' : continue
    print 'Scheme:', scheme

    try :
      domains, strands, fs, cs = nuskell.compile(input_crn, args.ts_dir+scheme)
    except ParseException as e:
      print 'cannot parse the translation scheme'
    #except RuntimeError as e:
    #  print 'cannot translate the CRN using this scheme'
    except SystemExit as e:
      print 'Scheme exits with:', e

    else : # It worked? ... print the numbers:
      print '# Number of Domains:', len(domains)
      print '# Lengths of Domains', [d.length for d in domains]
      print '# Number of Strands:', len(strands)
      print '# Lengths of Strands', [s.length for s in strands]
      print '# Number of formal species:', len(fs)
      print '# Strands in each formal complex:', [len(clx.strands) for clx in fs]
      print '# Number of constant species:', len(cs)
      print '# Strands in each constant complex:', [len(clx.strands) for clx in cs]

    print 

if __name__ == '__main__':
  main()
