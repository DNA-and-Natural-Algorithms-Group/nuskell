#!/usr/bin/env python

import os
import sys
import argparse

import nuskell
from nuskell.objects import Domain, Complex
from nuskell.parser import parse_pil_string

def get_args(parser) :
  parser.add_argument("-v", "--verbose", action='count', default=0,
      help="Print more output (-vv for extra debugging information)")
  return parser

def main(args) :
  """ 
  The Nuskell compiler.  
  """

  input_cplx = sys.stdin.readlines()
  input_cplx = "".join(input_cplx)
  print input_cplx

  ppil = parse_pil_string(input_cplx)

  def resolve_loops(loop):
    """ Return a sequence, structure pair from kernel format with parenthesis. """
    sequen = []
    struct = []
    for dom in loop :
      if isinstance(dom, str):
        sequen.append(dom)
        if dom == '+' :
          struct.append('+')
        else :
          struct.append('.')
      elif isinstance(dom, list):
        struct[-1] = '('
        old = sequen[-1]
        se, ss = resolve_loops(dom)
        sequen.extend(se)
        struct.extend(ss)
        sequen.append(old + '*' if old[-1] != '*' else old[:-1])
        struct.append(')')
    return sequen, struct

  domains = []
  for line in ppil :
    if line[0] == 'domain':
      raise Exception
    elif line[0] == 'complex':
      name = line[1]
      sequence, structure = resolve_loops(line[2])
      constant, concentration = None, float('inf')
      if len(line) > 3:
        i, c, u = line[3]
        constant = (i == 'constant')
        if u == 'M':
          concentration = float(c)
        elif u == 'mM':
          concentration = float(c)*1e-3
        elif u == 'uM':
          concentration = float(c)*1e-6
        elif u == 'nM':
          concentration = float(c)*1e-9
        elif u == 'pM':
          concentration = float(c)*1e-12
        else :
          raise ValueError('unknown unit for concentrations specified.')

      for e in range(len(sequence)):
        d = sequence[e]
        if d == '+': continue
        if d[-1] == '*' : 
          dname = d[:-1]
          dom = filter(lambda x: x.name == dname, domains)
          if len(dom) < 1 :
            dom.append(Domain(name=dname, sequence = ['N']))
            domains.append(dom[0])
          elif len(dom) > 1 :
            raise RuntimeError('Conflicting matches for domain specification', d)
          sequence[e] = dom[0].get_ComplementDomain(list('R'*dom[0].length))
        else :
          dname = d
          dom = filter(lambda x: x.name == dname, domains)
          if len(dom) < 1 :
            dom.append(Domain(name=dname, sequence = ['N']))
            domains.append(dom[0])
          elif len(dom) > 1 :
            raise RuntimeError('Conflicting matches for domain specification', d)
          sequence[e] = dom[0]

  cplx = Complex(sequence = sequence, structure = structure, name=name)

  for i in cplx.rotate :
    print i.name, '=', i.kernel


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser = get_args(parser)
  args = parser.parse_args()
  main(args)

