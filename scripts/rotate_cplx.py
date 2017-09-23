#!/usr/bin/env python

import os
import sys
import argparse

import nuskell
from nuskell.objects import NuskellDomain, NuskellComplex
from dsdobjects.parser import parse_kernel_string


def get_args(parser):
    parser.add_argument("-v", "--verbose", action='count', default=0,
                        help="Print more output (-vv for extra debugging information)")
    return parser


def main(args):
    """
    The Nuskell compiler.
    """

    input_cplx = sys.stdin.readlines()
    input_cplx = "".join(input_cplx)
    print input_cplx

    ppil = parse_kernel_string(input_cplx)

    def resolve_loops(loop):
        """ Return a sequence, structure pair from kernel format with parenthesis. """
        sequen = []
        struct = []
        for dom in loop:
            if isinstance(dom, str):
                sequen.append(dom)
                if dom == '+':
                    struct.append('+')
                else:
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

    domains = {'+' : '+'} # saves some code
    complexes = dict()
    for line in ppil:
        name = line[1]
        if line[0] == 'complex':
            sequence, structure = resolve_loops(line[2])

            # Replace names with domain objects.
            try :
                sequence = map(lambda d : domains[d], sequence)
            except KeyError:
                for e, d in enumerate(sequence):
                    if d not in domains :
                        #logging.warning("Assuming {} is a long domain.".format(d))
                        domains[d] = NuskellDomain(d, dtype='long')
                        cdom = ~domains[d]
                        domains[cdom.name] = cdom
                    sequence[e] = domains[d]

            constant, concentration = False, 0
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

            complexes[name] = NuskellComplex(sequence, structure, name=name)

        else :
            raise NotImplementedError('cannot interpret keyword:', line[0])

    for cplx in complexes.values():
        for i in cplx.rotate():
            print i.name, '=', i.kernel_string

    if False:
        from nuskell.objects import TestTube, TestTubeIO
        tt = TestTube()
        tt.add_complex(cplx, (None, None))
        TestTubeIO(tt).write_dnafile(sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = get_args(parser)
    args = parser.parse_args()
    main(args)
