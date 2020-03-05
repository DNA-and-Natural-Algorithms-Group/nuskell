#
#  nuskell/verifier/crn_pathway_equivalence.py
#  NuskellCompilerProject
#
# Copyright (c) 2009-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
from __future__ import absolute_import, division, print_function
from builtins import map

import logging
log = logging.getLogger(__name__)

import time
import string

from nuskell.verifier.basis_finder import find_basis
from nuskell.crnutils import genCRN, find_wastes

def pretty_crn(crn):
    for rxn in crn:
        yield '{} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]))

def printRxn(rxn, inter = {}):
    first = True
    for x in rxn[0]:
        if x[0] not in string.letters:
            if x in inter.keys() and inter[x]==[]:
                x = "w" + x
            else:
                x = "i" + x
        if not first:
            print("+",end='')
        else:
            first = False
        print(x,end='')
    print("->",end='')
    first = True
    for x in rxn[1]:
        if x[0] not in string.letters:
            if x in inter.keys() and inter[x]==[]:
                x = "w" + x
            else:
                x = "i" + x
        if not first:
            print("+",end='')
        else:
            first = False
        print(x,end='')
    print()

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

def passes_atomic_condition():
    pass

def passes_permissive_condition(fbasis, fbasis2, fbasis_raw, inter):
    # Every (formal) initial state can yield the same (formal) next state...

    def cartesian_product(l):
        if len(l) == 0:
            return []
        if len(l) == 1:
            return l[0]
        r = []
        for i in l[0]:
            for j in l[1]:
                r.append(i+j)
        return cartesian_product([r]+l[2:])

    interrev = {}
    for x in inter:
        for y in inter[x]:
            if y not in interrev:
                interrev[y] = [[x]]
            else:
                interrev[y].append([x])

    for rxn in fbasis:
        initial_states = cartesian_product(list(map(lambda x: interrev[x], rxn[0])))
        for initial in initial_states:
            initial = sorted(initial)
            flag = False
            for r in fbasis2:
                if r[0] == initial and r[1] == rxn[1]:
                    flag = True
                    break
            if not flag:
                log.info("Permissive test failed:")
                log.info("  Cannot get from {} to {}".format(initial, rxn[1]))
                log.info("Formal basis:")
                [log.info('    {}'.format(r)) for r in pretty_crn(basis)]
                return False
    return True

def passes_delimiting_condition(crn1, basis):
    # Every formal pathway must also be possible in the compiled CRN
    # Every compiled pathway must also be possible in the formal CRN
    flag = True
    for rxn in crn1:
        if rxn not in basis:
            reactants = {}
            for x in rxn[0]:
                if x in reactants:
                    reactants[x] += 1
                else:
                    reactants[x] = 1
            products = {}
            for x in rxn[1]:
                if x in products:
                    products[x] += 1
                else:
                    products[x] = 1
            if reactants != products:
                log.info("Delimiting condition failed: The formal reaction ...")
                log.info('  {} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1])))
                log.info("... is in the input CRN but not in the compiled CRN.")
                flag = False

    for rxn in basis:
        if rxn not in crn1:
            reactants = {}
            for x in rxn[0]:
                if x in reactants:
                    reactants[x] += 1
                else:
                    reactants[x] = 1
            products = {}
            for x in rxn[1]:
                if x in products:
                    products[x] += 1
                else:
                    products[x] = 1
            if reactants != products:
                log.info("Delimiting condition failed: The formal reaction ...")
                log.info('  {} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1])))
                log.info("... is in the compiled CRN but not in the input CRN.")
                flag = False

    return flag

def test(c1, c2, inter, integrated = False):
    """Test two CRNs for pathway equivalence.

    Args:
        c1: Tuple of formal CRN and formal species
        c2: Tuple of implementation CRN and formal species
        inter: partial interpretation of species
        integrated: Use integrated hybrid notion
    """
    (crn1, fs1) = c1
    (crn2, fs2) = c2
    crn1 = [[sorted(rxn[0]), sorted(rxn[1])] for rxn in crn1]
    crn2 = [[sorted(rxn[0]), sorted(rxn[1])] for rxn in crn2]
    crn2.sort()

    # Interpret waste species as nothing.
    wastes = find_wastes(crn2, fs2)
    fs2 = set(fs2) | wastes
    for x in wastes: 
        inter[x] = []

    # Remove trivial reactions.
    remove_target = []
    for [R, P] in crn2:
        if sorted(R) == sorted(P):
            remove_target.append([R, P])
    for r in remove_target:
        crn2.remove(r)

    species = set().union(*[set().union(*rxn) for rxn in crn2])

    log.info("Original CRN:")
    [log.info('    {}'.format(r)) for r in genCRN(crn1, rates = False)]
    log.info("")
    log.info("Compiled CRN with partial interpretation of ({}) species:".format(len(species)))
    [log.info('    {}'.format(r)) for r in genCRN(crn2, interpretation = inter, rates = False)]
    log.info("")

    t1 = time.time()
    log.debug('Formal species to find basis: {}'.format(fs2))
    basis = find_basis(crn2, fs2, True, inter if integrated else None)

    if basis == None: # irregular or nontidy
        log.info("Pathway equivalence: could not find formal basis.")
        return False

    if integrated: # integrated hybrid
        (fbasis_raw, fbasis) = basis 

        for i in range(len(fbasis_raw)):
            fbasis_raw[i][0].sort()
            fbasis_raw[i][1].sort()
        fbasis_raw = remove_duplicates(fbasis_raw)
        for i in range(len(fbasis)):
            fbasis[i][0].sort()
            fbasis[i][1].sort()
        fbasis = remove_duplicates(fbasis)

        # bisimulation test
        fbasis2 = []
        for [initial, final] in fbasis_raw:
            def collapse(l):
                l2 = []
                for x in l:
                    if x in inter.keys():
                        y = inter[x]
                    else:
                        y = [x]
                    l2 += y
                return l2
            r = [sorted(initial), sorted(collapse(final))]
            fbasis2.append(r)
    else: # compositional hybrid
        fbasis_raw = basis

        for i in range(len(fbasis_raw)):
            fbasis_raw[i][0].sort()
            fbasis_raw[i][1].sort()
        fbasis_raw = remove_duplicates(fbasis_raw)

        # bisimulation test
        fbasis2 = []
        fbasis = []
        for [initial, final] in fbasis_raw:
            def collapse(l):
                l2 = []
                for x in l:
                    l2 += inter.get(x, [x])
                return l2
            r = [sorted(initial), sorted(collapse(final))]
            fbasis2.append(r)
            r = [sorted(collapse(initial)), sorted(collapse(final))]
            fbasis.append(r)
        fbasis = remove_duplicates(fbasis)

    # TODO : the following is not strictly correct because it tests for
    #        strong bisimulation instead of weak bisimulation.
    if passes_permissive_condition(fbasis, fbasis2, fbasis_raw, inter):
        pass

    flag = passes_delimiting_condition(crn1, fbasis)

    log.info("Basis of the compiled CRN:")
    [log.info('    {}'.format(r)) for r in pretty_crn(fbasis)]

    t2 = time.time()
    log.debug("Elapsed time: {}".format(t2-t1))
    return flag

if __name__ == "__main__":
    import sys
    import argparse
    from nuskell import __version__
    from nuskell.crnutils import parse_crn_file, split_reversible_reactions, find_wastes, genCRN

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help="print verbose output. -vv increases verbosity level.")
    parser.add_argument("--crn-file", action='store', metavar='</path/to/file>',
            help="""Read a CRN from a file.""")
    parser.add_argument("--formal-species", nargs = '+', default = [], 
            action = 'store', metavar = '<str>', 
            help="""List formal species in the CRN.""")
    parser.add_argument("--fuel-species", nargs = '+', default = [], 
            action = 'store', metavar = '<str>', 
            help="""List fuel species in the CRN.""")
    parser.add_argument("--non-modular", action = 'store_true',
            help="""Do not optimize using CRN modules.""")
    args = parser.parse_args()

    if interactive:
        print("Enter formal species along with its interpretation:")
        print("(e.g. i187 -> A + B)")
        print("When done, press ctrl + D.")
        for line in sys.stdin:
            z = list(map(lambda x: x.strip(), line.split("->")))
            y1 = z[0]
            y2 = list(map(lambda x: x.strip(), z[1].split("+")))
            if y1[0] == "i" or y1[0] == "w": y1 = y1[1:]
            inter[y1] = y2
            fs2.add(y1)
        print()

    def remove_const(crn, const):
        for rxn in crn:
            for x in const:
                while x in rxn[0]:
                    rxn[0].remove(x)
                while x in rxn[1]:
                    rxn[1].remove(x)
        return crn

    log.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s %(message)s')
    ch.setFormatter(formatter)
    if args.verbose == 0:
        ch.setLevel(logging.WARNING)
    elif args.verbose == 1:
        ch.setLevel(logging.INFO)
    elif args.verbose == 2:
        ch.setLevel(logging.DEBUG)
    elif args.verbose >= 3:
        ch.setLevel(logging.NOTSET)
    log.addHandler(ch)

