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

from collections import Counter

from nuskell.verifier.basis_finder import find_basis, NoFormalBasisError
from nuskell.crnutils import genCRN, assign_crn_species

def pretty_crn(crn):
    for rxn in crn:
        yield '{} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]))

def passes_atomic_condition():
    # JDW (2019):
    # For every formal species A, there exists an implementation species a
    # s.t. m(a) = {| A |}
    pass

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

def passes_delimiting_condition(fcrn, icrn):
    """Checks delimiting condition. 

    Note: requires sorted reactants and products (see cleanup function).
    """
    # JDW (2019): 
    # The interpretation of any implementation reaction is either trivial or a
    # valid formal reaction.
    def trivial(rxn):
        return Counter(rxn[0]) == Counter(rxn[1])

    # Every formal reaction must also be possible in the implementation CRN
    # Every implementation reaction must also be possible in the formal CRN
    flag = True
    for rxn in fcrn:
        if rxn not in icrn:
            if not trivial(rxn):
                log.info("Delimiting condition failed: The formal reaction ...")
                log.info('  {} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1])))
                log.info("... is in the input CRN but not in the compiled CRN.")
                flag = False

    for rxn in icrn:
        if rxn not in fcrn:
            if not trivial(rxn):
                log.info("Delimiting condition failed: The formal reaction ...")
                log.info('  {} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1])))
                log.info("... is in the compiled CRN but not in the input CRN.")
                flag = False

    return flag

def passes_permissive_condition(fcrn, icrn, inter):
    """Checks permissive condition. 

    Note: requires sorted reactants and products (see cleanup function).

    """
    # JDW (2019):
    # If there exists a formal reaction r out of formal state S and an interpretation m(S') = S:
    # there exists an implementation reaction r' s.t. m(r') = r and r' must be possible in S'.

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

    for rxn in fcrn:
        initial_states = cartesian_product(list(map(lambda x: interrev[x], rxn[0])))
        for initial in initial_states:
            initial = sorted(initial)
            for [iR, iP] in icrn:
                if iR == initial and iP == rxn[1]:
                    break
            else: # didn't find a matching reaction ...
                log.info("Permissive test failed:")
                log.info("  Cannot get from {} to {}".format(initial, rxn[1]))
                log.info("Formal basis:")
                [log.info('    {}'.format(r)) for r in pretty_crn(basis)]
                return False
    return True

def cleanup(basis):
    for i in range(len(basis)):
        basis[i][0].sort()
        basis[i][1].sort()
    return remove_duplicates(basis)

def test(c1, c2, inter, integrated = False):
    """Test two CRNs for pathway equivalence.

    Args:
        c1: Tuple of formal CRN and formal species
        c2: Tuple of implementation CRN and species corresponding to formal species
        inter: partial interpretation of species
        integrated (bool, optional): Use integrated hybrid notion (True) 
                    or compositional hybrid notion (False). Defaults to False.

    """
    assert isinstance(inter, dict)

    (crn1, fs1) = c1
    (crn2, fs2) = c2
    crn1 = [[sorted(rxn[0]), sorted(rxn[1])] for rxn in crn1]
    crn2 = [[sorted(rxn[0]), sorted(rxn[1])] for rxn in crn2]
    crn2.sort()

    fs2 = set(fs2)
    # Interpret waste species as nothing.
    i, wastes, nw = assign_crn_species(crn2, fs2)
    fs2 = fs2 | wastes
    for x in wastes: 
        inter[x] = []

    # Remove trivial reactions.
    remove_target = []
    for [R, P] in crn2:
        if sorted(R) == sorted(P):
            remove_target.append([R, P])
    for r in remove_target:
        crn2.remove(r)


    log.info("Original CRN:")
    [log.info('    {}'.format(r)) for r in genCRN(crn1, rates = False)]
    log.info("")
    log.info(f"Compiled CRN with {len(set().union(*[set().union(*rxn) for rxn in crn2]))} species and {len(crn2)} reactions. After partial interpretation:")
    [log.info(f'    {r}') for r in genCRN(crn2, interpretation = inter, rates = False)]
    log.info("")

    try:
        log.debug('Formal species to find basis: {}'.format(fs2))
        fbasis_raw, fbasis_int = find_basis(crn2, fs2, 
                                           modular = True,
                                           interpretation = inter if integrated else None)
        fbasis_raw = cleanup(fbasis_raw)
        fbasis_int = cleanup(fbasis_int)
    except NoFormalBasisError as err:
        log.info("Could not find formal basis: {}".format(err))
        return False

    def collapse(l):
        l2 = []
        for x in l:
            l2 += inter.get(x, [x])
        return l2
 
    if integrated: # integrated hybrid
        fbasis2 = []
        for [initial, final] in fbasis_raw:
            r = [sorted(initial), sorted(collapse(final))]
            fbasis2.append(r)
    else: # compositional hybrid
        assert fbasis_int == []
        fbasis2 = []
        for [initial, final] in fbasis_raw:
            r = [sorted(initial), sorted(collapse(final))]
            fbasis2.append(r)
            r = [sorted(collapse(initial)), sorted(collapse(final))]
            fbasis_int.append(r)
        fbasis_int = remove_duplicates(fbasis_int)

    #
    # TODO: the following tests are for strong bisimulation instead of weak bisimulation.
    #
    if not passes_delimiting_condition(crn1, fbasis_int):
        return False

    if not passes_permissive_condition(fbasis_int, fbasis2, inter):
        return False

    return True

if __name__ == "__main__":
    import sys
    import argparse
    from nuskell import __version__
    from nuskell.crnutils import parse_crn_file, split_reversible_reactions, genCRN

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

