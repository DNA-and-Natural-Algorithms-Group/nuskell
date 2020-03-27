#!/usr/bin/env python
#
#  nuskell/verifier/verifier.py
#  NuskellCompilerProject
#
# Copyright (c) 2009-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
from __future__ import absolute_import, division, print_function
from builtins import map, filter

import logging
log = logging.getLogger(__name__)

from collections import Counter

from nuskell.crnutils import genCRN

class NoFormalBasisError(Exception):
    pass

class BasisFinderError(Exception):
    pass

def pretty_crn(crn):
    for rxn in crn:
        yield '{} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]))

def collapse(l):
    global inter
    l2 = []
    for x in l:
        l2 += inter.get(x, [x])
    return l2

def formal(s, fs):
  """Returns all formal species in the initial state S. """
  return list(filter(lambda x: x in fs, s))

def intermediate(s, fs):
  """Returns all non-formal species in the initial state S. """
  return list(filter(lambda x: x not in fs, s))

def next_state(s, rxn):
    """Applies reaction rxn to state s """
    s = s[:]
    for x in rxn[0]:
        s.remove(x)
    s += rxn[1]
    return sorted(s)

def minimal_initial_state(pathway):
    initial = []
    current = []
    for rxn in pathway:
        for r in rxn[0]:
            if r in current:
                current.remove(r) 
            else:
                initial.append(r)
        for r in rxn[1]:
            current.append(r)
    return sorted(initial)

def contained(a, b):
    """True if list a is contained in list b, False otherwise.

    Args: 
      a (List[str]): A list of elements, eg. reactands
      b (List[str]): A list of elements, eg. an initial state

    Returns:
      True/False
    """
    def contained_helper(a, b):
      if len(a) == 0:
        return True
      if len(b) < len(a):
        return False
      if a[0] == b[0]:
        return contained_helper(a[1:], b[1:])
      else:
        return contained_helper(a, b[1:])

    a = sorted(a)
    b = sorted(b)
    return contained_helper(a, b)

def final_state(p, S):
    c = S[:]
    for rxn in p:
        for x in rxn[0]:
            if x not in c: 
                return None
            c.remove(x)
        for x in rxn[1]:
            c.append(x)
    return sorted(c)

def linearcheck(path, S, nonw):
    init = S[:]
    for rxn in path:
        c1 = list(filter(lambda x: x in nonw, init))
        if len(c1) > 1: 
            return False

        for r in rxn[0]:
            assert r in init
            init.remove(r)
        for p in rxn[1]:
            init.append(p)

    c1 = list(filter(lambda x: x in nonw, init))
    if len(c1)>1: 
        return False
    return True

def setminus(a, b):
    a = a[:]
    b = b[:]
    T = []
    for x in a:
        if x in b:
            b.remove(x)
        else:
            T.append(x)
    return T

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

def formal_state(S, fs):
    return len(formal(S, fs)) == len(S)

def decompose(path, fs):
    def helper(p1, p2, remaining, fs):
        i1 = minimal_initial_state(p1)
        i2 = minimal_initial_state(p2)
        if not formal_state(i1, fs) or not formal_state(i2, fs):
            return []
        if len(remaining) == 0:
            if len(p1) > 0 and len(p2) > 0:
                # valid decomposition
                return [[final_state(p1, i1), final_state(p2, i2)]]
            else:
                return []
        else:
            return helper(p1 + [remaining[0]], p2, remaining[1:], fs) \
                 + helper(p1, p2 + [remaining[0]], remaining[1:], fs)
    l = remove_duplicates(sorted(helper([], [], path, fs)))
    return l

def formal_closure(p, fs):
    initial = minimal_initial_state(p)
    S = initial[:]
    T = formal(S, fs)
    assert len(S) == len(T) # STB: only defined for semiformal paths!
    for r in p:
        S = next_state(S, r)
        f = formal(S, fs)
        T = T + setminus(f, T)
    return sorted(T)

def regular_final_state(p, fs):
    initial = minimal_initial_state(p)
    S = initial[:]
    state_list = [S]
    flag = True
    T = []
    RFS = []
    j = 0
    for r in p:
        S = next_state(S, r)
        state_list.append(S)
        f = formal(S, fs)
        if not contained(f, initial):
            flag = False
        if flag: j += 1
    i = 0
    for r in reversed(p):
        i += 1
        (R,P) = r
        f = formal(state_list[-i], fs)
        T = T + setminus(f, T)
        if j >= len(p) - i and \
           len(formal(setminus(state_list[-(i+1)],R),fs)) == 0:
            T.sort()
            if T not in RFS:
                RFS.append(T)
    return sorted(RFS)

def width(pathway):
    S = minimal_initial_state(pathway)
    w = len(S)
    for rxn in pathway:
        S = next_state(S, rxn)
        if len(S) > w: w = len(S)
    return w

def enumerate_pathways(p, w_max, i_max, crn, fs, nonw = None):
    """Enumerate every path reaction in the CRN 
    Args:
        p = path
    """
    global inter # Integrated hybrid stuff ...
    global ret, ccheck, ebasis

    w_p = width(p)
    if w_p > w_max: 
        log.debug("Return: Path reached maximum width: {} > {}".format(w_p, w_max))
        return

    initial = minimal_initial_state(p)
    if not formal_state(initial, fs): 
        log.debug("Return: Not a formal initial state: {}".format(initial))
        return

    if len(initial) > i_max: 
        log.debug("Return: Initial state too big: {} > {}".format(len(initial), i_max))
        return

    if nonw is not None and not linearcheck(p, initial, nonw): 
        log.debug("Return: Linearcheck failed: init = {}, path = {}".format(initial, p))
        return

    final = final_state(p, initial)
    DFS = decompose(p, fs)
    for [p1, p2] in DFS:
        if formal_state(p1, fs) or formal_state(p2, fs):
            log.debug("Return: strong decomposable: init = {}, path = {}".format(initial, p))
            return
    fc = formal_closure(p, fs)
    RFS = regular_final_state(p, fs)
    sig = (initial, final, w_p, fc, DFS, RFS)

    if sig in ret: 
        return
    ret.append(sig)

    if len(DFS) == 0:
        if final not in ccheck:
            ccheck.append(final)
            if not tidy(final, crn, fs):
                log.info("The given system is not tidy.")
                log.info(" - from initial state: {}".format(initial))
                [log.info('    ' + r) for r in pretty_crn(p)]
                raise NoFormalBasisError("The given system is not tidy.")
        if len(p) > 0 and formal_state(final, fs):
            if inter: # integrated hybrid theory
                p1 = list(map(lambda rxn: [collapse(rxn[0]), collapse(rxn[1])], p))
                fsp = collapse(fs)
                initial1 = minimal_initial_state(p1)
                final1 = final_state(p1, initial1)
                RFS1 = regular_final_state(p1, fsp)
                if final1 not in RFS1:
                    log.info("The given system is not regular.")
                    log.info(" - from initial state: {}".format(initial))
                    [log.info('    ' + r) for r in pretty_crn(p)]
                    raise NoFormalBasisError("The given system is not regular.")
            elif final not in RFS:
                log.info("The given system is not regular.")
                log.info(" - from initial state: {}".format(initial))
                [log.info('    ' + r) for r in pretty_crn(p)]
                raise NoFormalBasisError("The given system is not regular.")
            ebasis.append(p)

    perm = False # turn on if you want a randomized/shuffled search
    if perm:
        import copy, random
        crn_perm = copy.copy(crn) 
        random.shuffle(crn_perm)
        tmp = crn_perm
    else :
        tmp = crn
    for r in crn:
        enumerate_pathways(p + [r], w_max, i_max, crn, fs, nonw = nonw)
  
def tidy(S, crn, fs):
    """BFS to test if the crn is *strongly(?)* tidy.

    Let p be a pathway with a formal initial state and T its final state. Then, a
    pathways p' is said to be a closing pathway of p if p' can occor in T and
    results in a formal state. A closing pathway is strong if its reactios do not
    consume any formal species.

    A CRN is weakly tidy if every pathway with a formal initial state has a
    closing pathway. A CRN is strongly tidy if every pathway with a formal
    initial state has a strong closing pathway.

    Args:
      S (List[str]) : Initial state S
      crn (LoL[[[reactands],[products]],...]): a list of irreversible reactions
      fs (List[str]) : A list of formal states

    # STW2016 - T3.5: Suppose C is a tidy and regular CRN and F is its formal
    # basis.  1) For any formal pathway q in F, there exists a formal pathway p in
    # C whose initial and final states are equal to those of q, such that p can be
    # interpreted as q.  2) Any formal pathway p in C can be interpreted as some
    # pathway q in F.
    # -- appearently C doesn't need to be tidy here ...

    """
    # BFS. May have to eventually rewrite to impose a width bound.
    queue = [intermediate(S, fs)]
    mem = [queue[0]]
    while len(queue) > 0:
        S = queue[0]
        queue = queue[1:]
        if len(S) == 0: 
            return True
        for rxn in crn:
            if len(rxn[0]) > 0 and contained(rxn[0], S):
                nS = intermediate(next_state(S, rxn), fs)
                if nS not in mem:
                    mem.append(nS)
                    queue.append(nS)
    return False

signatures = []
def enumerate_basis(crn, fs, nonw = None):
    """
    The width is the maximum number of species present at any timepoint during
    the reactions in a pathway.

    Args:
        crn: The CRN.
        fs: Formal species (and waste species).
        nonw: Reactants that are not formal species ("non-wastes') and which
            are supposed to pass the linearcheck: There is at most one
            non-waste in reactants and products for every reaction.
    """

    # Different branching factors are for optimizations 
    #   (and explained in 5.6 of SWS's Master's thesis)
    # b: the branching factor
    # bf: max over |intermediate(R)| of (R,P) in crn
    # br = [] # set of (|formal(R)|,|intermediate(R)|)
    b  = max(max(len(r),len(p)) for [r,p] in crn)
    bf = max(len(intermediate(r,fs)) for [r,p] in crn)
    br = [[len(formal(r,fs)), len(intermediate(r,fs))] for [r,p] in crn]
    br = remove_duplicates(br)
    log.debug('Branching Factors: b = {}, bf = {}, br = {}'.format(b, bf, br))

    global ret, ccheck, ebasis
    ccheck = [] # A cache to look-up which final states are tidy
    w_max, i_max = 0, 0
    while True:
        log.debug("Current bounds: w_max = {} i_max = {}".format(w_max, i_max))
        ret = []
        ebasis = []

        # If successfull, returns variables: 
        #  - ebasis: enumerated bases (if sucessfull)
        #  - ret: A list of signatures
        enumerate_pathways([], w_max, i_max, crn, fs, nonw = nonw)
        log.debug('After enumation of pathways: basis = {}'.format(ebasis))

        signatures = ret
        current_w = 0
        current_i = 0
        for (i, f, w, fc, dfs, rfs) in signatures:
            log.debug('Sig: {}'.format([i, f, w, fc, dfs, rfs]))
            if len(dfs) == 0:
                if w > current_w:
                    current_w = w
                if len(i) > current_i:
                    current_i = len(i)

        w_t = current_w * bf + b 
        i_t = 0
        for [x, y] in br:
            i_t = max(i_t, current_i * y + x)
        if w_t <= w_max and i_t <= i_max: 
            break
        w_max = w_t
        i_max = i_t
    
    fbasis_raw = []
    fbasis_int = [] # interpretation
    for p in ebasis:
        initial = minimal_initial_state(p)
        final = final_state(p, initial)
        r = [sorted(initial), sorted(final)]
        fbasis_raw.append(r)

        if inter: # integrated hybrid theory
            p1 = list(map(lambda rxn: [collapse(rxn[0]), collapse(rxn[1])], p))
            initial = minimal_initial_state(p1)
            final = final_state(p1, initial)
            r = [sorted(initial), sorted(final)]
            fbasis_int.append(r)

    fbasis_int = remove_duplicates(sorted(fbasis_int))
    fbasis_raw = remove_duplicates(sorted(fbasis_raw))
    return fbasis_raw, fbasis_int

def find_modules(crn, intermediates):
    # STW: Namely, we can partition the CRN into disjoint subsets that do
    # not share intermediate species with one another, find the formal basis
    # of each subset and then take the union of the found fomral bases.
    def ancestor(x):
        if parent[x] != x:
            parent[x] = ancestor(parent[x])
        return parent[x]

    # Get all the divisions
    parent   = {i:i  for i in intermediates}
    division = {i:[] for i in intermediates}
    # First, determine which intermediates occur together...
    for rxn in crn:
        t = list(filter(lambda x: x in intermediates, rxn[0] + rxn[1]))
        if len(t) > 1:
            # if you got multiple intermediates, parent should point 
            # to the same "ancestor"
            z = ancestor(t[0])
            for x in t[1:]:
                y = ancestor(x)
                if z != y:
                    parent[y] = z
    i = 0
    for rxn in crn:
        t = list(filter(lambda x: x in intermediates, rxn[0] + rxn[1]))
        if len(t) > 0:
            z = ancestor(t[0])
            division[z].append(rxn)
        else: # formals, wastes only
            division[i] = [rxn]
            i += 1

    return list(filter(lambda x: len(x) > 0, division.values()))


def find_basis(crn, fs, optimize = True, interpretation = None):
    """Finds all formal reactions in a CRN.
    
    STW_16 - Def12: The set of prime pathways in a given CRN is called elementary
    basis.  The fromal basis is the set of (initial state, final state) pairs of
    the pathways in the elementatry basis.

    Args:
      crn (List[Lists]): a CRN as returned from the nuskell.parser module
      fs (List[str]): a list of formal species in the CRN
      optimize (bool): Chop the CRN into modules and then find the basis 
                       separately for each of those modules.
      interpretation (?): A second interpretation crn, *only* used in the "integrated
                  hybrid" approach

    Returns:
        basis_raw, basis_int: The formal basis found without and with an 
            interpretation dictionary. If the dictionary is None, 

    """
    global inter 
    inter = interpretation if interpretation is not None else dict()

    if optimize:
        # Optimization needs to identify non-waste species 
        log.debug("Optimization: Identifying modules in the implementation CRN.")
        species = set().union(*[set().union(*rxn) for rxn in crn])
        intermediates = species - fs

        divs = find_modules(crn, intermediates)
        log.info("Divided the implementation CRN into {} modules.".format(len(divs)))

        basis_raw = []
        basis_int = []
        for e, mod in enumerate(divs, 1):
            log.info("Verifying module {}:".format(e))
            [log.info('    {}'.format(r)) for r in genCRN(mod, rates = False)]
            # Optimization using linear structure
            reactants = set().union(*[set(r) for [r, p] in mod])
            wastes = intermediates - reactants 
            nonwastes = intermediates - wastes
            log.debug('Nonwastes (non-formal reactants): {}'.format(nonwastes))

            linear = True
            for [r,p] in mod:
                r1 = list(filter(lambda x: x in nonwastes, r))
                p1 = list(filter(lambda x: x in nonwastes, p))
                if len(r1) > 1 or len(p1) > 1: 
                    linear = False
                    break
            log.debug('Found monomolecular substructre: {}'.format(linear))
            b_raw, b_int = enumerate_basis(mod, fs, nonw = nonwastes if linear else None)

            log.debug("Formal basis of the current module:")
            [log.debug('    {}'.format(r)) for r in pretty_crn(b_raw)]
            log.debug("Formal basis of the current module (interpreted):")
            [log.debug('    {}'.format(r)) for r in pretty_crn(b_int)]
            log.info('')
            basis_raw += b_raw
            basis_int += b_int
    else:
        basis_raw, basis_int = enumerate_basis(crn, fs)

    return basis_raw, basis_int

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
    parser.add_argument("--integrated", action = 'store_true',
            help="""Use interpretation when finding formal basis.""")
    args = parser.parse_args()

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

    crn, fs = parse_crn_file(args.crn_file)
    crn = split_reversible_reactions(crn)
      
    crn = [rxn[:2] for rxn in crn]

    if args.fuel_species:
        crn = remove_const(crn, args.fuel_species)

    fs = set(fs.keys())
    if args.formal_species:
        fs &= set(args.formal_species)

    log.info('Input formal species: {}'.format(fs))
    log.info('Input fuel species: {}'.format(args.fuel_species))
    log.info('Input CRN (without fuel species):')
    [log.info('    ' + r) for r in pretty_crn(crn)]

    wastes = find_wastes(crn, fs)
    log.info('{} waste species are treated as formal. ({})'.format(len(wastes), ' '.join(wastes)))

    linter = {} # local inter (do not use global inter)
    for x in fs: linter[x] = [x]
    for x in wastes: linter[x] = []

    # TODO: Waste species have to be part of the formal species, otherwise the
    # algorithm cannot find the formal base.
    try:
        basis_raw, basis_int = find_basis(crn, fs | wastes, 
                                optimize = not args.non_modular,
                                interpretation = linter if args.integrated else None)
    except NoFormalBasisError as err:
        print("Could not find formal basis: {}".format(err))
        raise SystemExit

    print("Formal basis:")
    for r in pretty_crn(basis_raw):
        print(r)
    print("Interpretation of formal basis:")
    for r in genCRN(basis_raw, rates = False, interpretation = linter):
        print(r)
    if args.integrated:
        print("Integrated hybrid basis:")
        for r in pretty_crn(basis_int):
            print(r)
