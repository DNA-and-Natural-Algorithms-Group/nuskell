#
#  nuskell/verifier/crn_pathway_equivalence.py
#  NuskellCompilerProject
#
# Copyright (c) 2009-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
from __future__ import absolute_import, division, print_function
#from builtins import map

import logging
log = logging.getLogger(__name__)

import sys
import time
import string

from nuskell.verifier.basis_finder import find_basis

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

def findWastes(crn, formal):
    """Returns waste species of a CRN.

    Waste species are all non-formal species that are only products, but never
    react. A non-waste is a formal species, or a species that is involved as a
    reactant in a reaction that involves a non-waste species.
    """
    species = set().union(*[set().union(*rxn) for rxn in crn])
    nonwastes = set(formal)
    while True:
        # Add x to non-waste if in any reaction x is an reactant while there
        # are non-wastes in the reaction. Reset the outer loop if you found a
        # new non-waste species.
        flag=False
        for x in species:
            if x in nonwastes: 
                continue
            for rxn in crn:
                if x in rxn[0] and (len(nonwastes & set(rxn[1]+rxn[0])) > 0):
                    nonwastes.add(x)
                    flag = True
                    break
        if not flag: break
    return species - nonwastes

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

def test(c1, c2, inter, integrated = False, interactive = False, verbose = False):
    """Test two CRNs for pathway equivalence.

    Args:
        c1: Tuple of formal CRN and formal species
        c2: Tuple of implementation CRN and formal species
        inter: partial interpretation of species
        integrated: Use integrated hybrid notion
        verbose: Print verbose information.
    
    """
    (crn1, fs1) = c1
    (crn2, fs2) = c2
    crn1 = [[sorted(rxn[0]), sorted(rxn[1])] for rxn in crn1]
    crn2 = [[sorted(rxn[0]), sorted(rxn[1])] for rxn in crn2]
    crn2.sort()

    # Interpret waste species as nothing.
    wastes = findWastes(crn2, fs2)
    fs2 = (set(fs2)).union(wastes)
    for x in wastes: inter[x] = []

    # Remove trivial reactions
    remove_target = []
    for [R, P] in crn2:
        if sorted(R) == sorted(P):
            remove_target.append([R, P])
    for r in remove_target:
        crn2.remove(r)

    species = set().union(*[set().union(*rxn) for rxn in crn2])
    if verbose :
        print("# Number of species :", len(species))
        print("Original CRN:")
        for [R,P] in crn1:
            print("   {} -> {}".format(' + '.join(R),' + '.join(P)))
        print()
        print("Compiled CRN with partial interpretation of species:")
        print(inter)
        for [R,P] in crn2:
            R = ['({})'.format(','.join(inter.get(r, [r]))) for r in R]
            P = ['({})'.format(','.join(inter.get(p, [p]))) for p in P]
            print("   {} -> {}".format(' + '.join(R),' + '.join(P)))
        print()

    if interactive:
        print("Enter formal species along with its interpretation:")
        print("(e.g. i187 -> A + B)")
        print("When done, press ctrl + D.")
        for line in sys.stdin:
            z = map(lambda x: x.strip(), line.split("->"))
            y1 = z[0]
            y2 = map(lambda x: x.strip(), z[1].split("+"))
            if y1[0] == "i" or y1[0] == "w": y1 = y1[1:]
            inter[y1] = y2
            fs2.add(y1)
        print

    t1 = time.time()
    basis = find_basis(crn2, fs2, True, inter if integrated else None)

    if basis == None: # irregular or nontidy
        print("Pathway equivalence: could not find formal basis.")
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
                    if x in inter.keys():
                        y = inter[x]
                    else:
                        y = [x]
                    l2 += y
                return l2
            r = [sorted(initial), sorted(collapse(final))]
            fbasis2.append(r)
            r = [sorted(collapse(initial)), sorted(collapse(final))]
            fbasis.append(r)
        fbasis = remove_duplicates(fbasis)

    # TODO : the following is not strictly correct because it tests for
    #        strong bisimulation instead of weak bisimulation.
    # permissive test
    interrev = {}
    for x in inter.keys():
        for y in inter[x]:
            if y not in interrev.keys():
                interrev[y] = [[x]]
            else:
                interrev[y].append([x])
    for rxn in fbasis:
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
        initial_states = cartesian_product(map(lambda x: interrev[x], rxn[0]))
        for initial in initial_states:
            initial = sorted(initial)
            flag = False
            for r in fbasis2:
                if r[0] == initial and r[1] == rxn[1]:
                    flag = True
                    break
            if not flag:
                print("Permissive test failed:")
                print("  Cannot get from ",initial," to",rxn[1])
                print("Formal basis found was:")
                for [r,p] in fbasis_raw:
                    print(r, "->", p)
                return False
    # permissive test end
    basis = fbasis

    log.info("Basis of the compiled CRN:")
    [log.info('    ' + r) for r in pretty_crn(basis)]

    # delimiting test
    flag = True
    for rxn in crn1:
        if rxn not in basis:
            reactants = {}
            for x in rxn[0]:
                if x in reactants.keys():
                    reactants[x] += 1
                else:
                    reactants[x] = 1
            products = {}
            for x in rxn[1]:
                if x in products.keys():
                    products[x] += 1
                else:
                    products[x] = 1
            if reactants != products:
              if verbose :
                print("Error : The formal pathway")
                print("    ",end='')
                printRxn(rxn)
                print(" is in the input CRN but not in the compiled CRN.")
              flag = False
    for rxn in basis:
        if rxn not in crn1:
            reactants = {}
            for x in rxn[0]:
                if x in reactants.keys():
                    reactants[x] += 1
                else:
                    reactants[x] = 1
            products = {}
            for x in rxn[1]:
                if x in products.keys():
                    products[x] += 1
                else:
                    products[x] = 1
            if reactants != products:
              if verbose :
                print("Error : The formal pathway")
                print("    ",end='')
                printRxn(rxn)
                print(" is in the compiled CRN but not in the input CRN.")
              flag = False


    t2 = time.time()
    if verbose :
      print("Elapsed time :", t2-t1)
    return flag
