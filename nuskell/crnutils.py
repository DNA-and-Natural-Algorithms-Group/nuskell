#
#  nuskell/crnutils.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, print_function, division
from builtins import map, dict

import logging
log = logging.getLogger(__name__)

from collections import namedtuple, Counter
from crnsimulator import parse_crn_string as pcs
from crnsimulator import parse_crn_file as pcf
from crnsimulator.crn_parser import CRNParseError
from dsdobjects.utils import natural_sort

Reaction = namedtuple('reaction', 'reactants products k_fwd k_rev')

def parse_crn_string(crnstring):
    # A wrapper for the crnsimulator version
    return post_process(pcs(crnstring, process = False), defaultrate = 1, defaultmode = None, defaultconc = None)

def parse_crn_file(crnfile):
    # A wrapper for the crnsimulator version
    return post_process(pcf(crnfile, process = False), defaultrate = 1, defaultmode = None, defaultconc = None)

def post_process(crn, defaultrate = None, defaultmode = 'initial', defaultconc = 0):
    """Process a parsed CRN.
    Taken and modified from crnsimulator, drop for Python 3.x
    """
    def flint(inp):
        return int(inp) if float(inp) == int(float(inp)) else float(inp)

    def remove_multipliers(species):
        flat = []
        for s in species:
            if len(s) == 1:
                flat.append(s[0])
            elif len(s) == 2:
                ss = [s[1]] * int(s[0])
                flat.extend(ss)
        return flat

    new = []
    species = dict()
    for line in crn:
        if line[0] == 'concentration':
            spe = line[1][0]
            ini = 'initial' if line[2][0][0] == 'i' else 'constant'
            num = line[3][0]
            species[spe] = (ini, flint(num))
            continue
        elif len(line) == 3:
            # No rate specified
            t, r, p = line
            r = remove_multipliers(r)
            p = remove_multipliers(p)
            if t == 'reversible':
                new.append(Reaction(r, p, defaultrate, defaultrate))
            elif t == 'irreversible':
                new.append(Reaction(r, p, defaultrate, 0))
            else:
                raise CRNParseError('Wrong CRN format!')
        elif len(line) == 4:
            t, r, p, k = line
            r = remove_multipliers(r)
            p = remove_multipliers(p)
            if t == 'reversible':
                assert len(k) == 2
                new.append(Reaction(r, p, flint(k[0]), flint(k[1])))
            elif t == 'irreversible':
                assert len(k) == 1
                new.append(Reaction(r, p, flint(k[0]), 0))
            else:
                raise CRNParseError('Wrong CRN format!')
        else:
            raise CRNParseError('Wrong CRN format!')
        for s in r + p:
            if s not in species:
                species[s] = (defaultmode, defaultconc)
    return new, species

def crn_to_standard(crn):
    """
    At the moment, nuskell is not consistent about what format to use for
    representing a reaction. This is supposed to change, but in the meantime,
    if you want to stick to your standard, use this function for conversion
    *before* using crnutils, and then convert back if you must ...
    """
    new = []
    for rxn in crn:
        if isinstance(rxn, Reaction):
            nrxn = rxn
        else:
            log.debug('Reaction: {}'.format(rxn))
            # Unify reactants and products.
            if isinstance(rxn[0], Counter):
                assert isinstance(rxn[1], Counter)
                reactants = rxn[0].elements()
                products = rxn[1].elements()
            else:
                assert isinstance(rxn[0], list)
                assert isinstance(rxn[1], list)
                reactants = rxn[0]
                products = rxn[1]

            # Unify reaction rates.
            if len(rxn) == 2: # irreversible, no rates specified.
                nrxn = Reaction(reactants, products, 1, 0)
            elif len(rxn) == 3: # rates specified
                if isinstance(rxn[2], list) and len(rxn[2]) == 2: # reversible
                    if rxn[2][0] is None:
                        kf = 1
                    else:
                        assert isinstance(rxn[2][0], (float, int)) and rxn[2][0]
                        kf = rxn[2][0] 
                    if rxn[2][1] is None:
                        kr = 1
                    else:
                        assert isinstance(rxn[2][1], (float, int))
                        kr = rxn[2][1] 
                elif isinstance(rxn[2], list) and len(rxn[2]) == 1: # irreversible
                    if rxn[2][0] is None:
                        kf = 1
                    else:
                        assert isinstance(rxn[2][0], (float, int)) and rxn[2][0]
                        kf = rxn[2][0] 
                    kr = 0
                else:
                    if rxn[2] is None:
                        kf = 1
                    else:
                        assert isinstance(rxn[2], (float, int)) and rxn[2]
                        kf = rxn[2]
                    kr = 0
                nrxn = Reaction(reactants, products, kf, kr)
        new.append(nrxn)
    return new

def genCON(species):
    for sp, data in natural_sort(species.items()):
        if None in data: continue
        yield '{} @{} {}'.format(sp, 'initial' if data[0][0] == 'i' else 'constant', data[1])

def genCRN(crn, reversible = True, rates = True, interpretation = None):
    """ Pretty printing of CRNs.

    Args:
        crn (list of reactiontuples): A CRN in list of list format.
        reversible (bool, optional): True to combine reversible reactions into one
            line.  False to split reversible reactions into irreversible reactions.
        rates (bool, optional): 

    Note:
      Order of reactants/products may differ from the order in the crn argument.
    """
    if reversible:
        try:
            pcrn = combine_reversible_reactions(crn)
        except ValueError as err:
            crn = crn_to_standard(crn)
            pcrn = combine_reversible_reactions(crn)
    else:
        pcrn = split_reversible_reactions(crn)

    def inter(l):
        l2 = []
        for x in l:
            l2 += interpretation.get(x, [x])
        return l2

    if interpretation:
        icrn = []
        for rxn in pcrn:
            R = inter(rxn.reactants)
            P = inter(rxn.products)
            icrn.append(Reaction(R, P, rxn.k_fwd, rxn.k_rev))
        pcrn = icrn

    if rates:
        for rxn in pcrn:
            if rxn.k_rev != 0:
                yield '{} <=> {} [kf = {:g}, kr = {:g}]'.format(' + '.join(rxn.reactants), ' + '.join(rxn.products), rxn.k_fwd, rxn.k_rev)
            else:
                yield '{} -> {} [k = {:g}]'.format(' + '.join(rxn.reactants), ' + '.join(rxn.products), rxn.k_fwd)
    else:
        for rxn in pcrn:
            if rxn.k_rev != 0:
                yield '{} <=> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]))
            else:
                yield '{} -> {}'.format(' + '.join(rxn[0]), ' + '.join(rxn[1]))

def split_reversible_reactions(crn):
    """
    Replace every reversible reaction with the two corresponding irreversible
    reactions.
    """
    new = []
    for [r, p, kf, kr] in crn:
        new.append(Reaction(r, p, kf, 0))
        if kr != 0:
            new.append(Reaction(p, r, kr, 0))
    return new

def combine_reversible_reactions(crn):
    """Condense two irreversible reactions into the corresponding reversible reactions. """
    new_crn = []
    removed = []
    for rxn in crn:
        if rxn in removed:
            continue
        [r, p, kf, kr] = rxn
        assert isinstance(r, list) and isinstance(p, list)

        for rxn2 in crn:
            [r2, p2, kf2, kr2] = rxn2
            if sorted(r) == sorted(p2) and sorted(p) == sorted(r2):
                if kr or kr2:
                    raise ValueError('Reaction specified twice!')
                elif sorted(r) == sorted(p):
                    removed.append(rxn2)
                else:
                    removed.append(rxn2)
                    kr = kf2
                break
        new_crn.append(Reaction(r, p, kf, kr))
    return new_crn

def removeSpecies(crn, fuel):
    log.debug('Removing species: {}'.format(fuel))
    assert isinstance(fuel, list)
    # Remove all fuel species, keep rates untouched.
    crn = [Reaction([s for s in rxn.reactants if s not in fuel], 
                    [s for s in rxn.products if s not in fuel], 
                    rxn.k_fwd, rxn.k_rev) for rxn in crn]
    return crn

def assign_crn_species(crn, fs):
    """ Returns a bunch of properties for a given CRN.

    Waste species are all non-formal species that are only products, but never
    react. A non-waste is a formal species, or a species that is involved as a
    reactant in a reaction that involves a non-waste species.

    Return:
        set(): intermediates
        set(): wastes
        set(): nonwastes
    """
    assert isinstance(fs, set)
    species = set().union(*[set().union(*rxn[:2]) for rxn in crn])
    nonwastes = fs.copy()
    wastes = set()
    while True:
        # Add x to non-waste if in any reaction x is an reactant while there
        # are non-wastes taking part in the reaction. Reset the outer loop if 
        # a new non-waste species is found.
        flag = False
        for x in species:
            if x in nonwastes:
                continue
            for [R, P] in crn:
                if x in R and len(nonwastes & set(R + P)):
                    nonwastes.add(x)
                    flag = True
                    break
        if not flag: 
            break
    # A formal species cannot be considered waste.
    wastes = species - nonwastes
    # A non-formal waste species
    nonwastes -= fs
    # An interemdiate species 
    intermediates = species - fs - wastes
    return intermediates, wastes, nonwastes
