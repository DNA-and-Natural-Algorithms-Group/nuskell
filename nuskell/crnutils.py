#
#  nuskell/crnutils.py
#  NuskellCompilerProject
#
import logging
log = logging.getLogger(__name__)

from itertools import chain
from dsdobjects.utils import natural_sort
from .dsdcompiler.crn_parser import (Reaction, 
                                     parse_crn_string, 
                                     parse_crn_file,
                                     CRNParseError)
from crnverifier.utils import assign_crn_species

class CRNerror(Exception):
    pass

def genCON(species):
    for sp, data in natural_sort(species.items()):
        if None in data: continue
        yield '{} @{} {}'.format(sp, 'initial' if data[0][0] == 'i' else 'constant', data[1])

def interpret(l, inter):
    """ Replace species with their interpretation. """
    return list(chain(*[inter.get(x, [x]) for x in l]))

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
        pcrn = combine_reversible_reactions(crn)
    else:
        pcrn = split_reversible_reactions(crn)

    if interpretation:
        icrn = []
        for rxn in pcrn:
            R = interpret(rxn.reactants, interpretation)
            P = interpret(rxn.products, interpretation)
            icrn.append(Reaction(R, P, rxn.k_fwd, rxn.k_rev))
        pcrn = removeDuplicates(removeTrivial(icrn))

    for rxn in pcrn:
        R = natural_sort(rxn.reactants)
        P = natural_sort(rxn.products)
        if rxn.k_rev == 0:
            rate = ' [k = {:g}]'.format(rxn.k_fwd) if rates else ''
            yield '{} -> {}{}'.format(' + '.join(R), ' + '.join(P), rate)
        else:
            if natural_sort([R, P]) == [R, P]: 
                rate = ' [kf = {:g}, kr = {:g}]'.format(rxn.k_fwd, rxn.k_rev) if rates else ''
                yield '{} <=> {}{}'.format(' + '.join(R), ' + '.join(P), rate)
            else:
                rate = ' [kf = {:g}, kr = {:g}]'.format(rxn.k_rev, rxn.k_fwd) if rates else ''
                yield '{} <=> {}{}'.format(' + '.join(P), ' + '.join(R), rate)

def split_reversible_reactions(crn):
    """
    Replace every reversible reaction with the two corresponding irreversible
    reactions.
    """
    new = []
    for [r, p, kf, kr] in crn:
        if sorted(r) == sorted(p):
            raise CRNerror('Trivial reaction found!')
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

        if sorted(r) == sorted(p):
            raise CRNerror('Trivial reaction found!')

        for rxn2 in crn:
            [r2, p2, kf2, kr2] = rxn2
            if sorted(r) == sorted(p2) and sorted(p) == sorted(r2):
                if kr or kr2:
                    raise CRNerror('Reaction specified twice!')
                else:
                    removed.append(rxn2)
                    kr = kf2
                break
        new_crn.append(Reaction(r, p, kf, kr))
    return new_crn

def removeSpecies(crn, fuel):
    log.debug('Removing species: {}'.format(fuel))
    assert isinstance(fuel, (list, set))
    # Remove all fuel species, keep rates untouched.
    crn = [Reaction([s for s in rxn.reactants if s not in fuel], 
                    [s for s in rxn.products if s not in fuel], 
                    rxn.k_fwd, rxn.k_rev) for rxn in crn]
    return crn

def removeTrivial(crn):
    log.debug('Removing trivial reactions.')
    # Remove all fuel species, keep rates untouched.
    return [rxn for rxn in crn if sorted(rxn.reactants) != sorted(rxn.products)]

def removeDuplicates(crn):
    log.debug('Removing duplicate reactions.')
    # Remove all fuel species, keep rates untouched.
    seen = []
    new = []
    for rxn in crn:
        if rxn not in seen:
            new.append(rxn)
        seen.append(rxn)
    return new

