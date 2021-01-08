#
#  nuskell/crnutils.py
#  NuskellCompilerProject
#
import logging
log = logging.getLogger(__name__)

from itertools import chain
from natsort import natsorted
from .dsdcompiler.crn_parser import (Reaction, 
                                     parse_crn_string, 
                                     parse_crn_file,
                                     CRNParseError)
from crnverifier.utils import assign_crn_species

class CRNerror(Exception):
    pass

def interpret(l, inter):
    """ Replace species with their interpretation. """
    assert all(isinstance(v, list) for v in inter.values())
    return list(chain(*[inter.get(x, [x]) for x in l]))

def genCON(species):
    for sp, data in natsorted(species.items()):
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
        pcrn = combine_reversible_rxns(crn)
    else:
        pcrn = split_reversible_rxns(crn)

    if interpretation:
        icrn = []
        for rxn in pcrn:
            R = interpret(rxn.reactants, interpretation)
            P = interpret(rxn.products, interpretation)
            icrn.append(Reaction(R, P, rxn.k_fwd, rxn.k_rev))
        pcrn = remove_duplicate_rxns(remove_trivial_rxns(icrn))

    for rxn in pcrn:
        R = natsorted(rxn.reactants)
        P = natsorted(rxn.products)
        if rxn.k_rev == 0:
            rate = ' [k = {:g}]'.format(rxn.k_fwd) if rates else ''
            yield '{} -> {}{}'.format(' + '.join(R), ' + '.join(P), rate)
        else:
            if natsorted([R, P]) == [R, P]: 
                rate = ' [kf = {:g}, kr = {:g}]'.format(rxn.k_fwd, rxn.k_rev) if rates else ''
                yield '{} <=> {}{}'.format(' + '.join(R), ' + '.join(P), rate)
            else:
                rate = ' [kf = {:g}, kr = {:g}]'.format(rxn.k_rev, rxn.k_fwd) if rates else ''
                yield '{} <=> {}{}'.format(' + '.join(P), ' + '.join(R), rate)

def split_reversible_rxns(crn):
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

def combine_reversible_rxns(crn):
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

def remove_species(crn, fuel):
    # Remove all fuel or waste species, keep rates untouched.
    assert isinstance(fuel, (list, set))
    log.debug('Removing species: {}'.format(fuel))
    crn = [Reaction([s for s in rxn.reactants if s not in fuel], 
                    [s for s in rxn.products if s not in fuel], 
                    rxn.k_fwd, rxn.k_rev) for rxn in crn]
    return crn

def remove_trivial_rxns(crn):
    # Remove all trivial reacitons, keep rates untouched.
    log.debug('Removing trivial reactions.')
    return [rxn for rxn in crn if sorted(rxn.reactants) != sorted(rxn.products)]

def remove_duplicate_rxns(crn, reversible = True):
    # Remove all duplicate reactions, keep rates untouched.
    log.debug('Removing duplicate reactions.')
    new, seen = [], set()
    for rxn in split_reversible_rxns(crn):
        R = tuple(sorted(rxn.reactants))
        P = tuple(sorted(rxn.products))
        if (R, P) in seen:
            continue
        new.append(rxn)
        seen.add((R,P))
    return combine_reversible_rxns(new) if reversible else new

def cleanup_rxns(crn, reversible = True):
    # Remove all duplicate reactions, keep rates untouched.
    log.debug('Removing duplicate reactions.')
    new, seen = [], set()
    for rxn in split_reversible_rxns(crn):
        R = tuple(sorted(rxn.reactants))
        P = tuple(sorted(rxn.products))
        if R == P:
            continue
        if (R, P) in seen:
            continue
        new.append(rxn)
        seen.add((R,P))
    return combine_reversible_rxns(new) if reversible else new

