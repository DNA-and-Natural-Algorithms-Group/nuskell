from collections import namedtuple
from crnsimulator import parse_crn_string as pcs
from dsdobjects.utils import natural_sort
from crnsimulator.crn_parser import CRNParseError

Reaction = namedtuple('reaction', 'reactants products k_fwd k_rev')

def parse_crn_string(crnstring):
    # A wrapper for the crnsimulator version
    return post_process(pcs(crnstring, process = False), defaultrate = 1, defaultmode = None, defaultconc = None)

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

def genCON(species):
    for sp, data in natural_sort(species.items()):
        if None in data: continue
        yield '{} @{} {}'.format(sp, 'initial' if data[0][0] == 'i' else 'constant', data[1])

def genCRN(crn, reversible = True, rates = True):
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
                else:
                    removed.append(rxn2)
                    kr = kf2
                break
        new_crn.append(Reaction(r, p, kf, kr))
    return new_crn

def removeSpecies(crn, fuel):
    # Remove all fuel species, keep rates untouched.
    crn = [Reaction(filter(lambda s: s not in fuel, rxn.reactants), filter(lambda s: s not in fuel, rxn.products), rxn.k_fwd, rxn.k_rev) for rxn in crn]
    return crn

