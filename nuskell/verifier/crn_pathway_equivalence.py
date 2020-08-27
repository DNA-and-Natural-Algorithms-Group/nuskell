#
#  nuskell/verifier/crn_pathway_equivalence.py
#  NuskellCompilerProject
#
# Copyright (c) 2009-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
import logging
log = logging.getLogger(__name__)

from collections import Counter

# TODO: It would be nicer to have clean_crn as part of crnutils ...
from nuskell.crnutils import genCRN, assign_crn_species
from nuskell.verifier.basis_finder import find_basis, NoFormalBasisError, clean_crn
from nuskell.verifier.crn_bisimulation_equivalence import test as get_crn_bisimulation

def test(c1, c2, inter, compositional = False, integrated = False):
    """ Test two CRNs for pathway equivalence.

    Args:
        c1 (list, list): Tuple of formal CRN and formal species.
        c2 (list, list): Tuple of implementation CRN and signal species.
        inter (dict): An interpretation of fs2 in terms of fs1.
        compositional (bool, optional): Use compositional hybrid notion.
            Defaults to False.
        integrated (bool, optional): Use integrated hybrid notion.
            Defaults to False.

    Returns:
        True: if formal basis of crn1 and formal basis of crn2 are equivalent.
        False: otherwise.
    """
    DEVELMODE = True
    # Use the *interpretation upfront* mode of integrated hybrid and the new
    # compositional hybrid implementation.
    if DEVELMODE:
        (integrated, integrated2) = (False, integrated) if integrated else (False, False)
        (compositional, compositional2) = (False, compositional) if compositional else (False, False)

    # Preprocess arguments, just to make sure ...
    crn1, fs1 = clean_crn([rxn[0:2] for rxn in c1[0]]), set(c1[1])
    if integrated2: 
        log.warning('Using integrated "interpretation upfront" mode.')
        crn2, fs2 = clean_crn([rxn[0:2] for rxn in c2[0]], inter = inter), set(c1[1])
    else:
        crn2, fs2 = clean_crn([rxn[0:2] for rxn in c2[0]]), set(c2[1])
        # Note: if inter provides interpretations for non-signal species, that's ok here.
        # They will be used as formal species when finding the formal basis.
        assert all(x in inter for x in fs2)

    assert isinstance(inter, dict)
    # Interpret waste species as nothing.
    intermediates, wastes, reactive_waste = assign_crn_species(crn2, fs2)
    if len(reactive_waste):
        log.warning(f'Reactive waste species detected: {reactive_waste}')
    if len(wastes):
        log.warning(f'{len(wastes)} waste species are treated as formal: ({", ".join(wastes)})')
        for x in wastes: 
            inter[x] = []

    log.debug("Formal CRN:")
    [log.debug('    {}'.format(r)) for r in genCRN(crn1, rates = False)]
    log.debug("")
    log.debug(f"Implementation CRN with {len(set().union(*[set().union(*rxn) for rxn in crn2]))} species and {len(crn2)} reactions.")
    [log.debug(f'    {r}') for r in genCRN(crn2, rates = False)]
    log.debug("")
    if integrated:
        log.debug(f"Implementation CRN after partial interpretation:")
        [log.debug(f'    {r}') for r in genCRN(crn2, interpretation = inter, rates = False)]
        log.debug("")

    try:
        log.debug(f'Formal species to find basis: {fs2 | wastes}')
        fbasis_raw, fbasis_int = find_basis(crn2, fs2 | wastes, modular = True,
                                            interpretation = inter if integrated else None)
    except NoFormalBasisError as err:
        log.info("Could not find formal basis: {}".format(err))
        return False

    log.info(f"Raw formal basis:")
    [log.info(f'    {r}') for r in genCRN(fbasis_raw, rates = False)]
    log.info(f"Interpreted formal basis:")
    [log.info(f'    {r}') for r in genCRN(fbasis_int, rates = False)]

    if compositional:
        return sorted(crn1) == sorted(clean_crn(fbasis_raw, inter = inter))
    elif compositional2:
        # Currently, compositional and integrated use the same raw basis,
        # because the raw basis after interpretation is the interpreted basis.
        fcrn = [[Counter(part) for part in rxn] for rxn in crn1]
        icrn = [[Counter(part) for part in rxn] for rxn in fbasis_raw]

        interC = dict()
        for k, v in inter.items():
            interC[k] = Counter(v)

        v, i = get_crn_bisimulation(fcrn, icrn, fs1, 
                                    interpretation = interC, 
                                    permissive = 'whole-graph')
        return v
    elif integrated:
        # If you use the standard integrated implementation, then you need
        # to use the interpreted formal base.
        return sorted(crn1) == sorted(clean_crn(fbasis_int))

    # Pure CRN pathway decomposition or *new* integrated implementation.
    return sorted(crn1) == sorted(fbasis_raw)

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

