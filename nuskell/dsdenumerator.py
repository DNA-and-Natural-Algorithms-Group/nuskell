#
# Copyright (c) 2010-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
#  nuskell/dsdenumerator.py
#  NuskellCompilerProject
#
import logging
log = logging.getLogger(__name__)

import gc
from itertools import chain
from peppercornenumerator import enumerate_pil
from peppercornenumerator.enumerator import UNI_REACTIONS
from peppercornenumerator.reactions import branch_3way, branch_4way

from .ioutils import write_pil, load_pil, get_domains
from .objects import NuskellComplex, NuskellMacrostate, NuskellReaction, SingletonError
from .objects import show_memory

class DSDenumerationError(Exception):
    pass

def enumerate_modules(modules, interpretation, solution, reactions, args, prefix = 'm'):
    """ Enumerate all modules, but replaces wildcard species with other signal species.
    """
    seen = set()
    mcomplexes, mreactions = [], []
    for e, module in enumerate(modules, 1):
        # first, replace history complexes with their interpretation!
        for cplx in list(module.values()):
            for k, v in interpretation.items():
                # k, v = A_1_, {A:1}
                if (cplx.name in v) and k != cplx.name:
                    newc = solution[k] 
                    module[k] = newc
                    if cplx.name in module:
                        del module[cplx.name]
                    del newc
            del cplx

        mc, mr = enumerate_solution(module, args, named = solution, prefix = prefix)

        # after enumeration, make sure there were no new 'm' species found.
        for mcplx in list(mc.values()):
            for scplx in solution.values():
                if scplx == mcplx:
                    if scplx.name != mcplx.name:
                        print(f'{scplx.name}, {mcplx.name}, {prefix=}')
                        print(f'{scplx}, {mcplx}')
                        print(f'{e}, {list(map(str, module))=}')
                        print(f'{interpretation=}')
                        print(f'{list(map(str, solution))=}')
                    assert scplx.name == mcplx.name
                    del scplx
                    break
            else:
                raise DSDenumerationError(f'Module complex {mcplx} not found in overall solution!')
            del mcplx
        log.debug(f'Module {e}:\n' + write_pil(mc, mr))
        mcomplexes.append(mc)
        mreactions.append(mr)
        seen |= set(mr)
        del module

    assert set(reactions).issuperset(seen)
    mr = set(reactions) - seen
    if len(mr):
        mc = {x.name: x for rxn in mr for x in chain(rxn.reactants, rxn.products)}
        log.debug(f'Module (crosstalk):\n' + write_pil(mc, mr))
        mcomplexes.append(mc)
        mreactions.append(mr)
    seen.clear()
    return mcomplexes, mreactions

def enumerate_solution(complexes, args, named = None, molarity = 'nM', prefix = 'i'):
    """
    """
    assert all(isinstance(x, NuskellComplex) for x in complexes.values())
    tmp_pil = write_pil(complexes, None, fh = None, molarity = molarity)

    # We want to pass also the named complexes ...
    if named is not None:
        tmp_pil += "\n# Named complexes ...\n"
        domains = get_domains(complexes.values())
        for cx in named.values():
            if cx.name in complexes:
                del cx
                continue
            if not all(d in domains for d in cx.domains):
                continue
            tmp_pil += "{:s} = {:s} @c 0 nM\n".format(cx.name, cx.kernel_string)
            del cx
        del domains
    log.debug(tmp_pil)

    kwargs = get_peppercorn_args(args)
    enum_obj, enum_pil = enumerate_pil(tmp_pil,
                                       detailed = args.enum_detailed,
                                       condensed = not args.enum_detailed,
                                       complex_prefix = prefix,
                                       enumconc = molarity,
                                       **kwargs)
    del enum_obj
    # Now we get the new NuskellComplex objects.  If you enumerate multiple
    # times, e.g. because you enumerate some modules separately, then you want
    # to make sure that the same complexes have the same name.
    reactions = set()

    cxs, rms, det, con = load_pil(enum_pil)
    if args.enum_detailed:
        for name, cx in cxs.items(): # The new complexes.
            try:
                obj = NuskellComplex(None, None, name = name)
            except SingletonError as err:
                obj = NuskellComplex(list(cx.sequence), list(cx.structure), name = name)
                del err
            complexes[obj.name] = obj
            del obj, cx
        for drxn in det:
            reactants = [complexes[s.name] for s in drxn.reactants]
            products = [complexes[s.name] for s in drxn.products]
            try:
                obj = NuskellReaction(reactants, products, drxn.rtype)
            except SingletonError as err:
                obj = err.existing
                del err
            obj.rate_constant = (drxn.rate_constant)
            reactions.add(obj)
            del drxn, obj, reactants, products
    else:
        for name, rm in rms.items():
            cx = rm.representative
            assert name == cx.name
            try:
                obj = NuskellComplex(None, None, name = name)
            except SingletonError as err:
                obj = NuskellComplex(list(cx.sequence), list(cx.structure), name = name)
                del err
            complexes[obj.name] = obj
            del obj, rm, cx
 
        for crxn in con:
            # We extract take objects with the correct names.
            reactants = [complexes[s.name] for s in crxn.reactants]
            products = [complexes[s.name] for s in crxn.products]
            try:
                obj = NuskellReaction(reactants, products, crxn.rtype)
            except SingletonError as err:
                obj = err.existing
                del err
            obj.rate_constant = crxn.rate_constant
            reactions.add(obj)
            del crxn, obj, reactants, products
    cxs.clear()
    rms.clear()
    det.clear()
    con.clear()
    return complexes, reactions

def interpret_species(complexes, reactions, fspecies, prune = True):
    """Get an interpretation dictionary.
    
    If a :obj:`NuskellComplex()` sequence contains a wildcard, then this
    function will find all matching complexes with history domains, rename
    them, and return them in form of a *partial interpretation* dictionary,
    mapping implementation species to (multisets of) formal species.  Complexes
    may have at most *one wildcard domain*, which corresponds to exactly *one
    unpaired long history domain*.

    Args:
        complexes (dict[name] = obj): A dictionary of complex names mapping to the object.
        rections (list[obj]): A list of reaction objects.
        fspecies (list[str], optional): A list of complex names that are potential
            regular-expression complexes.
        prune (bool, optional): Remove all formal species with wildcard domains
            from the network, for which there exists an logically equivalent
            species without the wildcard. Defaults to True.
          
    Example:
      - It is possible to specify sthg like:
        A = "? a b c" | ". . . ."
        B = "a ? b + c* a*" | "( . . + . )"

      - It is not possible to specify sthg like:
        A = "? a b ?" | "( . . )"
        A = "* a b c" | "* . . ."
        A = "? a ? *" | "( . ) *"
        A = "? a ? x + z* x* f* " | "? ? ? ( + . ) ."
        A = "* a * t" | "* . * ."

    Returns:
        dict[impl.name] = Counter([fs.name]): Interpretation dictionary.
        dict[name] = obj: complexes (after pruning)
        list[obj]: reactions (after pruning)
    """
    # Make sure that complexes and reactions point to the same objects
    assert all(id(x) in map(id, complexes.values()) for rxn in reactions \
                                                    for x in chain(rxn.reactants, rxn.products))

    def patternMatch(x, y, ignore='?'):
        """Matches two complexes if they are the same, ignoring history domains.

        Note: The strand order of the second complex changes to the strand order of
        the first complex, if there is a rotation under which both complexes are
        patternMatched.

        Args:
          x (NuskellComplex()) : A nuskell :obj:`NuskellComplex()` object.
          y (NuskellComplex()) : A nuskell :obj:`NuskellComplex()` object.

        Returns: True/False
        """
        if len(list(x.sequence)) != len(list(y.sequence)):
            return False

        def pM_check(pMx, pMy):
            """Recursively parse the current sequences and structures.

            Args:
              pMx [seqX,strX]: A list of two lists (sequence, structrure)
              pMy [seqY,strY]: A list of two lists (sequence, structrure)

            Returns: True/False
            """
            if len(pMx[0]) == 0:
                return True

            if (pMx[0][0] != ignore and pMy[0][0] != ignore) and \
                    (pMx[0][0] != pMy[0][0] or pMx[1][0] != pMy[1][0]):
                return False
            return pM_check([pMx[0][1:], pMx[1][1:]],
                            [pMy[0][1:], pMy[1][1:]])

        pMx = [list(map(str, x.sequence)), list(map(str, x.structure))]
        pMy = [list(map(str, y.sequence)), list(map(str, y.structure))]
        if pM_check(pMx, pMy):
            return True

    def get_matching_complexes(regex, hist):
        """ Find all matching complexes. """
        matching = []
        for cplx in complexes.values():
            if regex.name == cplx.name:
                continue
            elif patternMatch(regex, cplx, ignore = hist):
                matching.append(cplx)
        return matching

    # Find and rename signal species with history domains.
    interpretation = dict()
    need_to_prune = False
    for fs in fspecies:
        if fs not in complexes:
            log.warning(f'No complex found with name of formal species: {fs}')
            continue
        cplx = complexes[fs] # Get corresponding NuskellComplex object.
        hist = [d for d in map(str, cplx.sequence) if d[0] == 'h']
        if hist:
            if len(hist) > 1:
                raise DSDenumerationError('No support for multiple history domains')
            [hist] = hist
            matches = get_matching_complexes(cplx, hist)
            if matches:
                need_to_prune = True
                for e, m in enumerate(matches, 1):
                    newname = fs + '_' + str(e) + '_'
                    log.debug(f'Changing intermediate name {m.name} to signal name {newname}.')
                    del complexes[m.name] # Remove intermediate name from the dictionary.
                    m.name = newname
                    interpretation[m.name] = [fs]
                    complexes[m.name] = m
                del complexes[fs] # Remove wildcard complex from the dictionary.
            else:
                # NOTE: We cannot simply remove the domain, because we would need to
                # enumerate the network again and remove the domain everywhere! So
                # unless we enumerate up-front with a history-pruned species, this
                # gets us into trouble.
                interpretation[cplx.name] = [fs]
        else:
            interpretation[cplx.name] = [fs]

    for rxn in reactions:
        if rxn.name != rxn.auto_name:
            log.debug(f'Updating reaction name from {rxn.name} to {rxn.auto_name}.')
            rxn.name = rxn.auto_name

    log.debug('Interpretation: \n' + '\n'.join(
              [f'{k}: {v}' for k, v in interpretation.items()]))
    log.debug('New complexes: \n' + '\n'.join(
              [f'{k}: {v}' for k, v in complexes.items()]))
    log.debug('Reactions: \n' + '\n'.join([f'{rxn}' for rxn in reactions]))

    if prune and need_to_prune:
        log.debug('Pruning the network.')
        # Get rid of all reactions with history wildcards. Start with a set
        # of produce molecules and see what species emerge from reactions
        # consuming these molecules.
        # Alternative: enumerate again using history-replaced species.
        fuels = [x for x in complexes.keys() if x[0] == 'f']
        [prev, total] = [set(), set(interpretation.keys()) | set(fuels)]
        log.debug('Prev {}, Total {}'.format(prev, total))
        while prev != total:
            prev = set(total) # force a copy?
            for rxn in reactions:
                r = set([x.name for x in rxn.reactants])
                p = set([x.name for x in rxn.products])
                log.debug(f'R {r}, P {p}')
                if r.intersection(total) == r: 
                    total |= p
                #log.debug(f'Total = {total}')
        assert set(map(str, complexes.values())).issuperset(total)
        # Now filter all reactions that are possible from the pruned state space ...
        new_reactions = [r for r in reactions if set(x.name for x in r.reactants).issubset(total)]
        # and remove all the left-over complexes from the graph.
        new_complexes = {k: v for k, v in complexes.items() if k in total}
        reactions.clear() # A horror, but it needs to be done ...
        complexes.clear() # A horror, but it needs to be done ...
        reactions = new_reactions
        complexes = new_complexes

    return interpretation, complexes, reactions

def get_peppercorn_args(args):
    """Transfer options to self._enumerator object.
    Do NOT change default values here. These are supposed to be the defaults of
    peppercorn!  Defaults for nuskell or any other script using this library are
    set with the argparse object of your script, e.g. nuskell: scripts/nuskell.
    """
    kwargs = dict()

    kwargs['max_complex_size'] = args.max_complex_size
    kwargs['max_complex_count'] = args.max_complex_count
    kwargs['max_reaction_count'] = args.max_reaction_count
    kwargs['reject_remote'] = args.reject_remote
    kwargs['max_helix'] = not args.no_max_helix
    if args.ignore_branch_3way:
        if branch_3way in UNI_REACTIONS[1]:
            UNI_REACTIONS[1].remove(branch_3way)
        log.info('No 3-way branch migration.')
    if args.ignore_branch_4way:
        if branch_4way in UNI_REACTIONS[1]:
            UNI_REACTIONS[1].remove(branch_4way)
        log.info('No 4-way branch migration.')

    kwargs['release_cutoff_1_1'] = args.release_cutoff_1_1
    kwargs['release_cutoff_1_2'] = args.release_cutoff_1_2
    if args.release_cutoff:
        kwargs['release_cutoff'] = args.release_cutoff
    kwargs['k_slow'] = args.k_slow
    kwargs['k_fast'] = args.k_fast

    return kwargs
