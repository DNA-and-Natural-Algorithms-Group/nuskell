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

from peppercornenumerator import enumerate_pil
from peppercornenumerator.objects import PepperComplex
from dsdobjects import clear_memory, DSDDuplicationError
from dsdobjects.core import DSD_Complex

from .ioutils import write_pil, load_pil
from .objects import NuskellReaction

class DSDenumerationError(Exception):
    pass

def enumerate_modules(modules, interpretation, solution, reactions, args):
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

        mc, mr = enumerate_solution(module, args, prefix = 'm')
        # after enumeration, make sure there were no new 'm' species found.
        for mcplx in list(mc.values()):
            for scplx in solution.values():
                if scplx == mcplx:
                    assert scplx.name == mcplx.name
                    break
            else:
                raise DSDenumerationError(f'Module complex {mcplx} not found in overall solution!')

        log.debug(f'Module {e}:\n' + write_pil(mc, mr))
        mcomplexes.append(mc)
        mreactions.append(mr)
        seen |= set(mr)

    assert set(reactions).issuperset(seen)
    mr = set(reactions) - seen
    if len(mr):
        mc = {x.name: x for rxn in mr for x in rxn.reactants + rxn.products}
        log.debug(f'Module (crosstalk):\n' + write_pil(mc, mr))
        mcomplexes.append(mc)
        mreactions.append(mr)
    return mcomplexes, mreactions

def enumerate_solution(complexes, args, prefix = 'i'):
    """
    # TODO: pass arguments for peppercorn!
    """
    PepperComplex.PREFIX = prefix

    # A wrapper for enumeration? -> dsdenumerate.py
    tmp_pil = write_pil(complexes, None, fh = None,
                        molarity = args.concentration_units)

    # Memory management.
    backupCM = DSD_Complex.MEMORY
    # Since the enumerator uses DSDObjects, we need to clear the memory of
    # known complexes, domains, etc. Otherwise we get duplication errors.
    clear_memory() 

    enum_obj, enum_pil = enumerate_pil(tmp_pil, 
                                       detailed = args.enum_detailed, 
                                       condensed = not args.enum_detailed, 
                                       enumconc = args.concentration_units)
    
    clear_memory() 
    # Now we get the new NuskellComplex objects.  If you enumerate multiple
    # times, e.g. because you enumerate some modules separately, then you 
    # want to make sure that the same complexes have the same name.
    dom, clx, rms, det, con = load_pil(enum_pil)

    complexes = dict()
    reactions = []
    if args.enum_detailed:
        #assert all([x in clx for x in complexes]) # w.o. renaming. 
        for name, obj in clx.items(): # The new complexes.
            if obj.canonical_form in backupCM:
                # That should take care of MEMORY management.
                if obj.name != backupCM[obj.canonical_form].name:
                    obj.name = backupCM[obj.canonical_form].name
            complexes[obj.name] = obj
        reactions = det
    else:
        # assert all([x in rms for x in complexes]) # w.o. renaming
        for name, rms in rms.items():
            obj = rms.canonical_complex
            if obj.canonical_form in backupCM:
                # That should take care of MEMORY management.
                if obj.name != backupCM[obj.canonical_form].name:
                    obj.name = backupCM[obj.canonical_form].name
            complexes[obj.name] = obj
 
        for rxn in con:
            # We extract take objects with the correct names.
            reactants = [DSD_Complex.MEMORY[x.canonical_form] for x in rxn.reactants]
            products = [DSD_Complex.MEMORY[x.canonical_form] for x in rxn.products]
            try:
                new = NuskellReaction(reactants, products, rxn.rtype, rxn.rate)
            except DSDDuplicationError as err:
                # NOTE: well those reactions don't distinguish complexes and
                # macrotates, unfortunately, so let's overwrite it.
                new = NuskellReaction(reactants, products, 
                                      rxn.rtype, rxn.rate, 
                                      memorycheck = False)
                err.existing = new
            reactions.append(new)
    # Now at last make sure that we remain 
    DSD_Complex.MEMORY.update(backupCM)
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
                                                     for x in rxn._reactants + rxn._products)

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
        if len(x.sequence) != len(y.sequence):
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
        """Find all matching complexes. """
        regseq = regex.sequence
        regstr = regex.structure

        matching = []
        for cplx in complexes.values():
            if regex.name == cplx.name:
                continue
            elif patternMatch(regex, cplx, ignore=hist):
                matching.append(cplx)
        return matching

    # Find and rename signal species with history domains.
    need_to_prune = False
    interpretation = dict()
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
                    del complexes[m.name] # Remove intermediate name from the dictionary.
                    m.name = fs + '_' + str(e) + '_'
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
                r = set(map(str, rxn.reactants))
                p = set(map(str, rxn.products))
                #log.debug(f'R {r}, P {p}')
                if r.intersection(total) == r: 
                    total |= p
                #log.debug(f'Total = {total}')
        assert set(map(str, complexes.values())).issuperset(total)
        # Now filter all reactions that are possible from the pruned state space ...
        reactions = [r for r in reactions if set(map(str, r.reactants)).issubset(total)]
        # and remove all the left-over complexes from the graph.
        complexes = {k: v for k, v in complexes.items() if k in total}
    return interpretation, complexes, reactions

def set_peppercorn_args(enum, args):
    """Transfer options to self._enumerator object.
    Do NOT change default values here. These are supposed to be the defaults of
    peppercorn!  Defaults for nuskell or any other script using this library are
    set with the argparse object of your script, e.g. nuskell: scripts/nuskell.
    """

    if hasattr(args, 'max_complex_size'):
        enum.max_complex_size = args.max_complex_size
    else:
        enum.max_complex_size = 20

    if hasattr(args, 'max_complex_count'):
        enum.max_complex_count = args.max_complex_count
    else:
        enum.max_complex_count = 1000

    if hasattr(args, 'max_reaction_count'):
        enum.max_reaction_count = args.max_reaction_count
    else:
        enum.max_reaction_count = 5000

    if hasattr(args, 'reject_remote'):
        enum.reject_remote = args.reject_remote
    else:
        enum.reject_remote = False

    if hasattr(args, 'ignore_branch_3way') and args.ignore_branch_3way:
        if reactions.branch_3way in enum.FAST_REACTIONS:
            enum.FAST_REACTIONS.remove(reactions.branch_3way)

    if hasattr(args, 'ignore_branch_4way') and args.ignore_branch_4way:
        if reactions.branch_4way in enum.FAST_REACTIONS:
            enum.FAST_REACTIONS.remove(reactions.branch_4way)

    if hasattr(args, 'release_cutoff_1_1'):
        enum.release_cutoff_1_1 = args.release_cutoff_1_1
    else:
        enum.release_cutoff_1_1 = 6

    if hasattr(args, 'release_cutoff_1_N'):
        enum.release_cutoff_1_N = args.release_cutoff_1_N
    else:
        enum.release_cutoff_1_N = 6

    if hasattr(args, 'release_cutoff'):
        if args.release_cutoff is not None:
            enum.release_cutoff_1_1 = args.release_cutoff
            enum.release_cutoff_1_N = args.release_cutoff
            enum.release_cutoff = args.release_cutoff
