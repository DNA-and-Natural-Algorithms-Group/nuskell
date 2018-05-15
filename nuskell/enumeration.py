# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@caltech.edu)
#
# Preprocessing and interface to peppercorn enumerator
#

from nuskell.objects import TestTube, NuskellDomain, NuskellComplex, NuskellReaction
from nuskell.objects import clear_memory, DSDDuplicationError
from peppercornenumerator.objects import PepperDomain, PepperComplex

from peppercornenumerator import Enumerator
import peppercornenumerator.reactions as reactions

from peppercornenumerator.condense import PepperCondensation


class TestTubePeppercornIO(object):
    """ A Wrapperclass to communicate between ``Nuskell'' and ``Peppercorn''.

    This class reads a ``Nuskell'' TestTube() object, translates it into a
    ``Peppercorn'' Enumerator() object, enumerates, and updates both the
    Enumerator() and the TestTube() object with the new reactions.

    Alternatively, it may simply read a ``Peppercorn'' Enumerator() object and
    return a ``Nuskell'' TestTube() object. For example, one might start with an
    enumerated DSD network, transfer it to a TestTube() object for sequence design
    and then load it back into the Enumerator().

    Both objects refer to the same complexes and domains.

    Args:
      rename(int, optional): Name enumerated species automatically starting with
        this number. Defaults to None (use peppercorn enumeration number).

    """

    condensed = True
    interruptible = False

    def __init__(self, testtube=None, enumerator=None, pargs=None, nargs=None,
                 rename=None, prefix='e', init_memory=None):

        if not (bool(testtube) != bool(enumerator)):  # NOT XOR
            raise ValueError(
                "Need to specify either TestTube() or Enumerator(), but not both!")

        # Autoname resets Peppercorn() complex names and replace it with TestTube()
        # complex names. Two dictionaries store the mapping of names in TestTube()
        # and Enumerator() objects.
        PepperComplex.PREFIX = prefix

        # If True, the enumerator has been called and the TestTube is up-to-date.
        # only flips if used by the enumerate() function in *this* instance.
        self._processed = False

        # initialize the memory with complex names, just in case they get enumerated.
        self._init_memory = init_memory

        if testtube:
            self._testtube = testtube
            self._enumerator = self.testtube_to_enumerator(self._testtube)
        elif enumerator:
            raise NotImplementedError
            self._enumerator = enumerator
            self._testtube = self.enumerator_to_testtube(self._enumerator)

        set_peppercorn_args(self._enumerator, pargs)
        self._enumerator.interruptible = TestTubePeppercornIO.interruptible

    @property
    def testtube(self):
        if not self._processed and self._enumerator._resting_sets is not None:
            raise ValueError(
                'Enumerator has to be called within the IO object.')
        return self._testtube

    @property
    def enumerator(self):
        return self._enumerator

    def enumerate(self):
        # You may call self._enumerator.enumerate() directly, but then the testtube
        # object will not be updated afterwards. So this wrapper makes sure
        # everythings done properly, and signals that using
        # self._processed=True
        self._enumerator.enumerate()
        self._testtube += self.enumerator_to_testtube(self._enumerator)
        self._processed = True

    def testtube_to_enumerator(self, testtube):  # Does not add reactions!
        """Initialize the peppercorn enumerator object.

        Args:
          testtube <nuskell.objects.TestTube()>: complexes, strands and domains
          tN_to_eO <dict()>: Mapping from testtube names to enumerator Objects()

        Returns:
          Enumerator <peppercorn.enumerator.Enumerator()>
        """

        # Need to clear memory, otherwise DL_Domains.MEMORY will be full with
        # NuskellDomains.
        clear_memory()

        # Translate to peppercorn domains
        domains = {}
        for d in testtube.domains:
            domains[d.name] = PepperDomain(name=d.name, length=d.length)


        # Translate to peppercorn complexes
        complexes = {}
        for cplx in testtube.complexes:
            pepperseq = []
            for d in cplx.sequence:
                if d == '+':
                    pepperseq.append('+')
                else:
                    pepperseq.append(domains[d.name])

            complexes[cplx.name] = PepperComplex(pepperseq, cplx.structure[:], name=cplx.name)

        if self._init_memory:
            for cplx in self._init_memory:
                edible = True
                pepperseq = []
                for d in cplx.sequence:
                    if d == '+':
                        pepperseq.append('+')
                    else:
                        if d.name in domains:
                            pepperseq.append(domains[d.name])
                        else :
                            edible = False
                            break
                if edible:
                    try :
                        tmp = PepperComplex(pepperseq, cplx.structure[:], cplx.name)
                    except DSDDuplicationError:
                        pass

        return Enumerator(complexes.values())  # domains, strands, complexes)

    def enumerator_to_testtube(self, enumerator):
        """Merge Enumerator() complexes into the TestTube() object.

        Note that the enumerator cannot generate new domains or new strands, it
        only finds new complexes. This function adds all complexes that are
        present in resting sets in the system. Transient states are ignored.
        """
        condensed = TestTubePeppercornIO.condensed

        if condensed:
            enumCG = PepperCondensation(enumerator)
            enumCG.condense()
            reactions = enumCG.condensed_reactions
        else:
            reactions = enumerator.reactions

        clear_memory()

        # Translate to peppercorn domains
        domains = {}
        for d in enumerator.domains:
            domains[d.name] = NuskellDomain(name=d.name, length=d.length)

        if condensed:
            enum_complexes = enumCG.resting_set_representatives
        else:
            # resting_complexes + transient_complexes
            enum_complexes = enumerator.complexes

        complexes = {}
        for cplx in enum_complexes:
            nuskellseq = []
            for d in cplx.sequence:
                if d == '+':
                    nuskellseq.append('+')
                else:
                    nuskellseq.append(domains[d.name])

            complexes[cplx.name] = NuskellComplex(nuskellseq, cplx.structure[:], name=cplx.name)

        rxns = []
        for r in reactions:
            react = []
            for rs in r.reactants:
                react.append(NuskellComplex.MEMORY[rs.canonical_form])
            prod = []
            for rs in r.products:
                prod.append(NuskellComplex.MEMORY[rs.canonical_form])
            rxns.append(NuskellReaction(react, prod, rate=r.rate, rtype=r.rtype))

        return TestTube(complexes=complexes.values(), reactions=rxns)


def set_peppercorn_args(enum, args):
    """Transfer options to self._enumerator object.

    Do NOT change default values here. These are supposed to be the defaults of
    peppercorn!  Defaults for nuskell or any other script using this library are
    set with the argparse object of your script, e.g. nuskell: scripts/nuskell.

    """

    if hasattr(args, 'verbose'):
        import logging
        logger = logging.getLogger()
        if args.verbose == 1:
            logger.setLevel(logging.INFO)
        elif args.verbose == 2:
            logger.setLevel(logging.DEBUG)
        elif args.verbose >= 3:
            logger.setLevel(logging.NOTSET)

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
        enum.remote_migration = not args.reject_remote
    else:
        enum.remote_migration = True

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

    if hasattr(args, 'no_max_helix'):
        enum.max_helix_migration = not args.no_max_helix
    else:
        enum.max_helix_migration = True

    if hasattr(args, 'LEGACY_UNZIP'):
        raise DeprecationWarning

    if hasattr(args, 'MAX_COMPLEX_SIZE'):
        raise DeprecationWarning

    if hasattr(args, 'k_slow'):
        enum.k_slow = args.k_slow
    else:
        enum.k_slow = 0.0

    if hasattr(args, 'k_fast'):
        enum.k_fast = args.k_fast
    else:
        enum.k_slow = 0.0

    return
