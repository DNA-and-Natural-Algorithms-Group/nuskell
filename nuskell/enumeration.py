# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@caltech.edu)
#
# Preprocessing and interface to peppercorn enumerator
#

from nuskell.objects import TestTube, Complex, Reaction, reset_names

from peppercornenumerator.objects import PepperDomain, PepperComplex

from peppercornenumerator import Enumerator
import peppercornenumerator.reactions as reactions
from peppercornenumerator.condense import condense_resting_states


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
                 rename=None, prefix='e'):

        if not (bool(testtube) != bool(enumerator)):  # NOT XOR
            raise ValueError(
                "Need to specify either TestTube() or Enumerator(), but not both!")

        # Autoname resets Peppercorn() complex names and replace it with TestTube()
        # complex names. Two dictionaries store the mapping of names in TestTube()
        # and Enumerator() objects.
        self._rename = rename
        self._enum_prefix = prefix
        self._enum_append = ''
        self._enumN_to_ttubeO = dict()
        self._ttubeN_to_enumO = dict()

        # If True, the enumerator has been called and the TestTube is up-to-date.
        # only flips if used by the enumerate() function in *this* instance.
        self._processed = False

        if testtube:
            self._testtube = testtube
            self._enumerator = self.testtube_to_enumerator(
                self._testtube)  # NOTE Complexes only!!!
        elif enumerator:
            raise NotImplementedError
            self._enumerator = enumerator
            self._testtube = self.enumerator_to_testtube(self._enumerator)

        set_peppercorn_args(self._enumerator, pargs)
        self._enumerator.interruptible = TestTubePeppercornIO.interruptible

    @property
    def testtube(self):
        if not self._processed and self._enumerator._resting_states is not None:
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

        # Translate to peppercorn domains
        domains = {}
        for n, d in testtube.domains.items():
            domains[n] = PepperDomain(name=n, length=d.length)

        # Translate to peppercorn complexes
        complexes = {}
        for cplx in testtube.complexes:
            pepperseq = []
            for d in cplx.sequence:
                if d == '+':
                    pepperseq.append('+')
                else:
                    pepperseq.append(domains[d.name])

            complex = PepperComplex(
                pepperseq, cplx.structure[:], name=cplx.name)
            complexes[cplx.name] = complex

            self._enumN_to_ttubeO[cplx.name] = cplx
            self._ttubeN_to_enumO[cplx.name] = complex

        complexes = complexes.values()

        return Enumerator(complexes)  # domains, strands, complexes)

    def _ecplx_rename(self, x):
        """ A function to rename enumerator species names to a format
        which is compatible with nuskell.objects
        """
        assert x[0].isdigit()
        assert x[-1].isdigit()
        if self._rename is None:
            return self._enum_prefix + x + self._enum_append
        else:
            n = self._enum_prefix + str(self._rename) + self._enum_append
            self._rename += 1
            return n

    def _get_reaction_networkx(self):
        """ Merge the reaction network into self.testtube object. """
        pass

    def enumerator_to_testtube(self, enumerator):
        """Merge Enumerator() complexes into the TestTube() object.

        Note that the enumerator cannot generate new domains or new strands, it
        only finds new complexes. This function adds all complexes that are
        present in resting states in the system. Transient states are ignored.
        """
        reset_names()
        condensed = TestTubePeppercornIO.condensed
        ttcomplexes = dict()

        domains = enumerator.domains
        # The enumerator cannot find new domains, so if there is a testtube object
        # already, then we might as well save some effort and use the testtube
        # domains directly
        if self._testtube:
            dNames = map(lambda x: x.name, domains)
            tDomains = self._testtube.domains
            assert all(map(lambda d: d in tDomains, dNames))
            domains = tDomains
        else:
            # translate domains to nuskell.Domain format, easy...
            raise NotImplementedError

        if condensed:
            enum_complexes = enumerator.resting_complexes
        else:
            # resting_complexes + transient_complexes
            enum_complexes = enumerator.complexes

        for cx in enum_complexes:
            if cx.name in self._enumN_to_ttubeO:
                ttcplx = self._enumN_to_ttubeO[cx.name]
                ttcomplexes[cx.name] = (ttcplx, None, None)
                continue
            domseq = []
            for d in cx.sequence:
                if d == '+':
                    domseq.append('+')
                else:
                    assert (d.name in domains)
                    domseq.append(domains[d.name])
            domstr = cx.structure
            cplxname = cx.name
            if cplxname in ttcomplexes:
                raise ValueError("Complex found in muliple resting states?")
            ttcplx = Complex(sequence=domseq, structure=domstr, name=cplxname)

            ttcomplexes[ttcplx.name] = (ttcplx, None, None)
            self._enumN_to_ttubeO[cx.name] = ttcplx
            self._ttubeN_to_enumO[ttcplx.name] = cx

        crn = dict()
        if condensed:
            condensed = condense_resting_states(enumerator,
                                                compute_rates=True, k_fast=enumerator.k_fast)
            reactions = condensed['reactions']
            for r in reactions:
                # NOTE: takes only the *first complex* of resting states as
                # representative
                react = []
                for rs in r.reactants:
                    react.append(self._enumN_to_ttubeO[rs.complexes[0].name])
                prod = []
                for rs in r.products:
                    prod.append(self._enumN_to_ttubeO[rs.complexes[0].name])
                rxn = Reaction(react, prod, rate=r.rate, rtype=r.rtype)
                crn[rxn.name] = rxn
        else:
            for r in enumerator.reactions:
                react = []
                for re in r.reactants:
                    react.append(self._enumN_to_ttubeO[re.name])
                prod = []
                for pr in r.products:
                    prod.append(self._enumN_to_ttubeO[pr.name])
                rxn = Reaction(react, prod, rate=r.rate, rtype=r.rtype)
                crn[rxn.name] = rxn

        return TestTube(complexes=ttcomplexes, reactions=crn)


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
        enum.max_complex_size = 100

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
    else:
        enum.release_cutoff = None

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
