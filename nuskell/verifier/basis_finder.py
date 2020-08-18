#!/usr/bin/env python
#
#  nuskell/verifier/verifier.py
#  NuskellCompilerProject
#
# Copyright (c) 2009-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (bad-ants-fleet@posteo.eu)
#
import logging
log = logging.getLogger(__name__)

from itertools import chain
from collections import Counter

from nuskell import __version__
from nuskell.crnutils import (parse_crn_file, 
                              parse_crn_string, 
                              split_reversible_reactions,
                              assign_crn_species)

class NoFormalBasisError(Exception):
    pass

class BasisFinderError(Exception):
    pass

def pretty_rxn(rxn):
    return '{} -> {}'.format(' + '.join(sorted(rxn[0])), 
                             ' + '.join(sorted(rxn[1])))

def pretty_crn(crn):
    for rxn in crn:
        yield pretty_rxn(rxn)

def clean_crn(crn, duplicates = True, trivial = True):
    """Takes a crn of Counter objects and removes trivail / duplicate reactions. """
    new = []
    seen = set()
    for [R, P] in crn:
        lR = sorted(R)
        lP = sorted(P)
        tR = tuple(lR)
        tP = tuple(lP)
        if trivial and tR == tP:
            continue
        if duplicates and (tR, tP) in seen:
            continue
        new.append([lR, lP])
        seen.add((tR, tP))
    return new

def interpret(l, inter):
    """ Replace species with their interpretation. """
    return list(chain(*[inter.get(x, [x]) for x in l]))

def formal(s, fs):
    """ Returns all formal species in state s. """
    assert isinstance(s, (list, chain, tuple))
    return [x for x in s if x in fs]

def intermediate(s, fs):
    """ Returns all non-formal species in state s. """
    assert isinstance(s, (list, chain, tuple))
    return [x for x in s if x not in fs]

def is_formal_state(s, fs):
    """ True if the state s is a formal state, False otherwise. """
    return len(intermediate(s, fs)) == 0

def is_subset(a, b):
    """ True if *a* is a subset of *b*, False otherwise. """
    assert isinstance(a, Counter) and isinstance(b, Counter)
    return a & b == a

class Path:
    """ A sequence of reactions.

    Can be initialized using a list of reactions in *list* or *Counter* format.
    Use __add__() to save compuatation time when concatinating two paths.

    TODO: Think about how we can save compuation by expressing path decomposion
    as a combination of path objects.
    """
    def __init__(self, path, fs = None):
        assert isinstance(path, list)
        if len(path) == 0:
            self.lpath = self.cpath = []
        else:
            assert len(path[0]) == 2
            if isinstance(path[0][0], Counter):
                self.cpath = path
                self.lpath = [[list(R.elements()), list(P.elements())] for [R, P] in path]
            else: 
                self.cpath = [[Counter(R), Counter(P)] for [R, P] in path]
                self.lpath = path
        self.fs = fs # Formal species

        # set by: self._get_minits()
        self._minit = None
        self._final = None
        # set by: self.width
        self._width = None
        # set by: self.formal_closure
        self._formc = None
        # set by: self.states()
        self._cstates = None
        # set by: decomposed_final_states
        self._dfinals = None
        # set by: is_linear()
        self._linear = None
        # TODO
        self._dfs = None
        self._rfs = None

    @property
    def pretty_path(self):
        return map(pretty_rxn, self.lpath)

    @property
    def S0(self):
        """ tuple: hashable initial state. """
        return tuple(sorted(self.minimal_initial_state.elements()))

    @property
    def minimal_initial_state(self):
        """ Counter: minimal multiset of species to complete this path. """
        if self._minit is None:
            self._get_minits()
        return self._minit

    @property
    def Sn(self):
        """ tuple: hashable final state. """
        return tuple(sorted(self.final_state.elements()))

    @property
    def final_state(self):
        """ Counter: minimal multiset of species after completing this path. """
        if self._final is None:
            self._get_minits()
        return self._final

    @property
    def is_semiformal(self):
        """ bool: True if path has a formal initial state, False otherwise. """
        return is_formal_state(self.minimal_initial_state.elements(), self.fs)

    @property
    def is_closed(self):
        """ bool: True if path has a formal final state, False otherwise. """
        return is_formal_state(self.final_state.elements(), self.fs)

    @property
    def is_formal(self):
        """ bool: True if path has formal initial and formal final state, False otherwise. """
        # semi-formal & closed = formal
        return is_formal_state(self.minimal_initial_state.elements(), self.fs) \
                and is_formal_state(self.final_state.elements(), self.fs)

    @property
    def is_intermediate(self):
        """ bool: True if path has intermediate initial state, False otherwise. """
        return len(formal(self.minimal_initial_state.elements(), self.fs)) == 0 and \
               len(intermediate(self.minimal_initial_state.elements(), self.fs)) >= 1

    @property
    def is_prime(self):
        """ bool: True if path is undecomposable, False otherwise. """
        return len(self.dfs) == 0

    def is_linear(self, nonwastes):
        """ True if path has monomolecular substructure, False otherwise. 
        Args: nonwastes (list): non-formal reactants. """
        if self._linear is None:
            for state in self.states(counter = False):
                c1 = [x for x in state if x in nonwastes]
                if len(c1) > 1:
                    self._linear = False
                    break
            else:
                self._linear = True
        return self._linear

    @property
    def width(self):
        """ int: the largest state when folloing the path from minimal initial state. """
        if self._width is None:
            if len(self) == 0:
                self._width = 0
            else:
                self._width = max([sum(s.values()) for s in self.states(counter = True)]) 
        return self._width

    @property
    def fc(self):
        """ tuple: hashable formal closure. """
        return tuple(sorted(self.formal_closure.elements()))

    @property
    def formal_closure(self):
        """ Counter: minimal multiset of formal species involved in a path. """
        if self._formc is None:
            self._formc = Counter()
            for Si in self.states(counter = True):
                self._formc |= self.formal(Si)
        return self._formc

    @property
    def dfs(self):
        """ frozenset: a wrapper for decomposed final states. """
        return self.decomposed_final_states()

    @property
    def rfs(self):
        """ frozenset: a wrapper for regular final states. """
        return self.regular_final_states() # if len(self.dfs) == 0 else None
 
    def decomposed_final_states(self, path = None):
        """ Calculate decomposed final states.

        Note: If self._dfinal is None, then this routine always uses the internal path,
            otherwise it will update the internal final states with the given path.
            This is a speedup used by the __add__ function.
        """
        if self.fs is None:
            raise BasisFinderError('DFS: no formal species given.')

        if self._dfinals is None:
            self._dfinals = [(Counter(), Counter())]
            assert path is None
            path = self.cpath
        elif path is None:
            path = []

        def next_state_C(S, rxn):
            """ Counter: applies reaction rxn to state S. """
            if not is_subset(rxn[0], S):
                raise BasisFinderError('Reaction {} cannot occur in state S {}.'.format(S, rxn))
            return S - rxn[0] + rxn[1]

        def next_final_states(t1, t2, R, P):
            try:
                tmp1 = t1 - self.formal(R) # consume formal species, if they are present.
                tmp1 = next_state_C(tmp1, [self.intermediate(R), P]) # skip formal reactants
            except BasisFinderError as err:
                tmp1 = None
            if t1 == t2:
                return [(tmp1, t2), None] if tmp1 is not None else [None, None]
            try:
                tmp2 = t2 - self.formal(R) # consume formal species, if they are present.
                tmp2 = next_state_C(tmp2, [self.intermediate(R), P]) # skip formal reactants
            except BasisFinderError as err:
                tmp2 = None
            if tmp1 is None and tmp2 is None:
                return None, None
            elif tmp2 is None:
                return (tmp1, t2), None
            elif tmp1 is None:
                return (t1, tmp2), None
            else:
                return (tmp1, t2), (t1, tmp2)

        debug = False
        new_finals = []
        for [R, P] in path:
            new_finals = []
            for (T1, T2) in self._dfinals:
                if debug: 
                     if (is_formal_state(T1.elements(), self.fs) and sum(T1.values()) != 0) or \
                        (is_formal_state(T2.elements(), self.fs) and sum(T2.values()) != 0):
                        raise BasisFinderError('DFS: found formal final state during decomposition!')
                n1, n2 = next_final_states(T1, T2, R, P)
                if n1 is not None and n1 not in new_finals[1:]:
                    new_finals.append((n1))
                if n2 is not None and n2 not in new_finals[1:]:
                    new_finals.append((n2))
            if len(new_finals) < 1:
                raise BasisFinderError('DFS: not a semi-formal pathway!')
            self._dfinals = new_finals

        if len(self._dfinals) == 1:
            [(final, _)] = self._dfinals
            assert sum(_.values()) == 0 and final == self.final_state
            return frozenset()
        else:
            DFS = set()
            for e, (T1, T2) in enumerate(self._dfinals):
                if e == 0:
                    assert sum(T2.values()) == 0
                    (final, _) = (T1, T2)
                    assert final == self.final_state
                else:
                    DFS.add(tuple(sorted(T1.elements())))
                    DFS.add(tuple(sorted(T2.elements())))
            return frozenset(DFS)

    def regular_final_states(self):
        """ Calculates the set of final states assuming that the path is regular.

        A regular path must have a *turning point* reaction.

        Say a reaction (R,P)_t transforms state S_{t-1} to S_t. Fi is short for
        formal(Si).  A *turning point reaction* must fulfill several conditions:
              - every Fi with i<t is contained in S0 (= F0)
              - every Fi with i>=t is contained in final state Fn
              - formal(S_{t-1} - R_t) == formal(S_t - P_t) == 0

        A note for the last condition: "If there are formal species at state
        t-1, then those formal species have to be consumed by the turning
        point reaction". Equivalently, "if there are formal species at state t,
        then those species must have been produced by the turning point
        reaction".
       
        As example, a reaction (A -> A + i) does not change the formal species
        between S_{t-1} and S_t, but may be considered turning point reaction
        (if and only if A is in S0 and Sn). A path (A -> i; i -> B) has two
        turning point reactions. However, (A -> A + i; i -> B) must only have
        one turning point reaction, that corresponds to the formal basis S0 =
        {A} and Sn = {A, B}.  Since the states all contain formal species 
        (S0: {A}, S1: {A, i}, S2: {A, B}), only the reaction that guarantees a
        complete turnover of formal species is considered a turning point 
        (A -> A + i).

        Sidenote: Because waste species are treated as formal, they need to be
        produced after the turning point!
        """
        RFS = set() # Regular final states

        # Find last possible turning point j.
        # S0 --rxn1--> S1 --rxn2--> S2 ...
        # every Si < j is subset of S0, but Sj is not!
        states = self.states()
        S0 = self.minimal_initial_state
        for j, Si in enumerate(states):
            if not is_subset(self.formal(Si), S0):
                break

        T = Counter() # Final state
        for i in range(len(self), 0, -1):
            T |= self.formal(states[i])
            if i > j: 
                # Collect all formal species after turning point.
                continue
            # Get potential turning point reaction.
            [R, P] = self.cpath[i-1] # rxn 1 at pos 0
            remainder = self.formal(states[i-1] - R)
            if sum(remainder.values()) == 0:
                t = tuple(sorted(T.elements()))
                RFS.add(t)
        return frozenset(RFS)

    @property
    def signature(self):
        """ tuple: the signature of a path. """
        i = self.S0
        t = self.Sn
        w = self.width
        f = self.fc
        d = self.dfs
        r = self.rfs
        return tuple([i, t, w, f, d, r])

    def formal(self, S):
        """ Counter: multiset of formal species in the state S. """
        if self.fs is None:
            raise BasisFinderError('Path: no formal species given.')
        return Counter({k: v for k,v in S.items() if k in self.fs})

    def intermediate(self, S):
        """ Counter: multiset of non-formal species in the state S. """
        if self.fs is None:
            raise BasisFinderError('Path: no formal species given.')
        return Counter({k: v for k,v in S.items() if k not in self.fs})

    def states(self, counter = True):
        """ A list of states from S0 to Sn. 
        Args: counter (bool): returns states as Counters if True, returns lists otherwise.
        """
        if self._cstates is None:
            self._cstates = []
            self._cstates.append(self.minimal_initial_state)
            for [R, P] in self.cpath:
                self._cstates.append(self._cstates[-1] - R + P)
        if counter:
            # Typically, it is much more readable to work with counters.
            return self._cstates
        else:
            # Sometimes it is a *lot* faster to work with lists.
            return [sorted(s.elements()) for s in self._cstates]

    def _get_minits(self, path = None):
        """ Initialize internal S0 and Sn for this path object. """
        if path is None:
            path = self.cpath
        if self._minit is None:
            self._minit = Counter() # initial state
        if self._final is None:
            self._final = Counter() # current state
        for [R, P] in path:
            self._minit += (R - self._final)
            self._final = (self._final - R) + P

    def __add__(self, other):
        """Path: the concatenation of two Path objects. """
        assert type(self) == type(other)
        combo = Path(self.lpath + other.lpath)
        combo.fs = self.fs # force copy?

        # set by: self._get_minits()
        if self._minit:
            combo._minit = self._minit + Counter()
            combo._final = self._final + Counter()
            combo._get_minits(path = other.cpath)
        # set by: decomposed_final_states
        if combo.is_semiformal and self._dfinals is not None:
            combo._dfinals = self._dfinals
            combo.decomposed_final_states(path = other.cpath)
        return combo

    def __len__(self):
        return len(self.lpath)

    def __eq__(self, other):
        return self.signature == other.signature

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.signature)

    def __str__(self):
        return '{}'.format('; '.join(self.pretty_path))

    def __repr__(self):
        return "Path(l={}, {})".format(len(self), self.signature)

def tidy(queue, crn, fs, TC = None, bound = None):
    """BFS to test if a pathway p has a strong closing pathway q.

    Typically, queue contains only one state: all intermediate species of the
    final state T of a pathay p. The function attempts to find a closing
    pathway q, which produces a formal state (without requiring formal species
    as reactants).

    Args:
        queue (list): A list of states (states must contain only intermediate species).
        crn (list): The CRN.
        fs (set): A set of formal species.
        TC (dict, optional): A Tidy-Check dictionary that cotains known tidy states.
        bound (int, optional): A width bound to limit the search for a pathway q.

    Returns:
        True: if a strong closing pathway exists
        False: if no strong closing pathway exists
        queue: a list of states that exceed the current bound
    """
 
    def next_state(s, rxn):
        # Applies reaction rxn to state s.
        s = list(s)
        for x in rxn[0]:
            s.remove(x)
        s += rxn[1]
        return s

    # Keep this. queue should not be a single state but a list of states.
    assert all(isinstance(x, (list, tuple)) for x in queue)
    assert all(intermediate(x, fs) == list(x) for x in queue)

    mem = queue[:]
    while len(queue) > 0:
        queue.sort(key = len)
        iS = queue[0] 
        if len(iS) == 0:
            return True
        elif TC and TC.get(iS, False) is True: 
            return True
        elif bound and len(iS) > bound:
            return queue
        queue = queue[1:]
        for e, rxn in enumerate(crn):
            try:
                nS = tuple(sorted(intermediate(next_state(iS, rxn), fs)))
                if nS not in mem: 
                    mem.append(nS)
                    queue.append(nS)
            except ValueError as err:
                # iS did not contain the required reactants.
                continue
    return False

def crn_properties(crn, fs):
    """ Get different CRN properties.

    Returns:
        b (int): the branching factor
        iR: max over |intermediate(R)| of (R,P) in crn
        fRiR: set of (|formal(R)|, |intermediate(R)|)
        nw: nonwaste species for linear check
    """
    i, w, nw = assign_crn_species(crn, fs)
    assert i == nw
    linear = True
    for [r, p] in crn:
        r1 = [x for x in r if x in i and x not in w]
        p1 = [x for x in p if x in i and x not in w]
        if len(r1) > 1 or len(p1) > 1: 
            linear = False
            break
    log.debug('Monomolecular substructure: {}'.format(linear))
    b  = max(max(len(r),len(p)) for [r, p] in crn)
    iR = max(len(intermediate(r, fs)) for [r, p] in crn) 
    fRiR = set((len(formal(r, fs)), len(intermediate(r, fs))) for [r, p] in crn)
    log.debug('Branching Factors: b = {}, iR = {}, fRiR = {}'.format(b, iR, fRiR))
    return b, iR, fRiR, nw if linear else None

def get_formal_basis(crn, fs, inter = None): # former "enumerate_basis"
    """ Find the formal basis of a CRN.

    Args:
        crn (List[[R], [P]]): A list of irreversible reactions. First a list of
            reactants, second a list of products.
        fs (set(str)): A set of formal species.
        interpretation (dict, optional): An interpretation dictionary for the 
            "integrated hybrid" mode. Defaults to None.

    Returns:
        basis_raw, basis_int: The formal basis found without and with an 
            interpretation dictionary. 
    """
    TidyCheck = dict() # A cache to remember which final states are Tidy. 
    b, b_r, fRiR, nonw = crn_properties(crn, fs)

    # Uses some constants: w_max, i_max, nonw, fs, but not crn (in case we want it random)
    def enumerate_elementary_basis(p, crn):
        """Enumerates all signatures and thereby the relevant paths.

        Note that the search is bounded by the external constants w_max, i_max, and that
        it relies on the definition of formal species (fs) and non-wastes (nonw).
        """
        if not p.is_semiformal:
            return
        if p.width > w_max: 
            return 
        if sum(p.minimal_initial_state.values()) > i_max: 
            return
        if nonw is not None and not p.is_linear(nonw):
            return
        for df in p.dfs:
            if is_formal_state(df, fs):
                return

        sig = p.signature
        if sig in signatures:
            return
        signatures.add(sig)

        (i, fin, w, fc, dfs, rfs) = sig
        if len(dfs) == 0: # Ok, we found a valid path? Check if the system is tidy and regular!
            log.debug(f"New undecomposable path: {p}")
            ifin = tuple(sorted(intermediate(fin, fs)))
            if ifin not in TidyCheck or TidyCheck[ifin] is not True:
                log.debug(f"Checking if tidy from state: {ifin} ...")
                queue = [ifin] if ifin not in TidyCheck else TidyCheck[ifin]
                TidyCheck[ifin] = tidy(queue, crn, fs, TidyCheck, bound = w_max)
                if TidyCheck[ifin] is False:
                    log.info("The given system is not tidy.")
                    log.info(" - from initial state: {}".format(p.S0))
                    [log.info('    {}'.format(r)) for r in p.pretty_path]
                    raise NoFormalBasisError("The given system is not tidy.")
                elif TidyCheck[ifin] is not True:
                    log.info(f"The system is not tidy from initial state {p.S0} and final state {p.Sn} given the current width bound: {w_max}.")
                log.debug(f"... done.")
            if len(p) > 0 and p.is_formal:
                if inter: # integrated hybrid theory
                    fs1 = set(interpret(fs, inter))
                    p1 = [[interpret(R, inter), interpret(P, inter)] for [R, P] in p.lpath]
                    p1 = Path(p1, fs1)
                    i1 = p1.S0
                    t1 = p1.Sn
                    if t1 not in p1.rfs:
                        log.info("The given system is not regular.")
                        log.info(" - from initial state: {}".format(p.S0))
                        [log.info('    {}'.format(r)) for r in p1.pretty_path]
                        raise NoFormalBasisError("The given system is not regular.")
                elif fin not in rfs:
                    log.info("The given system is not regular.")
                    log.info(" - from initial state: {}".format(p.S0))
                    [log.info('    {}'.format(r)) for r in p.pretty_path]
                    raise NoFormalBasisError("The given system is not regular.")
                ebasis.append(p)
                return # adding reaction must yield a strongly decomposable pathway.

        perm = False # turn on if you want a randomized/shuffled search
        if perm:
            import copy, random
            crn_perm = copy.copy(crn) 
            random.shuffle(crn_perm)
            tmp = crn_perm
        else :
            tmp = crn
        for r in crn:
            enumerate_elementary_basis(p + Path([r]), tmp)

    ############################################################################
    #          On w and w_max as criteria to stop enumeration:                 #
    # By STW2019 Theorem 4.2, the width of an undecomposable semi-formal path  #
    # w_p and the width of that path after removing the last reaction w_{p-1}  #
    # are related by:                                                          #
    #                       w_p >= w_{p-1} >= w_p - b                          #
    # Also, p-1 can be decomposed into at most b_r undecomposable pathways     #
    # (see Thm 4.10), and at least one of those pathways must have width:      #
    #                           w >= (w-b)/b_r                                 #
    # Now take an interval between w and w_max+1, where                        #
    #                           w_max = w * b_r + b,                           #
    # then any undecomposable semi-formal pathway with width v > w_max can be  #
    # decomposed into at least one semi-formal path that has width in this     #
    # interval or (if it has w >= wmax) it can be further decomposed to have   #
    # width in this interval.                                                  #
    # Conversely, if no undecomosable semi-formal path exists in that interval,#
    # there exists no semi-formal path with width greater than w_max.          #
    ############################################################################

    w_max, i_max = 0, 0
    while True:
        log.debug("Current bounds: w_max = {}, i_max = {}".format(w_max, i_max))
        signatures, ebasis = set(), [] # We actually enumerate signatures, but care about paths.
        enumerate_elementary_basis(Path([], fs), crn)
        log.debug('After enumation of pathways: len(basis) = {}'.format(len(ebasis)))
        [log.debug(f'    {p}') for p in ebasis]
        new_w, new_i = 0, 0
        for (i, f, w, fc, dfs, rfs) in signatures:
            # NEW: only *strongly* semiformal paths (no formal paths)
            if len(dfs) == 0 and not is_formal_state(f, fs):
                log.debug(f'w: {w}, S0: {i}, Sn: {f}') 
                new_w = max(new_w, w)
                new_i = max(new_i, len(i)) # Size of initial (formal) state.
        w_t = new_w * b_r + b
        i_t = 0
        for (fr, ir) in fRiR:
            i_t = max(i_t, new_i * ir + fr)
        if w_t <= w_max and i_t <= i_max:
            for ifin in TidyCheck:
                if TidyCheck[ifin] is not True:
                    if tidy(TidyCheck[ifin], crn, fs, TidyCheck, bound = new_w) is not True:
                        raise NoFormalBasisError(f"The given system is not tidy from state {ifin}.")
            break
        w_max = w_t
        i_max = i_t

    # Now translate elementary basis into formal basis:
    fbasis_raw = []
    fbasis_int = [] # interpretation
    for path in ebasis:
        S = path.minimal_initial_state
        T = path.final_state
        r = [sorted(S.elements()), sorted(T.elements())] 
        fbasis_raw.append(r)

        if inter: # integrated hybrid theory
            fs1 = set(interpret(fs, inter))
            p1 = [[interpret(R, inter), interpret(P, inter)] for [R, P] in path.lpath]
            p1 = Path(p1, fs1)
            S1 = p1.minimal_initial_state
            T1 = p1.final_state
            r = [sorted(S1.elements()), sorted(T1.elements())] 
            fbasis_int.append(r)
    return fbasis_raw, fbasis_int

def get_crn_modules(crn, intermediates):
    """ Partition the CRN into modules that do not share intermediate species. """
    def ancestor(x):
        if parent[x] != x:
            parent[x] = ancestor(parent[x])
        return parent[x]

    # Get all the divisions
    parent   = {i: i  for i in intermediates}
    division = {i: [] for i in intermediates}

    # First, determine which intermediates occur together...
    for rxn in crn:
        t = [x for x in rxn[0] + rxn[1] if x in intermediates]
        if len(t) > 1:
            # parent of all intermediates should point to the same "ancestor".
            z = ancestor(t[0])
            for x in t[1:]:
                y = ancestor(x)
                if z != y:
                    parent[y] = z
    i = 0 # TODO: test this function!
    for rxn in crn:
        t = [x for x in rxn[0] + rxn[1] if x in intermediates]
        if len(t) > 0:
            z = ancestor(t[0])
            division[z].append(rxn)
        else: # formals, wastes only
            division[i] = [rxn]
            i += 1
    return [v for v in division.values() if len(v) > 0]

def find_basis(crn, fs, modular = True, interpretation = None):
    """Finds all formal reactions in a CRN.
    
    STW2019 - Def13: The set of prime pathways in a given CRN is called the
    *elementary basis* of the CRN. The *formal basis* is the set of (initial
    state, final state) pairs of the pathways in the elementatry basis.

    Args:
      crn (List[[R], [P]]): A list of irreversible reactions. First a list of
            reactants, second a list of products.
      fs (set(str)): A set of formal species.
      modular (bool, optional): Chop the CRN into modules and then find the basis
            separately for each of those modules. Defaults to True.
      interpretation (dict, optional): An interpretation dictionary for the 
            "integrated hybrid" mode. Defaults to None.

    Returns:
        basis_raw, basis_int: The formal basis found without and with an 
            interpretation dictionary. 
    """
    if modular:
        intermediates, wastes, nonwastes = assign_crn_species(crn, fs)
        divs = sorted(get_crn_modules(crn, intermediates), key = lambda x: len(x))
        log.info(f"Divided the implementation CRN into {len(divs)} modules " + \
                f"with {[len(d) for d in divs]} reactions.")
        basis_raw = []
        basis_int = []
        for e, mod in enumerate(divs, 1):
            log.info("Verifying module {}:".format(e))
            [log.info(f'    {r}') for r in pretty_crn(mod)]
            b_raw, b_int = get_formal_basis(mod, fs, inter = interpretation)
            log.info("Formal basis of the current module:")
            [log.info('    {}'.format(r)) for r in pretty_crn(b_raw)]
            if interpretation is not None:
                log.info("Formal basis of the current module (integrated hybrid):")
                [log.info('    {}'.format(r)) for r in pretty_crn(b_int)]
            log.info('')
            basis_raw += b_raw
            basis_int += b_int
    else:
        basis_raw, basis_int = get_formal_basis(crn, fs, inter = interpretation)

    return clean_crn(basis_raw), clean_crn(basis_int)

def my_parse_crn(string, is_file = False):
    crn, species = parse_crn_file(string) if is_file else parse_crn_string(string)
    crn = split_reversible_reactions(crn)
    crn = [list(rxn[:2]) for rxn in crn]
    return crn, set(species.keys())

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help="print verbose output. -vv increases verbosity level.")
    parser.add_argument("--old", action = 'store_true',
            help="""Use original algo.""")
    parser.add_argument("--crn-file", action='store', metavar='</path/to/file>',
            help="""Read a CRN from a file.""")
    parser.add_argument("--formal-species", nargs = '+', default = [], 
            action = 'store', metavar = '<str>', 
            help="""List formal species in the CRN.""")
    parser.add_argument("--fuel-species", nargs = '+', default = [], 
            action = 'store', metavar = '<str>', 
            help="""List fuel species in the CRN.""")
    parser.add_argument("--non-modular", action = 'store_true',
            help="""Do not attempt to split into smaller CRN modules.""")
    parser.add_argument("--integrated", action = 'store_true',
            help="""Use interpretation when finding formal basis.""")
    parser.add_argument("--profile", action = 'store_true',
            help="""Get some code profiling information (requires statprof-smarkets).""")
    args = parser.parse_args()

    if args.profile:
        try:
            import statprof
        except ImportError as err:
            print('Cannot import statprof module.')
            args.profile = False

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

    crn, fs = my_parse_crn(args.crn_file, is_file = True)

    if args.fuel_species:
        crn = remove_const(crn, args.fuel_species)

    if args.formal_species:
        fs &= set(args.formal_species)

    log.info('Input formal species: {}'.format(fs))
    log.info('Input fuel species: {}'.format(args.fuel_species))
    log.info('Input CRN (without fuel species):')
    [log.info('    ' + r) for r in pretty_crn(crn)]

    i, wastes, nw = assign_crn_species(crn, fs)
    assert i == nw
    #wastes = set()
    log.info('{} waste species are treated as formal. ({})'.format(len(wastes), ' '.join(wastes)))

    inter = {}
    for x in fs: inter[x] = [x]
    for x in wastes: inter[x] = []

    try:
        if args.profile:
            statprof.start()
        basis_raw, basis_int = find_basis(crn, fs | wastes, 
                                modular = not args.non_modular,
                                interpretation = inter if args.integrated else None)
        if args.profile:
            statprof.stop()

        if basis_raw == []:
            print("The formal basis is an empty set of reactions.")
        else:
            print("Formal basis:")
            for r in pretty_crn(basis_raw):
                print(r)
            if args.integrated:
                print("Integrated hybrid basis:")
                for r in pretty_crn(basis_int):
                    print(r)
    except NoFormalBasisError as err:
        print("Could not find formal basis: {}".format(err))

    if args.profile:
        statprof.display()
