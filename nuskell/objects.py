# -*- coding: utf-8 -*-
#
# Written by Stefan Badelt (badelt@caltech.edu)
#
# nuskell.objects: shared between different components of nuskell
#
#from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import networkx as nx
from collections import Counter
from sympy import Symbol, sympify, Matrix

from dsdobjects import clear_memory
from dsdobjects import DL_Domain, DSD_Complex, DSD_Reaction
from dsdobjects import DSDObjectsError, DSDDuplicationError

#'http://www.github.com/bad-ants-fleet/crnsimulator'
from crnsimulator import writeODElib

from dsdobjects.parser import parse_kernel_file


class NuskellObjectError(Exception):
    """nuskell.objects error class."""

    def __init__(self, msg):
        self.message = msg
        super(NuskellObjectError, self).__init__(self.message)

class NuskellDomain(DL_Domain):
    """Nucleic acid domain sequence.

    Inherits from :obj:`dsdobjects.base_classes.DL_Domain()`.

    A domain is a sequence of consecutive nucleotides. Sequences of domains can
    form the same secondary structures as sequences of nucleotides.  Several
    options for specifying domain properties are allowed. Domains might have an
    explicit integer (bp) length, or may be designated as short or long. If the
    latter method is used, the code will use the relevant constant as the
    integer domain length.

    Globals:
      ID (int): An automatically assigned ID that is used to find automated
      species names.

    Args:
      name (str, optional): Name of this domain. If not specified, an automatic
        name is generated, defaults to ''
      prefix (str, optional): A prefix for automated naming of Domains. Has to be
        set if no name is provided. Defaults to 'd' for 'domain'. Usually, the
        prefix 't' is used for 'toehold domains', and 'h' for 'histoy domains'.
      dtype (bool, optional): One of *short* or *long*. Defaults to None, i.e. 
        guess from length parameter.
      length (int, optional): Set length of a domain. Defaults to None, i.e. set
        default length for *dtype* parameter.

    Raises:
      NuskellObjectError: Domain prefix must not be empty!
      NuskellObjectError: NuskellDomain prefix must not be empty!
      NuskellObjectError: NuskellDomain prefix must not end with a digit!

    """

    ID = 0          # ID is used to assign names automatically

    def __init__(self, name='', prefix='d', dtype=None, length=None):
        # Assign name
        if name == '':
            if prefix == '':
                raise NuskellObjectError('NuskellDomain prefix must not be empty!')
            elif prefix[-1].isdigit():
                raise NuskellObjectError('NuskellDomain must not end with a digit!')
            name = prefix + str(NuskellDomain.ID)
            NuskellDomain.ID += 1
        super(NuskellDomain, self).__init__(name, dtype, length)

    @property
    def complement(self):
        """
        Returns the complement of a Domain, initializes a new object if necessary.
        """
        if self._complement is None:
            cname = self._name[:-1] if self.is_complement else self._name + '*'
            if cname in DL_Domain.MEMORY:
                self._complement = DL_Domain.MEMORY[cname]
            else :
                self._complement = NuskellDomain(name=cname, dtype=self.dtype, length=self.length)
        return self._complement

    @property
    def is_complement(self):
        """
        Returns true if this domain is a complement (e.g. A* rather than A),
        false otherwise.
        """
        return self._name[-1:] == '*'

class NuskellComplex(DSD_Complex):
    """ A sequence and structure pair.

    Inherits from :obj:`dsdobjects.base_classes.Complex()`.

    Args:
        sequence (list): A domain-level or nucleotide-level sequence.
        structure (list): A domain-level or nucleotide-level dot-bracket notation.
        name (str, optional): Name of this domain. If not specified, an automatic
            name is generated.
        prefix (str, optional): A prefix for automatic naming of Domains.
            Defaults to 'cplx'.
        memorycheck (bool, optional): Use built-in memory checks. Defaults to True.

    """

    def __init__(self, sequence, structure, name='', prefix='cplx', memorycheck=True):
        super(NuskellComplex, self).__init__(sequence, structure, name, prefix, memorycheck)

    @staticmethod
    def clear_memory(memory=True, names=True, ids=True):
        if memory :
            DSD_Complex.MEMORY = dict()
        if names :
            DSD_Complex.NAMES = dict()
        if ids: 
            DSD_Complex.ID = dict()

class NuskellReaction(DSD_Reaction):
    """ A reaction pathway.

    Args:
        reactants (list): A list of :obj:`NuskellComplex()` objects.
        products (list): A list of :obj:`NuskellComplex()` objects.
        rtype (str, optional): Reaction type. Must be one of:
            'formal', 'condensed', 'open', 'bind11', 'bind21', 'branch-3way', 'branch-4way'.
        rate (flt, optional): Reaction rate.
        memorycheck (bool, optional): Use built-in memory checks. Defaults to True.
    """
    RTYPES = set(['formal', 'condensed', 'open', 'bind11', 'bind21', 'branch-3way', 'branch-4way'])

    def __init__(self, reactants, products, rtype=None, rate=None, memorycheck=True):
        super(NuskellReaction, self).__init__(reactants, products, 
                rtype=rtype, rate=rate, memorycheck=memorycheck)
        if self._rtype not in NuskellReaction.RTYPES:
            try:
                del DSD_Reaction.MEMORY[self.canonical_form]
            except KeyError:
                pass
            raise NuskellObjectError('Reaction type {} not supported! '.format(self.rtype) + 
                'Set supported reaction types using NuskellReaction.RTYPES')

    @staticmethod
    def clear_memory():
        DSD_Reaction.MEMORY = dict()

class TestTube(object):
    """A reaction network of nucleic acid complexes.

    **Description:**
      :obj:`TestTube()` objects are Nuskell's interface to enumerate and simulate
      nucleic acid systems.  Domain-level reaction networks are enumerated using
      the Python package `peppercornenumerator`_ and simulated using the Python package
      `crnsimulator`_. :obj:`TestTube()` objects to low-level data structures,
      e.g. to verify the equivalence between two :obj:`TestTube()` objects.

      Single or multiple :obj:`NuskellComplex()` and/or :obj:`NuskellReaction()` objects can be
      accessed, added and removed from the system.  :obj:`TestTube()` provides
      (optional) assert statements checking if :obj:`NuskellComplex()` and
      :obj:`NuskellDomain()` instances have been duplicated, but they might be
      time-consuming for large networks.

      Built-in functions can process domain-level networks to remove strands with
      wildcard-domains (also called history-domains).  Together with the
      remainder of the molecule, a species with a wildcard forms a
      regular-expression, matching every other species in the system that differs
      only by a single domain instead of ‘?’.  If there exists a species matching
      the regular expression, then the species with the wildcard domain and every
      enumerated reaction emerging from that species is be removed from the
      system, otherwise, the wildcard domain is replaced by a regular long
      domain.

    **Structure:**
      :obj:`TestTube()` is based on networkx.MultiDiGraph(), with two types of nodes: (a)
      :obj:`NuskellComplex()` and (b) :obj:`NuskellReaction()`.  The graph is bipartite, edges
      are directed (from reactants to prodcuts), they connect reactants to a
      reaction node and a reaction node to products.

      :obj:`TestTube()` provides an additional *concentration* attribute and *constant*
      atribute to :obj:`NuskellComplex()` nodes, as well as a rate attribute to
      :obj:`NuskellReaction()` nodes. These attributes are accessed when writing ODE
      systems and (optionally) updated after simulations. This means a
      :obj:`TestTube()` *can* be initialized without using :obj:`NuskellComplex()` and
      :obj:`NuskellReaction()` objects, but by using a consistent naming scheme.

    **Developers:**
      It is recommended to store all other node attributes (e.g. free energies,
      etc.) directly in the :obj:`NuskellComplex()` and :obj:`NuskellReaction()` objects.
      :obj:`TestTube()` does not provide an I/O interface for file formats. There is a
      separate :obj:`TestTubeIO()` object explicitly to parse and write compatible file
      formats (\*.pil, \*.dom, \*.dna, etc.).

    Args:
      complexes (:obj:`dict`, optional) =  A dictionary of complex names that stores
        a tuple of [0]: the respective :obj:`NuskellComplex()` object or :obj:`None`, [1]: the
        concentration, and [2]: boolean value to specify if the concentration
        stays constant during a simulation: True = constant; False = initial
        For example: complexes['A'] = (None, 100, True) is a new
        species called 'A', which has no :obj:`NuskellComplex()` object initialized, and
        remains constant at 100 [nM] concentration.
        Note: The constant attribute defaults to False if not specified, in which
        case the concentrations specify the initial state and might be updated
        after a simulation.

      reactions <dict(), optional> =  A dictionary of reaction names that stores
        either a :obj:`NuskellReaction()` object or a list [[reactants], [products], rate]
        Reactants and products have to be known complexes.

    Raises:
      NuskellObjectError: 'Wrong initialization of arguments for TestTube()'
      NuskellObjectError: 'Invalid Reaction format'


    .. _peppercornenumerator:
        http://www.github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator
    .. _crnsimulator:
        http://www.github.com/bad-ants-fleet/crnsimulator

    """

    # A global TestTube() variable to make (time-consuming) sanity-checks when
    # separate instances are subject to boolean or arithmetic opertions.
    sanitychecks = True
    warnings = True

    def __init__(self, complexes=None, reactions=None):
        self._RG = nx.MultiDiGraph()

        if complexes and isinstance(complexes, dict):
            for name, data in complexes.items():
                assert isinstance(data, list) and len(data) == 3, 'wrong format for complexes'
                if data[1] is not None :
                    assert isinstance(data[1], float), 'concentration must be of type float.'
                if data[2] is not None :
                    assert isinstance(data[2], bool), '"constant" must be of type bool.'

                cplx = data[0]
                conc = data[1]
                const = data[2]
                if isinstance(cplx, NuskellComplex):
                    self._RG.add_node(cplx, concentration=conc, constant=const)
                else:
                    raise NotImplementedError('please specify NuskellComplex objects.')
                    #self._RG.add_node(name, concentration=conc, constant=const)
        elif complexes and isinstance(complexes, list):
            for c in complexes:
                assert isinstance(c, NuskellComplex)
                self._RG.add_node(c, concentration=None, constant=None)

        if reactions:
            for react in reactions:
                assert isinstance(react, NuskellReaction)
                self._RG.add_node(react, rate=react.rate)
                for r in react.reactants:
                    assert self._RG.has_node(r)
                    self._RG.add_edge(r, react)
                for p in react.products:
                    assert self._RG.has_node(p)
                    self._RG.add_edge(react, p)

        #self._complexes = None
        self._domains = None
        self._strands = None

    @property
    def ReactionGraph(self):
        """:obj:`networkx.MultiDiGraph()`: bipartite reaction graph. """
        return self._RG

    @ReactionGraph.setter
    def ReactionGraph(self, RG):
        self._RG = RG
        self._domains = None
        self._strands = None

    @property
    def complexes(self):
        """list: a list of :obj:`NuskellComplex()` objects. """
        return [n for n in self._RG.nodes() if isinstance(n, NuskellComplex)]

    def has_complex(self, cplx):
        return self._RG.has_node(cplx)

    @property
    def reactions(self):
        """list: a list of :obj:`NuskellReaction()` objects. """
        return [n for n in self._RG.nodes() if isinstance(n, NuskellReaction)]

    def set_complex_concentration(self, cplx, concentration, constant):
        """
        Args:
          cplx (:obj:`NuskellComplex()`): complex.
          concentration (flt): concentration.
          constant (bool): True if the concentration is kept constant, False for
            initial concentration.
        """
        self._RG.node[cplx]['concentration'] = concentration
        self._RG.node[cplx]['constant'] = constant

    def get_complex_concentration(self, cplx):
        """flt, bool: First value is concentration, second is True if it is constant,
        False if it is variable."""
        return self._RG.node[cplx]['concentration'], self._RG.node[cplx]['constant']

    def selected_complexes(self, names):
        """list: a list of :obj:`NuskellComplex()` objects that correspond to specified names. """
        return [n for n in self._RG.nodes() if isinstance(n, NuskellComplex)
                and n.name in names]

    # TODO: exclude? or filter outside?
    def present_complexes(self, exclude=[], th=0):
        """Returns a list of :obj:`NuskellComplex()` objects with occupancy greater than a threshold.

        Args:
          exclude (list, optional): Exclude particular complex names.
          th (flt, optional): Specify the threshold to consider complexes as
            *present* in solution. Defaults to 0.
        """
        return [n for n, att in self._RG.node.items() 
                if isinstance(n, NuskellComplex) and 
                    att['concentration'] > th and n.name not in exclude]

    def interpret_species(self, species, prune=True):
        """Get an interpretation dictionary.

        If a :obj:`NuskellComplex()` sequence contains a wildcard, then this function will find
        all matching complexes, and return those as interpretation.  Regex-nodes
        may have at most *one wildcard* per complex, a wildcard corresponds to
        exactly *one unpaired long domain*.

        Args:
          species (list[str], optional): A list of complex names that are potential
            regular-expression complexes.
          prune (bool, optinal): True: Remove all wildcards from the network. If a
            matching complex has been found, then the regex-complex and all emerging
            complexes are removed, if no matching complex has been found, then the
            wildcard domain is replaced by a regular long domain. False: Return the
            interpretation without udating the reaction network. Defaults to True.

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
          [:obj:`dict()`]: Interpretation of signal species: dict['A_i'] = Counter('A':1)
        """

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

            pMx = [map(str, x.sequence), map(str, x.structure)]
            pMy = [map(str, y.sequence), map(str, y.structure)]
            if pM_check(pMx, pMy):
                return True

        def get_matching_complexes(regex, hist):
            """Find all matching complexes. """
            regseq = regex.sequence
            regstr = regex.structure

            matching = []
            for cplx in self.complexes :
                if regex.name == cplx.name :
                    continue
                elif patternMatch(regex, cplx, ignore=hist):
                    matching.append(cplx)
            return matching

        need_to_prune = False
        interpretation = dict()
        for fs in species:
            cplxs = self.selected_complexes([fs])
            if len(cplxs) == 0:
                print '====='
                print 'WARNING: No complex found with name of formal species:', fs
                print '====='
                continue
            else:
                assert len(cplxs) == 1, 'Duplicate complex names?'

            cplx = cplxs[0]
            hist = filter(lambda x: x[0] == 'h', map(str, cplx.sequence))
            if hist :
                assert len(hist) == 1, 'no support for multiple history domains'
                hist = hist[0]
                matches = get_matching_complexes(cplx, hist)
                if matches:
                    need_to_prune = True
                    for e, m in enumerate(matches, 1):
                        m.name = fs + '_' + str(e) + '_'
                        interpretation[m.name] = Counter([fs])
                    self.rm_complex(cplx, force=True)
                else:
                    # NOTE: We cannot simply remove the domain, because we would need to
                    # enumerate the network again and remove the domain everywhere! So
                    # unless we enumerate up-front with a history-pruned species, this
                    # gets us into trouble.
                    interpretation[cplx.name] = Counter([fs])
            else:
                interpretation[cplx.name] = Counter([fs])

        if prune and need_to_prune:
            # Get rid of all reactions with history wildcards. Start with a set
            # of produce molecules and see what species emerge from reactions
            # consuming these molecules.
            # Alternative: enumerate again using history-replaced species.
            rxns = self.reactions
            [prev, total] = [set(), set(interpretation.keys() + map(str, self.present_complexes()))]
            while prev != total:
                prev = set(list(total)) # force a copy?
                for rxn in rxns:
                    self.rm_reaction(rxn)
                    r = map(str, rxn.reactants)
                    p = map(str, rxn.products)
                    if set(r).intersection(total) == set(r):
                        total = total.union(set(p))
            map(self.add_reaction, filter(lambda x: set(map(str, x.reactants)).intersection(
                    total) == set(map(str, x.reactants)), rxns))

            # Now remove all the left-over complexes from the graph.
            all_nodes = set(self.complexes)
            assert set(map(str, all_nodes)).issuperset(total)
            total = set([n for n in self.complexes if str(n) in total])
            remove = all_nodes.difference(total)
            self._RG.remove_nodes_from(remove)

        return interpretation

    def add_complex(self, cplx, conctup = (None, None), sanitycheck=True):
        """Add a complex to the TestTube.

        Args:
          cplx (:obj:`NuskellComplex()`): The complex object.
          (conc, const) (flt, bool): Concentration and True/False for constant or
            initial concentrations.
          sanitycheck (bool): True: Check if complex exists under a different name.
            This can be time consuming. Defaults to True.

        Note:
          A new complex resets TestTube.domains and TestTube.strands

        """
        (conc, const) = conctup
        assert isinstance(cplx, NuskellComplex), 'must be a NuskellComplex format'

        if self._RG.has_node(cplx): # Does not check for the name
            if conc is not None:
                if self._RG.node[cplx]['concentration'] is None:
                    self._RG.node[cplx]['concentration'] = conc
                else:
                    assert self._RG.node[cplx]['concentration'] == conc, \
                            "Conflicting complex concentrations"
            if const is not None:
                if self._RG.node[cplx]['constant'] is None:
                    self._RG.node[cplx]['constant'] = const
                else:
                    assert self._RG.node[cplx]['constant'] == const, \
                            "Conflicting complex concentrations"
        else:
            # NOTE: This might become inefficient at some point, but it has been
            # introduced to overcome issues with some translation schemes that
            # produce the same fuel strand multiple times.
            if sanitycheck and cplx.canonical_form in map(
                    lambda x: x.canonical_form, self.complexes):
                raise NuskellObjectError('trying to add duplicate complex: {} = {}'.format(
                    cplx.name, cplx.kernel_string))
            else:
                self._RG.add_node(cplx, concentration=conc, constant=const)
                self._domains = None
                self._strands = None

    def rm_complex(self, cplx, force=False):
        """Remove a Complex from the TestTube.

        Args:
          cplx (:obj:`NuskellComplex()`): The complex object.
          force (bool): True: remove complex and all reactions it is engaged in.
            False: Raise an Error if complex is engaged in a reaction.

        Raises:
          RuntimeError: Cannot remove a complex engaged in reactions.

        Note:
          Removing a complex resets TestTube.domains and TestTube.strands
        """
        if self._RG.has_node(cplx):
            if force:
                for (r, c) in self._RG.in_edges(cplx):
                    assert isinstance(r, NuskellReaction)
                    self.rm_reaction(r)
                for (c, r) in self._RG.out_edges(cplx):
                    assert isinstance(r, NuskellReaction)
                    self.rm_reaction(r)
            elif (self._RG.in_edges(cplx) or self._RG.out_edges(cplx)):
                raise NuskellObjectError(
                    "Cannot remove a complex engaged in reactions.")
            self._RG.remove_node(cplx)
            self._domains = None
            self._strands = None

    def add_reaction(self, react, sanitycheck=True):
        """Add a reaction to the TestTube.

        Args:
          react (:obj:`NuskellReaction()`): The *irreversible* reaction to be added.
          sanitycheck (bool): True: Check if reaction exists under a different name.
            This can be time consuming. Defaults to True.
        """

        assert isinstance(react, NuskellReaction) 

        if self._RG.has_node(react):
            assert self._RG.node[react]['rate'] == react.rate
        else:
            # NOTE: This might become inefficient at some point, but there might be
            # cases where reactions are duplicated, so we check if the very same
            # reaction exists as a different node:
            if sanitycheck and react.canonical_form in map(
                    lambda x: x.canonical_form, self.reactions):
                raise NuskellObjectError('trying to add duplicate reaction:', react.kernel_string)
            else:
                self._RG.add_node(react, rate=react.rate)
                for r in react.reactants:
                    assert self._RG.has_node(r)
                    self._RG.add_edge(r, react)
                for p in react.products:
                    assert self._RG.has_node(p)
                    self._RG.add_edge(react, p)

    def rm_reaction(self, react):
        """Remove a reaction from the TestTube.

        Args:
          react (:obj:`NuskellReaction()`): The reaction object to be removed.
        """
        if self._RG.has_node(react):
            self._RG.remove_node(react)

    @property
    def strands(self):
        """Return a list of strands present in the TestTube.

        A strand is a nucleic-acid molecule connected by a single covalent
        backbone. Strands are named automatically, and their names may change
        whenever a new Complex is added to the TestTube.

        Returns:
          [:obj:`dict()`]: strands[strand_1] = [Domain(X), Domain(Y), Domain(Z)]
        """
        if not self._strands:
            count = 0
            self._strands = dict()
            self._strand_names = dict()
            for cplx in self.complexes:
                for s in cplx.lol_sequence:
                    strand = tuple(map(str, s))
                    if strand not in self._strand_names:
                        name = 'strand_{}'.format(count)
                        count += 1
                        self._strand_names[strand] = name
                        self._strands[name] = s
        return self._strands

    @property
    def domains(self):
        """Return a dictionary of Domain Objects present in the TestTube.

        Returns:
          [:obj:`dict()`]: domains[Domain.name] = Domain
        """
        if not self._domains:
            self._domains = set()
            for cplx in self.complexes:
                for d in cplx.domains:
                    self._domains.add(d)
        return list(self._domains)

    def enumerate_reactions(self, args=None, condensed=True, rename=None, 
            prefix='e', init_memory=None):
        """Enumerate reactions using the *peppercorn* enumerator.
        Args:
          args(:obj:`argparse.ArgumentParser()`, optional): Arguments for *peppercorn*.
          condensed (bool, optional): Udate the reaction graph using *condensed* format.
        """

        from nuskell.enumeration import TestTubePeppercornIO

        TestTubePeppercornIO.condensed = condensed
        interface = TestTubePeppercornIO(testtube=self, enumerator=None,
                                         pargs=args, rename=rename, prefix=prefix,
                                         init_memory = init_memory)
        interface.enumerate()
        self.ReactionGraph = nx.compose(self.ReactionGraph, interface.testtube.ReactionGraph)

    def simulate_crn(self, odename, sorted_vars=None, unit='M'):
        oR = dict()
        conc = dict()
        ode = dict()
        for r in self._RG.nodes_iter():
            if isinstance(r, NuskellComplex):
                concentration = self._RG.node[r]['concentration']
                const = self._RG.node[r]['constant']
                if concentration == float('inf'):
                    concentration = 100 * 1e-9
                elif concentration is None:
                    concentration = 0.

                if unit == 'M':
                    pass
                elif unit == 'mM':
                    concentration *= 1e3
                elif unit == 'uM':
                    concentration *= 1e6
                elif unit == 'nM':
                    concentration *= 1e9
                else:
                    raise NuskellObjectError(
                        'Concentration unit not supported', unit)

                conc[str(r)] = concentration
                continue

            rate = 'k' + str(len(oR.keys()))
            if unit == 'M':
                oR[rate] = str(r.rate)
            elif unit == 'mM':
                if r.arity[0] > 1:
                    factor = r.arity[0] - 1
                    oR[rate] = str(float(r.rate) / (factor * 1e3))
                else:
                    oR[rate] = str(r.rate)
            elif unit == 'uM':
                if r.arity[0] > 1:
                    factor = r.arity[0] - 1
                    oR[rate] = str(float(r.rate) / (factor * 1e6))
                else:
                    oR[rate] = str(r.rate)
            elif unit == 'nM':
                if r.arity[0] > 1:
                    factor = r.arity[0] - 1
                    oR[rate] = str(float(r.rate) / (factor * 1e9))
                else:
                    oR[rate] = str(r.rate)
            else:
                raise NuskellObjectError(
                    'concentration unit not supported', unit)

            reactants = []
            for reac in self._RG.predecessors_iter(r):
                for i in range(self._RG.number_of_edges(reac, r)):
                    reactants.append(Symbol(str(reac)))

            products = []
            for prod in self._RG.successors_iter(r):
                for i in range(self._RG.number_of_edges(r, prod)):
                    products.append(Symbol(str(prod)))

            for x in reactants:
                if x in ode:
                    ode[x].append(['-' + rate] + reactants)
                else:
                    ode[x] = [['-' + rate] + reactants]

            for x in products:
                if x in ode:
                    ode[x].append([rate] + reactants)
                else:
                    ode[x] = [[rate] + reactants]

        if sorted_vars:
            assert len(sorted_vars()) == len(ode.keys())
            oV = map(Symbol, sorted_vars)
        else:
            oV = sorted(ode.keys(), key=lambda x: str(x))
            oC = map(lambda x: conc[str(x)], oV)

        # Sympy Symbol namespace
        ns = dict(zip(map(str, oV), oV))

        oM = []
        for dx in oV:
            sfunc = sympify(
                ' + '.join(['*'.join(map(str, xp)) for xp in ode[dx]]), locals=ns)
            ode[dx] = sfunc
            oM.append(sfunc)

        oM = Matrix(oM)
        oJ = None

        oFile, oname = writeODElib(
            oV, oM, jacobian=oJ, rdict=oR, concvect=oC, filename=odename)
        return oFile, oname

    def __add__(self, other):
        assert isinstance(other, TestTube)
        combined = TestTube()

        # global TestTube() variable
        if not TestTube.sanitychecks:
            if TestTube.warnings:
                print Warning('TestTube() - sanity checks turned off!')
            combined.ReactionGraph = nx.compose(
                self.ReactionGraph, other.ReactionGraph)

        elif len(other.complexes) > len(self.complexes):
            combined.ReactionGraph.add_nodes_from(
                other.ReactionGraph.nodes(data=True))
            combined.ReactionGraph.add_edges_from(
                other.ReactionGraph.edges(data=True))
            map(lambda c: combined.add_complex(c, self.get_complex_concentration(c),
                                               sanitycheck=True), self.complexes)
            map(lambda r: combined.add_reaction(r, sanitycheck=True), self.reactions)
        else:
            combined.ReactionGraph.add_nodes_from(
                self.ReactionGraph.nodes(data=True))
            combined.ReactionGraph.add_edges_from(
                self.ReactionGraph.edges(data=True))
            map(lambda c: combined.add_complex(c, other.get_complex_concentration(c),
                                               sanitycheck=True), other.complexes)
            map(lambda r: combined.add_reaction(r, sanitycheck=True), other.reactions)
        return combined

    def __radd__(self, other):
        # Reverse add is used for: sum([Testtube1, Testtube2, ...])
        if other == 0:
            return self
        else:
            return self.__add__(other)

class TestTubeIO(object):
    """A wrapper class to handle I/O of TestTube objects.

    Args:
      ttube (obj:`TestTube()`): A :obj:`TestTube()` object that should be initialized
        or written to a text format.

    """

    def __init__(self, ttube):
        assert isinstance(ttube, TestTube)
        self._testtube = ttube

    @property
    def testtube(self):
        """:obj:`TestTube()` property."""
        return self._testtube

    def write_pil_kernel(self, pil, unit='M', crn=None, ts=None):
        """Write the contents of :obj:`TestTube()` into a PIL file -- KERNEL notation).

        Args:
          pil (filehandle): A filehandle that the output is written to.
          unit (str, optional): Specify a unit of concentrations (M, mM, uM, nM, pM).
          crn (list[list], optional): a nuskell-style CRN expression
          ts (str, optional): name of the translation scheme

        Example:
          length d1 = 6
          length d2 = 4
          length h3 = 1
          cplx1 = h3 d1( d2( + )) @ initial 10 nM
        """
        pil.write("# File autogenerated by nuskell. ")

        if ts:
            pil.write("\n# - Translation Scheme: {}".format(ts))
        if crn:
            pil.write("\n# - Input CRN: \n")
            for rxn in crn:
                assert len(rxn) == 3
                if len(rxn[2]) == 2:
                    pil.write("#    {} <=> {}\n".format(
                        ' + '.join(rxn[0]), ' + '.join(rxn[1])))
                else:
                    pil.write("#    {} -> {}\n".format(
                        ' + '.join(rxn[0]), ' + '.join(rxn[1])))
        pil.write("#\n\n".format(crn))

        domains = self._testtube.domains

        def adjust_conc(conc, unit):
            units = ['M', 'mM', 'uM', 'nM', 'pM']
            # 0,  3,   6,   9,   12
            assert unit in units
            mult = units.index(unit) * 3
            return conc * (10**mult), unit

        # Print Domains
        pil.write("# Domain Specifications\n")
        seen = set()
        for d in sorted(domains, key=lambda x: x.name):
            if d.is_complement:
                dom = ~d
            else :
                dom = d
            if dom not in seen:
                pil.write("length {:s} = {:d}\n".format(dom.name, dom.length))
                seen.add(dom)

        pil.write("\n# Complex Specifications\n")

        # Print Complexes
        for cplx in sorted(self._testtube.complexes, key=lambda x: str(x)):
            pil.write("{:s} = ".format(cplx.name))
            seq = cplx.sequence
            sst = cplx.structure
            for i in range(len(seq)):
                if sst[i] == '+':
                    pil.write("{:s} ".format(str(sst[i])))
                elif sst[i] == ')':
                    pil.write("{:s} ".format(str(sst[i])))
                elif sst[i] == '(':
                    pil.write("{:s} ".format(str(seq[i]) + str(sst[i])))
                else:
                    pil.write("{:s} ".format(str(seq[i])))

            conc, const = self._testtube.get_complex_concentration(cplx)
            if const is True:
                pil.write(" @ constant {} {}".format(*adjust_conc(conc, unit)))
            elif const is False:
                pil.write(" @ initial {} {}".format(*adjust_conc(conc, unit)))
            pil.write(" \n")

    def load_pil_kernel(self, pilfile):
        """Parses a file written in PIL - KERNEL notation! """
        ppil = parse_kernel_file(pilfile)

        def resolve_loops(loop):
            """ Return a sequence, structure pair from kernel format with parenthesis. """
            sequen = []
            struct = []
            for dom in loop:
                if isinstance(dom, str):
                    sequen.append(dom)
                    if dom == '+':
                        struct.append('+')
                    else:
                        struct.append('.')
                elif isinstance(dom, list):
                    struct[-1] = '('
                    old = sequen[-1]
                    se, ss = resolve_loops(dom)
                    sequen.extend(se)
                    struct.extend(ss)
                    sequen.append(old + '*' if old[-1] != '*' else old[:-1])
                    struct.append(')')
            return sequen, struct

        # Do domains first, just in case...
        domains = {'+' : '+'} # saves some code
        for line in ppil :
            name = line[1]
            if line[0] == 'dl-domain':
                if line[2] == 'short':
                    (dtype, dlen) = ('short', None)
                elif line[2] == 'long':
                    (dtype, dlen) = ('long', None)
                else :
                    (dtype, dlen) = (None, int(line[2]))
                if name not in domains:
                    domains[name] = NuskellDomain(name, dtype = dtype, length = dlen)
                logging.info('Domain {} with length {}'.format(domains[name], len(domains[name])))
                cname = name[:-1] if domains[name].is_complement else name + '*'
                if cname in domains:
                    assert domains[cname] == ~domains[name]
                else :
                    domains[cname] = ~domains[name]

        complexes = {}
        for line in ppil :
            name = line[1]
            if line[0] == 'dl-domain':
                pass
            elif line[0] == 'complex':
                sequence, structure = resolve_loops(line[2])

                # Replace names with domain objects.
                try :
                    sequence = map(lambda d : domains[d], sequence)
                except KeyError:
                    for e, d in enumerate(sequence):
                        if d not in domains :
                            logging.warning("Assuming {} is a long domain.".format(d))
                            domains[d] = PepperDomain(d, 'long')
                            cdom = ~domains[d]
                            domains[cdom.name] = cdom
                        sequence[e] = domains[d]

                constant, concentration = False, 0
                if len(line) > 3:
                    i, c, u = line[3]
                    constant = (i == 'constant')
                    if u == 'M':
                        concentration = float(c)
                    elif u == 'mM':
                        concentration = float(c)*1e-3
                    elif u == 'uM':
                        concentration = float(c)*1e-6
                    elif u == 'nM':
                        concentration = float(c)*1e-9
                    elif u == 'pM':
                        concentration = float(c)*1e-12
                    else :
                        raise ValueError('unknown unit for concentrations specified.')

                complexes[name] = NuskellComplex(sequence, structure, name=name)
                self._testtube.add_complex(complexes[name], (concentration, constant))

            elif line[0] == 'reaction':
                rtype = line[1][0][0] if line[1] != [] and line[1][0] != [] else None
                rate = float(line[1][1][0]) if line[1] != [] and line[1][1] != [] else None
                if rate is None or rtype is None or rtype == 'condensed' :
                    r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
                    logging.warning("Ignoring input reaction without a rate: {}".format(r))
                    continue
                else :
                    r = "[{} = {:12g}] {} -> {}".format(
                            rtype, rate, ' + '.join(line[2]), ' + '.join(line[3]))
                    logging.warning("Ignoring input reaction: {}".format(r))
                    continue
                #try :
                #    reactants = map(lambda c : complexes[c], line[2])
                #    products  = map(lambda c : complexes[c], line[3])
                #except KeyError:
                #    logging.warning("Ignoring input reaction with undefined complex: {}".format(r))
                #    continue

                #reaction = PepperReaction(reactants, products, rtype=rtype, rate=rate)
                #reactions.append(reaction)

            elif line[0] == 'resting-state':
                logging.warning("Ignoring resting-state specification: {}".format(name))
            else :
                raise NotImplementedError('cannot interpret keyword:', line[0])

        return self._testtube

    def write_dnafile(self, fh, signals=[], crn=None, ts=None):
        """ Write a TestTube Object into VisualDSD \*.dna format.

        Note:
          This function assumes that toehold domains are named starting with a 't',
          history domains start with a 'h' and anti-sense domains end with '*'.

        Args:
          fh (filehandle): The function prints to this filehandle.
          signals (list[str], optional): A list of signal species.
          crn (list[list], optional): a nuskell-style CRN expression
          ts (str, optional): name of the translation scheme
        """

        def pair_table(ss, chars=['.']):
            """Return a secondary struture in form of pair table:
        
            Args:
              ss (str): secondary structure in dot-bracket format
              chars (list, optional): a list of characters that are ignored. Defaults to
                ['.']
        
            Example:
               ((..)). => [5,4,-1,-1,1,0,-1]
        
            Raises:
               NuskellObjectError: Too many closing brackets in secondary structure.
               NuskellObjectError: Too many opening brackets in secondary structure.
               NuskellObjectError: Unexpected character in sequence: "{}"
        
            Returns:
              [list]: A pair-table
            """
            stack = []
        
            pt = [-1] * len(ss)
        
            for i, char in enumerate(ss):
                if (char == '('):
                    stack.append(i)
                elif (char == ')'):
                    try:
                        j = stack.pop()
                    except IndexError as e:
                        raise NuskellObjectError(
                            "Too many closing brackets in secondary structure")
                    pt[i] = j
                    pt[j] = i
                elif (char == '+'):
                    pt[i] = '+'
                elif (char not in set(chars)):
                    raise NuskellObjectError(
                        "Unexpected character in sequence: '" + char + "'")
        
            if stack != []:
                raise NuskellObjectError(
                    "Too many opening brackets in secondary structure")
            return pt

        fh.write("(* File autogenerated by nuskell. ")

        if ts:
            fh.write("\n - Translation Scheme: {}".format(ts))
        if crn:
            fh.write("\n - Input CRN: \n")
            for rxn in crn:
                assert len(rxn) == 3
                if len(rxn[2]) == 2:
                    fh.write("    {} <=> {}\n".format(
                        ' + '.join(rxn[0]), ' + '.join(rxn[1])))
                else:
                    fh.write("    {} -> {}\n".format(
                        ' + '.join(rxn[0]), ' + '.join(rxn[1])))

        fh.write("*)\n\n".format(crn))

        fh.write("def Fuel = 20\n")
        fh.write("def Signal = 5\n\n")

        first = True
        for cplx in sorted(self._testtube.complexes, key=lambda x: x.name):

            if first:
                fh.write('( ')
                first = False
            else:
                fh.write('| ')

            if cplx.name in signals:
                fh.write("Signal * ")
            else:
                fh.write("constant Fuel * ")

            name = cplx.name
            sequ = cplx.sequence
            stru = cplx.structure

            ptab = pair_table(stru)

            dnaexpr = [[]]
            pos = 0
            for e, d in enumerate(ptab):
                if d == '+':
                    flag = 'top' if flag == 'bound' else flag
                    expr = 'cut'

                elif d == -1:
                    toe = '^' if sequ[e].name[0] == 't' else ''
                    if sequ[e].name[-1] == '*':
                        flag = 'bottom'
                        expr = sequ[e].name[:-1] + toe + '*'
                    elif sequ[e].name[0] == 'h':
                        flag = 'top'
                        expr = '_'
                    else:
                        flag = 'top'
                        expr = sequ[e].name + toe

                elif d > e:  # '('
                    flag = 'bound'
                    toe = '^' if sequ[e].name[0] == 't' else ''
                    if sequ[e].name[-1] == '*':
                        expr = sequ[e].name[:-1] + toe + '*'
                    elif sequ[e].name[0] == 'h':
                        raise NuskellObjectError(
                            'Unexpected bound history domain.')
                    else:
                        expr = sequ[e].name + toe

                    dnaexpr.append([])
                    pos += 1

                elif d < e:  # ')'
                    flag = 'bottom'
                    expr = None
                    pos -= 1
                    if pos < 0:
                        raise NuskellObjectError('too many closing base-pairs')
                    continue
                else:
                    raise NuskellObjectError('strange case:', e, d)

                if dnaexpr[pos] == []:
                    dnaexpr[pos] = [[flag, expr]]
                else:
                    dnaexpr[pos].append([flag, expr])

            # decode dnaexpr
            dnaflat = []
            for d in dnaexpr:
                for dd in d:
                    dnaflat.append(dd)

            # PRINT TO FILE
            close = None
            for e, d in enumerate(dnaflat):
                if d[1] == 'cut':
                    fh.write(close)
                    close = None
                    if e == len(dnaflat) - 1:
                        continue
                    if d[0] == 'bottom':
                        fh.write('::')
                    else:
                        fh.write(':')
                    continue

                if d[0] == 'bottom':
                    if close is None:
                        fh.write('{')
                        close = '}'
                    elif close == ']' or close == '>':
                        fh.write('{}{'.format(close))
                        close = '}'

                if d[0] == 'bound':
                    if close is None:
                        fh.write('[')
                        close = ']'
                    elif close == '}' or close == '>':
                        fh.write('{}['.format(close))
                        close = ']'

                if d[0] == 'top':
                    if close is None:
                        fh.write('<')
                        close = '>'
                    elif close == '}' or close == ']':
                        fh.write('{}<'.format(close))
                        close = '>'
                fh.write(" {} ".format(d[1]))
            if close:
                fh.write("{} (* {} *)\n".format(close, name))
            else:
                fh.write(" (* {} *)\n".format(name))

        fh.write(")\n")

