#
#  nuskell/objects.py
#  NuskellCompilerProject
#
from __future__ import absolute_import, division, print_function
from builtins import map, filter

import logging
log = logging.getLogger(__name__)

import networkx as nx
from collections import Counter
from dsdobjects import clear_memory
from dsdobjects import DL_Domain, DSD_Complex, DSD_Reaction, DSD_Macrostate
from dsdobjects import DSDObjectsError, DSDDuplicationError
from dsdobjects.utils import natural_sort
from dsdobjects.prototypes import Reaction as NuskellReaction
from dsdobjects.prototypes import Complex as NuskellComplex
import dsdobjects.objectio as oio
from peppercornenumerator import Enumerator
from peppercornenumerator.objects import PepperComplex, PepperDomain
from peppercornenumerator.utils import PeppercornUsageError

NuskellReaction.RTYPES.add('formal')

from nuskell import __version__
from nuskell.crnutils import genCRN, genCON

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

    def __new__(cls, name='', prefix='d', dtype=None, length=None):
        # The new method returns the present instance of an object, if it exists
        self = DL_Domain.__new__(cls)
        if name == '':
            if prefix == '':
                raise NuskellObjectError('NuskellDomain prefix must not be empty!')
            elif prefix[-1].isdigit():
                raise NuskellObjectError('NuskellDomain prefix must not end with a digit!')
            name = prefix + str(NuskellDomain.ID)
            NuskellDomain.ID += 1
        try:
            super(NuskellDomain, self).__init__(name, dtype, length)
        except DSDDuplicationError as e:
            other = e.existing
            if dtype and (other.dtype != dtype) :
                raise DSDObjectsError('Conflicting dtype assignments for {}: "{}" vs. "{}"'.format(
                    name, dtype, other.dtype))
            elif length and (other.length != length) :
                raise DSDObjectsError('Conflicting length assignments for {}: "{}" vs. "{}"'.format(
                    name, length, other.length))
            return e.existing

        self.nucleotides = None
        return self

    def __init__(self, name='', prefix='d', dtype=None, length=None):
        # Remove default initialziation to get __new__ to work
        pass

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

class NuskellMacrostate(DSD_Macrostate):
    pass

class TestTube(object):
    r"""A reaction network of nucleic acid complexes.

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
      only by a single domain instead of '?'.  If there exists a species matching
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

    def __init__(self, complexes = None, reactions = None):
        self._RG = nx.MultiDiGraph()

        if complexes and isinstance(complexes, dict):
            for name, data in complexes.items():
                assert isinstance(data, list) and len(data) == 4, 'wrong format for complexes'
                if data[1] is not None :
                    assert isinstance(data[1], float), 'concentration must be of type float.'
                if data[2] is not None :
                    assert isinstance(data[2], bool), '"constant" must be of type bool.'
                if data[3] is not None :
                    assert isinstance(data[3], str), '"ctype" must be of type string.'

                cplx = data[0]
                conc = data[1]
                const = data[2]
                ct = data[3]
                if isinstance(cplx, NuskellComplex):
                    self._RG.add_node(cplx, concentration = conc, constant = const, ctype = ct)
                else:
                    raise NotImplementedError('please specify NuskellComplex objects.')

        elif complexes and isinstance(complexes, list):
            for c in complexes:
                assert isinstance(c, NuskellComplex)
                self._RG.add_node(c, concentration = None, ctype = None, constant = None)

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
        return [n for n in self._RG.nodes if isinstance(n, NuskellComplex)]

    @property
    def signal_complexes(self):
        """list: a list of :obj:`NuskellComplex()` objects. """
        return [n for n in self._RG.nodes if \
                isinstance(n, NuskellComplex) and self._RG.nodes[n]['ctype'] == 'signal']

    @property
    def fuel_complexes(self):
        """list: a list of :obj:`NuskellComplex()` objects. """
        return [n for n in self._RG.nodes if \
                isinstance(n, NuskellComplex) and self._RG.nodes[n]['ctype'] == 'fuel']

    @property
    def unspecified_complexes(self):
        """list: a list of :obj:`NuskellComplex()` objects. """
        return [n for n in self._RG.nodes if \
                isinstance(n, NuskellComplex) and self._RG.nodes[n]['ctype'] is None]

    def has_complex(self, cplx):
        return self._RG.has_node(cplx)

    @property
    def reactions(self):
        """list: a list of :obj:`NuskellReaction()` objects. """
        return [n for n in self._RG.nodes() if isinstance(n, NuskellReaction)]

    # TODO: set complex concentrations in complex objects??
    def set_complex_concentration(self, cplx, concentration, constant):
        """
        Args:
          cplx (:obj:`NuskellComplex()`): complex.
          concentration (flt): concentration.
          constant (bool): True if the concentration is kept constant, False for
            initial concentration.
        """
        self._RG.nodes[cplx]['concentration'] = concentration
        self._RG.nodes[cplx]['constant'] = constant

    def get_complex_concentration(self, cplx):
        """flt, bool: First value is concentration, second is True if it is constant,
        False if it is variable."""
        return self._RG.nodes[cplx]['concentration'], self._RG.nodes[cplx]['constant']

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
        return [n for n in self._RG.nodes if isinstance(n, NuskellComplex) and \
                (self._RG.nodes[n]['concentration'] and \
                self._RG.nodes[n]['concentration'] > th) and n.name not in exclude]

    def interpret_species(self, species, prune = True):
        """Get an interpretation dictionary.

        If a :obj:`NuskellComplex()` sequence contains a wildcard, then this
        function will find all matching complexes, and return those as
        interpretation.  Regex-nodes may have at most *one wildcard* per
        complex, a wildcard corresponds to exactly *one unpaired long domain*.

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

            pMx = [list(map(str, x.sequence)), list(map(str, x.structure))]
            pMy = [list(map(str, y.sequence)), list(map(str, y.structure))]
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
                log.warning('No complex found with name of formal species: {}'.format(fs))
                continue
            else:
                assert len(cplxs) == 1, 'Duplicate complex names?'

            cplx = cplxs[0]
            hist = list(filter(lambda x: x[0] == 'h', map(str, cplx.sequence)))
            if hist :
                assert len(hist) == 1, 'no support for multiple history domains'
                hist = hist[0]
                matches = get_matching_complexes(cplx, hist)
                if matches:
                    need_to_prune = True
                    for e, m in enumerate(matches, 1):
                        m.name = fs + '_' + str(e) + '_'
                        interpretation[m.name] = Counter([fs])
                        self.ReactionGraph.nodes[m]['ctype'] = 'signal'
                    self.rm_complex(cplx, force = True)
                else:
                    # NOTE: We cannot simply remove the domain, because we would need to
                    # enumerate the network again and remove the domain everywhere! So
                    # unless we enumerate up-front with a history-pruned species, this
                    # gets us into trouble.
                    interpretation[cplx.name] = Counter([fs])
            else:
                interpretation[cplx.name] = Counter([fs])

        if prune and need_to_prune:
            log.debug('Pruning the network.')
            # Get rid of all reactions with history wildcards. Start with a set
            # of produce molecules and see what species emerge from reactions
            # consuming these molecules.
            # Alternative: enumerate again using history-replaced species.
            rxns = self.reactions
            [prev, total] = [set(), 
                    set(list(interpretation.keys()) + list(map(str, self.present_complexes())))]
            log.debug('Prev {}, Total {}'.format(prev, total))
            while prev != total:
                prev = set(list(total)) # force a copy?
                for rxn in rxns:
                    self.rm_reaction(rxn)
                    r = set(map(str, rxn.reactants))
                    p = set(map(str, rxn.products))
                    log.debug('R {}, P {}'.format(r, p))
                    if r.intersection(total) == r: 
                        total |= p
                    log.debug('Total = {}'.format(total))

            # Now add all reactions that are possible from the pruned state.
            [self.add_reaction(x) for x in rxns if set(map(str, x.reactants)).issubset(total)]

            # Now remove all the left-over complexes from the graph.
            all_nodes = set(self.complexes)
            assert set(map(str, all_nodes)).issuperset(total)
            total = set([n for n in self.complexes if str(n) in total])
            remove = all_nodes.difference(total)
            self._RG.remove_nodes_from(remove)

        return interpretation

    def add_complex(self, cplx, conctup = (None, None), ctype = None, sanitycheck = None):
        """Add a complex to the TestTube.

        Args:
          cplx (:obj:`NuskellComplex()`): The complex object.
          (conc, const) (flt, bool): Concentration and True/False for constant or
            initial concentrations.

        Note:
          A new complex resets TestTube.domains and TestTube.strands

        """
        if sanitycheck is not None:
            log.warning('TestTube.add_complex: Using deprecated optional argument: sanitycheck')

        (conc, const) = conctup
        assert isinstance(cplx, NuskellComplex), 'must be a NuskellComplex format'

        if self._RG.has_node(cplx): # Does not check for the name
            if cplx.name not in map(str, self.complexes):
                raise NuskellObjectError('Complex Name does not match canonical form: {}'.format(cplx))
            # NOTE: I am not happy about using add_complex to update concentrations...
            if conc is not None:
                if self._RG.nodes[cplx]['concentration'] is None:
                    self._RG.nodes[cplx]['concentration'] = conc
                else:
                    assert self._RG.nodes[cplx]['concentration'] == conc, \
                            "Conflicting complex concentrations"
            if const is not None:
                if self._RG.nodes[cplx]['constant'] is None:
                    self._RG.nodes[cplx]['constant'] = const
                else:
                    assert self._RG.nodes[cplx]['constant'] == const, \
                            "Conflicting complex concentrations"
            if ctype is not None:
                if self._RG.nodes[cplx]['ctype'] is None:
                    self._RG.nodes[cplx]['ctype'] = ctype
                else:
                    assert self._RG.nodes[cplx]['ctype'] == ctype, "Conflicting complex type"
        else:
            if cplx.canonical_form in map(lambda x: x.canonical_form, self.complexes):
                raise NuskellObjectError('trying to add duplicate complex: {} = {}'.format(
                    cplx.name, cplx.kernel_string))
            else:
                self._RG.add_node(cplx, concentration = conc, constant = const, ctype = ctype)
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
                for (r, c) in list(self._RG.in_edges(cplx)):
                    assert isinstance(r, NuskellReaction)
                    self.rm_reaction(r)
                for (c, r) in list(self._RG.out_edges(cplx)):
                    assert isinstance(r, NuskellReaction)
                    self.rm_reaction(r)
            elif (self._RG.in_edges(cplx) or self._RG.out_edges(cplx)):
                raise NuskellObjectError(
                    "Cannot remove a complex engaged in reactions.")
            self._RG.remove_node(cplx)
            self._domains = None
            self._strands = None

    def add_reaction(self, react, sanitycheck = None):
        """Add a reaction to the TestTube.

        Args:
          react (:obj:`NuskellReaction()`): The *irreversible* reaction to be added.
        """

        if sanitycheck is not None:
            log.warning('TestTube.add_complex: Using deprecated optional argument: sanitycheck')

        assert isinstance(react, NuskellReaction) 

        if self._RG.has_node(react):
            assert self._RG.nodes[react]['rate'] == react.rate
        else:
            if react.canonical_form in map(lambda x: x.canonical_form, self.reactions):
                # This cannot be true, ever!
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

    def enumerate_reactions(self, pargs, prefix = 'e', condensed = True):
        """Enumerate reactions using the *peppercorn* enumerator.
        Args:
          args(:obj:`argparse.ArgumentParser()`, optional): Arguments for *peppercorn*.
          condensed (bool, optional): Udate the reaction graph using *condensed* format.
        """
        PepperComplex.PREFIX = prefix

        if not condensed:
            raise NotImplementedError('detailed reaction networks currently not supported.')

        selfIO = TestTubeIO(self)
        pilfile = selfIO.write_pil()

        # Memory management.
        backupCM = DSD_Complex.MEMORY
        backupRM = DSD_Reaction.MEMORY
        clear_memory()

        cxs = {}
        for can, obj in backupCM.items():
            seq = []
            for dom in obj.sequence:
                if dom == '+':
                    seq.append(dom)
                else:
                    seq.append(PepperDomain(dom.name, dtype = dom.dtype, length = dom.length))
            cxs[obj.name] = PepperComplex(seq, obj.structure, obj.name)

        init_cplxs = [cxs[n.name] for n in self.complexes if n.concentration is None or n.concentration.value != 0]
        
        enum = Enumerator(init_cplxs)
        set_peppercorn_args(enum, pargs)
        enum.enumerate()
        enum.condense()

        out = enum.to_pil(detailed = False, condensed = True)

        clear_memory()
        selfIO.load_pil(out)

        for can, obj in backupCM.items():
            if can not in DSD_Complex.MEMORY:
                DSD_Complex.MEMORY[can] = obj
            else:
                assert DSD_Complex.MEMORY[can] == obj

        for can, obj in backupRM.items():
            if can not in DSD_Reaction.MEMORY:
                DSD_Reaction.MEMORY[can] = obj
            else:
                assert DSD_Reaction.MEMORY[can] == obj

    def __add__(self, other):
        assert isinstance(other, TestTube)
        combined = TestTube()
        combined.ReactionGraph = nx.compose(self.ReactionGraph, other.ReactionGraph)
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

    def write_pil(self, fh = None, unit = 'M', crn = None, fs = None, ts = None):
        """Write the contents of :obj:`TestTube()` into a PIL file -- KERNEL notation).

        Args:
          fh (filehandle): A filehandle that the output is written to or None.
          unit (str, optional): Specify a unit of concentrations (M, mM, uM, nM, pM).
          crn (list[list], optional): a nuskell-style CRN expression
          ts (str, optional): name of the translation scheme

        Example:
          length d1 = 6
          length d2 = 4
          length h3 = 1
          cplx1 = h3 d1( d2( + )) @ initial 10 nM
        """

        out = []
        def output_string(string):
            if fh is None:
                out.append(string)
            else :
                fh.write(string)

        # Write header #
        output_string("# File generated by nuskell-{}\n".format(__version__))
        output_string("#\n")
        if ts:
            output_string("# - Translation Scheme: {}\n".format(ts))
            output_string("#\n")
        if fs:
            output_string("# - Input concentration: \n")
            for ico in genCON(fs):
                output_string('#    {}\n'.format(ico))
            output_string("#\n")
        if crn:
            output_string("# - Input CRN: \n")
            for rxn in genCRN(crn):
                output_string('#    {}\n'.format(rxn))
            output_string("#\n")

        domains = self._testtube.domains
        def adjust_conc(conc, unit):
            units = ['M', 'mM', 'uM', 'nM', 'pM']
            # 0,  3,   6,   9,   12
            assert unit in units
            mult = units.index(unit) * 3
            return conc * (10**mult), unit

        # Print Domains
        output_string("\n# Domain Specifications\n")
        seen = set()
        for d in natural_sort(domains):
            if d.is_complement:
                dom = ~d
            else :
                dom = d
            if dom not in seen:
                output_string("length {:s} = {:d}\n".format(dom.name, dom.length))
                seen.add(dom)


        def print_cplxs(complexlist):
            for cplx in complexlist:
                output_string("{:s} = ".format(cplx.name))
                seq = cplx.sequence
                sst = cplx.structure
                for i in range(len(seq)):
                    if sst[i] == '+':
                        output_string("{:s} ".format(str(sst[i])))
                    elif sst[i] == ')':
                        output_string("{:s} ".format(str(sst[i])))
                    elif sst[i] == '(':
                        output_string("{:s} ".format(str(seq[i]) + str(sst[i])))
                    else:
                        output_string("{:s} ".format(str(seq[i])))

                conc, const = self._testtube.get_complex_concentration(cplx)
                if const is True:
                    output_string(" @constant {} {}".format(*adjust_conc(conc, unit)))
                elif const is False:
                    output_string(" @initial {} {}".format(*adjust_conc(conc, unit)))
                output_string(" \n")

        sc = natural_sort(self._testtube.signal_complexes)
        if len(sc):
            output_string("\n# Signal complexes ({})\n".format(len(sc)))
            print_cplxs(sc)

        fc = natural_sort(self._testtube.fuel_complexes)
        if len(fc):
            output_string("\n# Fuel complexes ({})\n".format(len(fc)))
            print_cplxs(fc)

        oc = natural_sort(self._testtube.unspecified_complexes)
        if len(oc):
            output_string("\n# Other complexes ({})\n".format(len(oc)))
            print_cplxs(oc)


        output_string("\n# Reactions ({})\n".format(len(self._testtube.reactions)))
        for rxn in natural_sort(self._testtube.reactions):
            output_string("reaction {:s}\n".format(rxn.full_string(unit, 's')))

        return ''.join(out)

    def load_pil(self, data, is_file = False):
        """Parses a file written in PIL notation! """
        oio.LogicDomain = NuskellDomain
        oio.Complex = NuskellComplex
        oio.Reaction = NuskellReaction
        oio.Macrostate = NuskellMacrostate

        doms, cplxs, rms, det, con = oio.read_pil(data, is_file)

        for name, rm in rms.items():
            self._testtube.add_complex(rm.canonical_complex) # concentration, ctype?

        for rxn in con:
            del DSD_Reaction.MEMORY[rxn.canonical_form] # avoid duplication warning 
            react = []
            for rm in rxn.reactants:
                react.append(rm.canonical_complex)
            prod = []
            for rm in rxn.products:
                prod.append(rm.canonical_complex)
            self._testtube.add_reaction(
                    NuskellReaction(react, prod, rate = rxn.rate, rtype = rxn.rtype))
        return

    def write_dnafile(self, fh, signals=[], crn=None, ts=None):
        r""" Write a TestTube Object into VisualDSD \*.dna format.

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

        fh.write("(* File generated by nuskell. ")

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

    if hasattr(args, 'no_max_helix'):
        enum.max_helix = not args.no_max_helix
    else:
        enum.max_helix = True

    if hasattr(args, 'k_slow'):
        enum.k_slow = args.k_slow
    else:
        enum.k_slow = 0.0

    if hasattr(args, 'k_fast'):
        enum.k_fast = args.k_fast
    else:
        enum.k_slow = 0.0

    return
