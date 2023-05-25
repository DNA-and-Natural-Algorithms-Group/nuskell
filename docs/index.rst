.. Nuskell documentation master file, created by
   sphinx-quickstart on Thu May 25 15:03:21 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Nuskell compiler documentation
==============================

.. note::
    The source for the Nuskell project is available on GitHub: 

    https://github.com/DNA-and-Natural-Algorithms-Group/nuskell

.. include:: ../README.md
    :parser: myst_parser.sphinx_

The Nuskell library is meant to help researchers building custom **DSD**
compilers.  Nuskell provides (i) a number of different **CRN-to-DSD**
translation schemes, (ii) functions that test different **CRN** equivalence
notions, as well as (iii) **DSD** reaction enumeration and **CRN-to-ODE**
simulations based on empirical (sequence-independent) DNA folding parameters.

----------------------
CRN-to-DSD translation
----------------------

DSD system requirements
-----------------------
Automated compilation of **DSD** system requires **DSD** systems to follow a
particular format.  First, all involved species and complexes need to be free
of pseudo-knots. Second, we distinguish **signal species** and **fuel
species**.  Signal species are at low concentrations and they present the
information (input/output) unit.  Fuel species are at high (ideally constant)
concentrations and they mediate the information transfer by consuming and/or
releasing signal species.  After compilation, every species in the formal CRN
corresponds to one signal species. Thus, all signal species must have the same
domain-level constitution and structure, but they need to be independent of
each other. A signal species may be a complex composed of multiple molecules.

**History domains** are common in many translation schemes. A history domain is
considered to be an inert domain of a signal species, but it is unique to the
reaction that has produced the signal species.  Hence, multiple species that
differ only by their history domains map to the same formal species. In the
translation scheme language, a history domain is a **wildcard**: ``?``. Together
with the remainder of the molecule, a species with a wildcard forms a
regular-expression, matching every other species in the system that differs only
by a single domain instead of ``?``. 

.. In order to produce a minimal domain-level system specification, Nuskell
.. automatically removes signal species that are specified using wildcard domains
.. after domain-level enumeration. In particular, if there exists a species
.. matching the regular expression, then the species with the wildcard domain and
.. every enumerated reaction emerging from that species is be removed from the
.. system, otherwise, the wildcard domain is replaced by a regular long domain.

The Nuskell programming language
----------------------------------
The ``Nuskell programming language`` is used to write ``translation schemes``,
i.e. design algorithms which can be interpreted by the Nuskell compiler. 
Translation schemes translate a **CRN** into a **DSD** system.  A library of existing
schemes can be found in the official Nuskell `repository`_, the links below
point to a tutorial to write your own translation scheme:

.. toctree::

  tlstutorial

Analysis of implementation networks
-----------------------------------
Besides the top-down interface of Nuskell to translate CRNs to DSD systems,
there also exists a modular, bottom-up interface where users can analyze,
simulate and verify handcrafted or alternatively designed DSD systems.

.. code-block:: bash

  nuskell --readpil zhang2007_catalyst.pil --verify bisimulation < formal_crn.in 

The option ``--readpil <file>`` tells Nuskell to load domain-level
specifications from a text file, as opposed to automated design via translation
schemes. The input format is a `variation` of the pepper internal language
(**PIL**) kernel notation which allows the specification of ``constant`` or
``initial`` concentrations in ``M``, ``mM``, ``uM``, ``nM``, ``pM``.

.. code-block:: none

  # Use '#' for comments.

  # Domains
  length d1  = 10
  length d2a =  6
 
  # Complexes               # Concentratios
  C = d4 d5                 @initial 2 nM
  OB = d1 d2a d2b d2c       @constant 100 nM

The concentration specification (e.g. ``@initial 10 nM``) is **optional, but
relevant** for both verification and simulation of DSD systems. Nuskell's
verification has to be provided with the information of which species correspond
to signal and fuel species. 

A complex with a **corresponding name in the formal CRN**, is **always**
interpreted as a **signal species**, independent of whether or not `constant` or
`initial` concentrations have been specified. Species that are not present in
formal CRN default to **fuel species** if: **(i)** they have no concentration specified,
or **(ii)** their concentration is higher than ``0``.  Variant (i) allows a compact
DSD system specification, which is equivalent to the format of a ``PIL`` file when
using option ``--pilfile``, and compatible with input for the
peppercornenumerator. 
Variant (ii) enables us to define named **intermediate complexes** as those which
are explicitly **initially not present**, i.e are followed by the ``@initial 0
nM`` tag.  Note, do **not** use ``@constant 0 nM`` to specify an intermediate
species, as the behavior of Nuskell is currently undefined and might change in
future versions.

The following ``PIL`` file shows a complete DSD system specification, including
initial concentrations for signal species, formal species and all enumerated
intermediate species:

.. code-block:: none

  #
  # Zhang, Turberfield, Yurke, Winfree (2007) 
  # "Engineering Entropy-Driven Reactions and Networks Catalyzed by DNA"
  #
  # A DSD implementation of the catalyst reaction (Figure 1A + 1D)
  # Note: Domain 2 is actually contains two toeholds (2a, 2b)
  #
  # CRN:
  #   C + S -> C + OB
  #   OB -> ROX
  #
  # verify:
  #   echo "C + S -> C + OB; OB -> ROX" | nuskell --readpil zhang2007_catalyst.pil --verify pathway bisimulation
  #     => not pathway equivalent
  #     => bisimulation equivalent
  #   echo "C -> C + OB; OB -> ROX" | nuskell --readpil zhang2007_catalyst.pil --verify pathway bisimulation
  #     => not pathway equivalent
  #     => bisimulation equivalent
  #
  # Coded by Stefan Badelt (badelt@caltech.edu)
  
  # Domains
  length d1  = 10
  length d2a =  6
  length d2b =  6
  length d2c = 12
  length d3  =  4
  length d4  = 16 
  length d5  =  6
  length d6  = 16
  
  # Species
  C = d4 d5                 @initial 2 nM     # defaults to fuel
  OB = d1 d2a d2b d2c       @initial 0 nM     # defaults to intermediate
  ROX = d1 d2a              @initial 0 nM     # defaults to intermediate
  S = d1 d2a( d2b( d2c( + d6 d3( d4( + d5* ) ) ) ) )  @initial 100 nM # defaults to fuel
  F = d2a d2b d2c d3 d4     @initial 100 nM   # defaults to fuel
  OR = d1( d2a( + d2b* ) )  @initial 100 nM   # defaults to fuel
  SB = d6 d3 d4             @initial 0 nM     # defaults to intermediate


.. _Background :

----------
Background
----------

DSD enumeration
---------------
The domain-level representation provides a more coarse-grained perspective on
nucleic acid folding than the single-nucleotide-level. At the nucleotide-level
every step is a base pair opening or closing reaction and the corresponding rate
can be calculated from the free energy change of a reaction. On the
domain-level, we consider a more diverse set of reactions in order to compensate
for the fine-grained details that can happen on the sequence level.  Nuskell
uses a domain-level reaction enumeration library [`peppercornenumerator`_], to
predict desired and, potentially, undesired reactions emerging from previously
compiled signal and fuel species. 

The general types of reactions are spontaneous **binding and unbinding of domains**,
**3-way branch migration**, **4-way branch migration** and **remote toehold
branch-migration**. Peppercorn's enumeration semantics are
justified based on the assumption that the DSD system is operated at
sufficiently low concentrations, such that unimolecular reactions always go to
completion before the next bimolecular interaction takes place.
Under the assumptions of low concentrations, a **condensed** CRN can be
calculated, with reactions that indicate just the eventual results after all
unimolecular reactions complete, and with rate constants systematically derived
from the detailed reaction network rate constants. 

There are **many** options available to adjust the semantics of reaction
enumeration, they are described in detail in [`Badelt et al. (2020)`_].

**More:**
 * Peppercornenumerator: `Badelt et al. (2020)`_

CRN equivalence 
---------------
The most fundamental requirement towards compilation of large scale DSD systems,
is verification. Every formal reaction is translated into multiple
implementation reactions. Thus, there are many possibilities to introduce
`bugs`, i.e. unwanted side reactions that alter the implemented algorithm. 
Nuskell supports currently variants of two case-by-case verification strategies
that compare formal CRNs with their implementations.  As intended, our approach
does not verify the general correctness of a particular scheme, but the
correctness of a particular implementation.

**Pathway decomposition equivalence.** 
The core idea is to represent each implementation trajectory as a
combination of independent pathways of reactions between formal species.
Pathway decomposition yields a set of pathways which are
indivisible (or `prime`) and are called the `formal basis` of a CRN.
The formal basis is unique for any valid implementation. Any two CRNs are said
to be equivalent if they have the same formal basis. Conveniently, a CRN
without intermediate species has itself as the formal basis, but it is worth
pointing out that this equivalence relation allows for the comparison of one
implementation with another implementation.  

**CRN bisimulation verification.** A CRN bisimulation is an
interpretation of the implementation CRN, where every implementation
species is mapped to a multiset of formal species. This often yields so-called
`trivial` reactions, where reactants and products do not change according
to the interpretation.  An
interpretation is only a bisimulation, if three conditions are fulfilled: **(i)
atomic condition** -- for every formal species there exists an
implementation species that interprets to it, **(ii) delimiting condition**
-- any reaction in the implementation is either trivial or a valid formal
reaction, and **(iii) permissive condition** -- for any initial condition in
the implementation CRN, the set of possible next non-trivial reactions is
exactly the same as it would be in the formal CRN.  CRNs are said to be
bisimulation equivalent, if the translation can be interpreted as an
implementation of that formal CRN. 

**More:**
 * CRN bisimulation: `Johnson et al. (2019)`_
 * Pathway decomposition equivalence: `Shin et al. (2019)`_

Simulations of DSD systems
--------------------------
Translate enumerated CRNs into ODEs

.. _API Reference:

-----------
Nuskell API
-----------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   framework
   ioutils
   objects
   crnutils
   dsdcompiler.objects
   dsdcompiler.compiler
   dsdcompiler.interpreter
   dsdenumerator
   crnverifier

.. _Developer Guidelines:

Developer Guidelines
====================

In order to ensure sustainability of the Nuskell compiler package, there
are a few rules for developers before submitting a pull request.

 * **Use Google docstring format** `[Google Docstring Guidelines] <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_ 

**Cite:**
 * Nuskell: `Badelt et al. (2017)`_
 * Peppercornenumerator: `Badelt et al. (2020)`_
 * Pathway decomposition equivalence: `Shin et al. (2019)`_
 * CRN bisimulation: `Johnson et al. (2019)`_


.. _repository: https://github.com/DNA-and-Natural-Algorithms-Group/nuskell
.. _peppercornenumerator: https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator
.. _crnsimulator: https://github.com/bad-ants-fleet/crnsimulator

.. _Badelt et al. (2017): http://dna.caltech.edu/DNAresearch_publications.html#NuskellCompiler
.. _Badelt et al. (2020): http://dna.caltech.edu/DNAresearch_publications.html#Peppercorn
.. _Shin et al. (2019): http://dna.caltech.edu/DNAresearch_publications.html#PathwayDecomposition
.. _Johnson et al. (2019): http://dna.caltech.edu/DNAresearch_publications.html#CRN-Bisimulation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
