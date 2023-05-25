Analysis of implementation networks
===================================
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



