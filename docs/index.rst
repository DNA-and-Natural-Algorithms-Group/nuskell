.. Nuskell documentation master file, created by
   sphinx-quickstart on Mon May  9 21:46:13 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Nuskell Documentation
=====================

------------
Introduction
------------

``Nuskell`` integrates biophysical modelling of nucleic acids with rigor of
modern compilers to design and characterize *toehold-mediated DNA strand
displacement* networks. It translates formal chemical reaction networks (CRNs)
into nucleic acid chemistry, and verifies wether a given nucleic acid network
implements a particular dynamic behavior. 

Quickstart
-----------------------
The following command translates a formal CRN into a DSD system using the
built-in translation scheme: ``qian2011_3D.ts``, verifies whether the
implementation is correct using bisimulation and prints the domain-level
specification in the ``*.pil`` file format.

.. code-block:: bash

  nuskell --ts qian2011_3D.ts --verify bisimulation --pilfile < formal_crn.in 


The ``Nuskell`` library functions enable researchers to make custom compilers.
``Nuskell`` provides (i) a number of different CRN-to-DSD translation schemes,
(ii) functions that test different CRN equivalence notions, as well as (iii)
reaction enumeration and CRN-to-ODE simulations based on emperical
(sequence-independent) DNA folding parameters.

.. code-block:: python

  from nuskell import translate, verify

  testtube = translate('A+B->C', scheme = 'soloveichik2010.ts')

  # Get the enumerated CRN
  testtube.enumerate_reactions()

  # Interpret the enumerated CRN, i.e. replace history species
  interpretation = testtube.interpret_species(['A','B','C'], prune=True)

  # Formulate reversible reactions as two irreversible reactions.
  fcrn = [['A','B'],['C']]
  vcrn = []
  for r in testtube.reactions:
    rxn = [map(str,r.reactants), map(str,r.products)]
    vcrn.append(rxn)

  v = verify(fcrn, vcrn, fs, method = 'bisimulation')
  
  if v :
    print("Input CRN and TestTube-Species are bisimulation equivalent.")
  else :
    print("Input CRN and TestTube-Species are not bisimulation equivalent.")


-----------------------------
Nuskell programming language
-----------------------------

**Language** Translation schemes provide instructions for ``Nuskell`` to translate a
chemical reaction network (CRN) into a domain-level strand displacement (DSD)
system.  The ``Nuskell`` programming language is the interface to formulate DSD
design principles as algorithms and to add a new translation scheme to the
library.

**DSD requirements** All CRN-to-DSD translation schemes need to compile to a
particular format of DSD systems. First, all involved species and complexes
need to be free of pseudo-knots. Second, we distinguish *signal species* and *fuel
species*.  Signal species are at low concentrations and they present the
information (input/output) unit.  Fuel species are at high (ideally constant)
concentrations and they mediate the information transfer by consuming and/or
releasing signal species.  After compilation, every species in the formal CRN
corresponds to one signal species. Thus, all signal species must have the same
domain-level constitution and structure, but they need to be independent of
each other. A signal species may be a complex composed of multiple molecules.

This section describes the ``syntax`` of translation schemes.  Once the user is
familiar with such simple schemes, the more general concepts of the language,
such as writing macros to support $n$-arity reactions, are straightforward to
implement.

**History domains** History domains are common in many translation
schemes. A history domain is considered to be an inert domain of a signal
species, but it is unique to the reaction that has produced the signal species.
Hence, multiple species that differ only by their history domains map to the
same formal species. In the translation scheme language, a history domain is a
wildcard: ``?``. Together with the remainder of the molecule, a species with a
wildcard forms a regular-expression, matching every other species in the system
that differs only by a single domain instead of ``?``. 

.. In order to produce a minimal domain-level system specification, ``Nuskell``
.. automatically removes signal species that are specified using wildcard domains
.. after domain-level enumeration. In particular, if there exists a species
.. matching the regular expression, then the species with the wildcard domain and
.. every enumerated reaction emerging from that species is be removed from the
.. system, otherwise, the wildcard domain is replaced by a regular long domain.

.. toctree::

  tlstutorial

--------------------------
Compiler API
--------------------------

Reference manual for the ``Nuskell`` compiler library. **Quickstart** introduces
high-level concepts that are supported by library. Base-level functions for 
advanced development are described in the **API Reference** section.


API Reference
-----------------------

.. toctree::
  :maxdepth: 2

  compiler
  interpreter
  verifier

..  enumeration
..  objects
..  parser

.. Nuskell Objects
.. -----------------------
.. 
.. Input Parsers
.. -----------------------
.. 
.. CRN-to-DSD Translation
.. -----------------------
.. 
.. CRN Equivalence
.. -----------------------
.. 
.. Output Formats
.. -----------------------


Developer Guidelines
============================

In order to ensure sustainability of the ``Nuskell`` compiler package, there
are a few rules for developers before submitting a pull request.

 * **Use two-whitespace indents.** 

   Don not use Tab-characters, the only white-space characters that should be
   used are ' ' and '\\n'.

 * **Do not exceed 80 characters per line** 

 * **Use Google docstring format** `[Google Dochstring Guidelines] <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_ 


The ``Nuskell`` `repository`_ can be found on ``GitHub``: 

https://github.com/bad-ants-fleet/nuskell

.. _repository: https://github.com/bad-ants-fleet/nuskell


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

