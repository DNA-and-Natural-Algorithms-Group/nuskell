.. Nuskell documentation master file, created by
   sphinx-quickstart on Mon May  9 21:46:13 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Nuskell Documentation
============================

--------------------------
Introduction
--------------------------

``Nuskell`` integrates biophysical modelling of nucleic acids with rigor of
modern compilers to design and characterize *toehold-mediated DNA strand
displacement* networks. It translates formal chemical reaction networks (CRNs)
into nucleic acid chemistry, and verifies wether a given nucleic acid network
implements a particular dynamic behavior. 

Currently, ``Nuskell`` supports three different I/O formats for domain-level
specifications of DSD networks: **\*.dom** (Nuskell format), **\*.pil** (Pepper
internal language format), **\*.dna** (VisualDSD format).

.. code-block:: bash

  nuskell [options] < formal_crn.in > domain-specification.out


The ``Nuskell`` library functions enable researchers to make custom compilers
for alternative design approaches. ``Nuskell`` provides (i) a number of
different CRN-to-DSD translation schemes, (ii) functions that test different
CRN equivalence notions, as well as (iii) reaction enumeration and CRN-to-ODE
simulations based on emperical (sequence-independent) DNA folding parameters.

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

  v = verify(fcrn, vcrn, notion = 'bisimulation')
  
  if v :
    print("Input CRN and TestTube-Species are bisimulation equivalent.")
  else :
    print("Input CRN and TestTube-Species are not bisimulation equivalent.")


--------------------------
The Nuskell Language
--------------------------

**Translation schemes** provide instructions for Nuskell to translate a
chemical reaction network (CRN) into a domain-level strand displacement (DSD)
system.  The Nuskell programming language is the interface to add a new
translation scheme to the library.

This section describes the philosophy and the ``syntax`` of translation
schemes. 


.. toctree::

  tlstutorial

--------------------------
Compiler API
--------------------------

Reference manual for the ``Nuskell`` compiler library. **Quickstart** introduces
high-level concepts that are supported by library. Base-level functions for 
advanced development are described in the **API Reference** section.

Quickstart
-----------------------

*	Translate a formal CRN into a TestTube object:
  
  .. code-block:: python

    from nuskell import translate

    # Test a simple bimolecular reaction
    crn = 'A + B -> C + D'
    tls = 'soloveichik2010.ts'

    solution, _ = translate(foo, tls)


API Reference
-----------------------

.. toctree::
  :maxdepth: 2

  compiler
  parser
  interpreter
  verifier

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
are a few rules that developers have to respect before submitting a pull request.
Please make sure that your editor or IDE supports these settings.

 * **Use two-whitespace indents.** 

   Tab-characters are strictly forbidden, as they default to four-whitespace
   indents. In fact, the only white-space characters that should be used are 
   ' ' and '\\n'.

 * **Do not exceed 80 characters per line** 

 * **Format docstrings according to [Google Docstring Guidelines]**

 * other than that, ... try to follow the guides for coding as much as possible,
   e.g. Google, Hitchhiker's, ...
   

The main repository of ``Nuskell`` can be found at ``GitHub``:

https://github.com/bad-ants-fleet/nuskell


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

