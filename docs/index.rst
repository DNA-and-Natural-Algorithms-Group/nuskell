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

Please also consider citing the dependencies of nuskell:
 * Peppercornenumerator: `Badelt et al. (2020)`_
 * Pathway decomposition equivalence: `Shin et al. (2019)`_
 * CRN bisimulation: `Johnson et al. (2019)`_

----------------------
CRN-to-DSD translation
----------------------

.. toctree::
   intro
   tlstutorial
   analysis
   peppercorn
   crn_equivalence

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

Please document code using the `[Google Docstring Guidelines] <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_.
If possible, provide unittests that demonstrate the functionalty of any new code contribution.

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
