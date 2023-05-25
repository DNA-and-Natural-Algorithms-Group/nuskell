DSD enumeration
===============
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

.. _Badelt et al. (2020): http://dna.caltech.edu/DNAresearch_publications.html#Peppercorn
.. _peppercornenumerator: https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator
