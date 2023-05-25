CRN equivalence
===============
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

.. _Shin et al. (2019): http://dna.caltech.edu/DNAresearch_publications.html#PathwayDecomposition
.. _Johnson et al. (2019): http://dna.caltech.edu/DNAresearch_publications.html#CRN-Bisimulation
