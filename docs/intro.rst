DSD system requirements
=======================
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


