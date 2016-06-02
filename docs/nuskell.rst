The nuskell programming language
=================

The translation scheme of nuskell scripts can contain the following functions
Code has to be written as a list of functions, with ';' delimiter.

------------------
**Keywords**
------------------

**global**

**class**

**function**

**module**
module main(crn) = sum(map(rxn, crn))
    where
        crn = irrev_reactions(crn)

**macro**

------------------
**Expressions**
------------------

**where { ...; ... }**
the keyword where expects a list of assignments. The format is
either ... where { t=func(t); l=func2(l); ... } or where t=func(t)

**if ... then ... (elif ... then ...) else ...**

------------------
**Built-in functions**
------------------

**tail(x)**



builtin base_level functions:
tail, complement, infty, unique, flip, rev_reactions, irrev_reactions

builtin sematics:
if, or, and, dna, where, list, uminus

trailer:
apply, index, attribute

builtin header functions:
sum, tail, range, reverse, rxn_degree, unirxn, birxn, map, map2



------------------
**parsing grammar**
------------------
.. automodule:: nuskell.parser.ts_parser
  :members:


