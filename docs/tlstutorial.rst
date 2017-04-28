The Nuskell programming language
=================================

**Translation schemes** provide instructions for Nuskell to translate a
chemical reaction network (CRN) into a domain-level strand displacement (DSD)
system.  The Nuskell programming language is the interface to add a new
translation scheme to the library.

This section describes the philosophy and the ``syntax`` of translation
schemes. 

.. Nuskell is a functional programming language, 

-----------------------------------
Tutorial script 1 - Fist Steps
-----------------------------------

This basic script translates a particular input CRN into a domain-level strand
displacement system. First, it defines formal and fuel complexes, second, it
defines gates, i.e. instructions to design fuel complexes specific for an
reaction of formal species. Last, the main function applies the gate
construction to every reaction in the CRN and unifies all fuel strands.

.. code-block:: none

  # -----------------------------------------------------------------------------
  # Translate formal reactions with two reactants and two products.
  # Lakin et. al (2012) "Abstractions for DNA circuit design." [Figure 5]
  # -----------------------------------------------------------------------------
  
  # Define a global short toehold domain
  global toehold = short();
  
  # Write a class to define domains and structure of signal species
  # ? is a wildcard for a history domain.
  class formal(s) = "? t f" | ". . ."
    where { t = toehold ; f = long() };
  
  # Write a class to produce fuel complexes for bimolecular reactions
  class bimol_fuels(r, p) = 
    [ "a t i + b t k + ch t c + dh t d + t* dh* t* ch* t* b* t* a* t*" 
    | "( ( . + ( ( . + (  ( . + (  ( . + )   )  )   )  )  )  )  )  . ",
      "a t i" | " . . . ", "t ch t dh t" | ". . . . ." ]
    where {
      a = r[0].f; 
      b = r[1].f;
      c = p[0].f; ch = long();
      d = p[1].f; dh = long();
      i = long(); k = long();
      t = toehold };
  
  # Write a module that applies the fuel production to every reaction
  module rxn(r) = sum(map(infty, fuels))
    where fuels = 
      if len(r.reactants) != 2 or len(r.products) != 2 then
        abort('Reaction type not implemented')
      else 
        bimol_fuels(r.reactants, r.products);
  
  # Write the module *main* that applies *rxn* to the crn.
  module main(crn) = sum(map(rxn, crn)) 
    where crn = irrev_reactions(crn);


-----------------------------------
Tutorial script 2 - Generalization
-----------------------------------

-----------------------------------
Tutorial script 3 - Optimization
-----------------------------------

---------------
Nuskell Syntax
---------------

* Keywords:

  **global**
  
  **class**
  
  **function**
  
  **module**
  
  **macro**

* Expressions:

**where { ...; ... }**
the keyword where expects a list of assignments. The format is
either ... where { t=func(t); l=func2(l); ... } or where t=func(t)

**if ... then ... (elif ... then ...) else ...**

----------------------
Built-In Functions
----------------------

**tail(x)**

builtin base_level functions:
tail, complement, infty, unique, flip, rev_reactions, irrev_reactions

builtin sematics:
if, or, and, dna, where, list, uminus

trailer:
apply, index, attribute

builtin header functions:
sum, tail, range, reverse, rxn_degree, unirxn, birxn, map, map2

