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
    # This script translates reactions with two reactants and two products into a 
    # DSD network. The scheme has been originally proposed by 
    #   - Lakin et al. (2011) "..."
    #
    # As this script is part of the Nuskell tutorial, try it directly, for example:
    #
    # echo "A + B -> C + D" | nuskell --ts tutorial_1.ts
    #
    # -----------------------------------------------------------------------------
    
    # Start with defining a global short domain:
    global toehold = short();

    # Write a class to define domains and structure of formal species:
    class formal(s) = "h t x" | ". . ."
      where {
        h = long();
        t = toehold;
        x = long();
      }
    ;

    # Write a class to define domains and structure of fuel species:
    class formal(s) = "h t x" | ". . ."
      where {
        h = long();
        t = toehold;
        x = long();
      }
    ;


    #	Write a module that applies the gate :
    module rxn(r) =
      if len(r.reactants) == 0 then
        infty(i) + infty(l) + infty(t) + sum(map(infty,b))
        where {
          i = signal();
          [l, t, b] = maingate([i], r.products)}

      else
        infty(l) + infty(t) + sum(map(infty,b))
          where 
            [l, t, b] = maingate(r.reactants, r.products);

    # Finally, write a module that applies ``rxn`` to the crn.
    module main(crn) = sum(map(rxn, crn))
      where crn = irrev_reactions(crn)
    ;

Here is a VisualDSD image of the compiled DSD system.

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

