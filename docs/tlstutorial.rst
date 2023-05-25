.. _Nuskell programming language:

The `Nuskell` programming language
==================================


**Translation schemes** are design algorithms that translate a chemical
reaction network (CRN) into a domain-level strand displacement (DSD) system.  
The `Nuskell` programming language is inspired by the functional programming
language *Haskell* and provides DSD specific classes, functions and macros to
generalize translations for arbitrary CRNs.  
This section describes the ``syntax`` of the `Nuskell` programming language in
order to add new translation schemes to the scheme library.
A library of existing schemes can be found in the official Nuskell `repository`_.

------
Syntax
------

Every translation scheme consists of variable assignments in the form of:

.. code-block:: none

  name = value

and functions in the form of:
    
.. code-block:: none

  declarator function(arg, ...) = value ;

The ``where`` keyword allows for more verbose formulations, the character ``#``
is used for (inline) comments:

.. code-block:: none

  # This is a comment.
  declarator function(arg, ...) = result
    where 
      result = value;                     # delimiter ';' closes the function 

  declarator function1(arg, ...) = result1 + result2
    where {
      result1 = value;                    # delimiter ';' separates two assignments
      result2 = function2(value)          # no delimiter ';' before closing '}'
    };

Conditional statements can be written using ``if``, ``then``, ``else``,
``elseif`` keywords.  The operators ``and``, ``or``, ``*``, ``/``, ``+``,
``-``, ``==``, ``!=``, ``<``, ``>``, ``>=``, ``<=`` are supported and logically
equivalent to their implementation in `Python`:

.. code-block:: none

  # Operators are treated as in Python.
  declarator function1(arg1, arg2, arg3, arg4) = result 
    where
      result = 
        if arg1 < arg2 then
          value1
        elseif arg3 and arg4 then
          value2
        else
          value3 ;

Note that **white space** formatting is optional, all statements above can be written
on a single line. However, the use of ``{}`` and ``;`` as delimiters is obligatory.
Every function assignment has to be closed by a ``;``. A ``where`` statement can
be followed by a single assignment, or a list of assignments ``{}`` with a ``;`` as
delimiter.

Built-in functions
------------------
The `Nuskell` language provides a number of built-in functions.

* ``y = short()`` -- returns a new toehold domain.
* ``y = long()``  -- returns a new branch-migration domain.

* ``y = infty(x)`` -- returns a set with one fuel from input complex ``x`` (assigns ``infinite`` concentration)
* ``y = empty``    -- return an empty set of fuels.

* ``y = tail(x)`` -- return a list ``x`` without its first element.
* ``y = flip(x)`` -- return a transposed matrix ``x``. Similar to `Python`'s zip().

* ``y = rev_reactions(crn)``   -- return input ``crn`` such that corresponding irreversible reactions are combined to reversible reactions
* ``y = irrev_reactions(crn)`` -- return input ``crn`` such that reversible reactions are split into two irreversible reactions

* ``y = print(m)`` -- print message ``m``, return nothing(!)
* ``y = abort(m)`` -- exit with message ``m``, return nothing(!)

* ``y = len(x)`` -- returns the length of list ``x``
* ``y = sum(x)`` -- return the sum of list ``x``. Often used to sum over sets of fuels.
* ``y = range(x)`` -- returns a list ``[0 .. x-1]``
* ``y = reverse(x)`` -- returns list ``x`` in reverse 
* ``y = map(f,x)`` -- applies function ``f(x)`` to every element in list ``x``
* ``y = map2(f,y, x)`` -- applies function ``f(y,x)`` to every element in list ``x``

* ``y = birxn(x)`` --
* ``y = unirxn(x)`` --
* ``y = rxn_degree(x,r)`` -- 

Note that even though the functions ``print(m)`` and ``abort(m)`` do not have a
return value, the `Nuskell` language syntax of requires them to be formulated
within an assignment. For example:

.. code-block:: none

  # print() and abort() have no return value.
  declarator function1(arg1, arg2) = result 
    where {
      void = print('Computing result:');
      result = 
        if arg1 < arg2 then
          abort('Error:', arg1, '<', arg2)
        else
          arg2 - arg1;
      void = print('Returning value:', result)
    };

Function declarators
--------------------
Translation schemes can use a variety of function
declarators to indicate function return values:

* ``function`` -- a recursive definition of a function. As an example, some of
  the built-in functions described above are implemented within the `Nuskell`
  language: 

  .. code-block:: none

    function len(x) = 
      if x == [] then 
        0 
      else 
        1 + len(tail(x)) ; # tail(x) returns list x without the first element.

    function sum(x) = if len(x) == 0 then empty elseif len(x) == 1 then x[0] else x[0] + sum(tail(x)) ;
    function map(f, x) = if len(x) == 0 then [] else [f(x[0])] + map(f, tail(x)) ;

* ``class`` -- returns a domain-level complex or a list of domain-level
  complexes. Domain-level complexes are specified as a tuple of sequence and
  structure, for example: ``"a b a*" | "( . )"`` denotes a single strand with
  three domains forming a hairpin loop. By convention ``a*`` denotes a domain
  complementary to ``a``. On the other hand, ``["a b a*" | ". . .", "a b a*" |
  "( . )"]`` is a list of two molecules, which differ in their secondary
  structure, but not in their sequence.

  .. code-block:: none
  
    class get_complexes() = ["a b a*" | "( . )", "a b a*" | ". . ."]
      where {
        a = short(); # Note that "a*" is implicitly assigned, as the complement of "a"
        b = long()
      };
 
    # The required class "formal(s)" must read one argument (a formal species)
    # and returns a single domain-level complex.  The wildcard "?" can be used to 
    # specify history domains, enabling a many to one mapping from singal to
    # formal species.
    class formal(s) = "? t f" | ". . ."
      where { 
        t = short(); 
        f = long() };

    # In most cases, a class will translate a list of reactant signal species "r" 
    # and product signal species "p" into a domain-level complex: 
    class binary_fuel_complexes(r, p) = 
      [ "a t i + b t k + ch t c + dh t d + t* dh* t* ch* t* b* t* a* t*" 
      | "( ( . + ( ( . + (  ( . + (  ( . + )   )  )   )  )  )  )  )  . ",
        "a t i" | " . . . ", "t ch t dh t" | ". . . . ." ]
      where {
        a = r[0].f; # The domain f of the first signal species in list r
        b = r[1].f;
        c = p[0].f; ch = long();
        d = p[1].f; dh = long();
        i = long(); k = long();
        t = short() };
 
* ``macro`` -- has the same return value as ``class``, but is used to denote
  partial domain-level complexes, i.e. a ``class`` can employ ``macros`` to
  generalize translation schemes on the domain-level for arbitrary CRNs.

  .. code-block:: none
  
    # TODO
    class get_fuels() = []
      where {
        [l, p, q] = flip(map(chen2D_O, reverse(prod)), 2);
        [l, p, q] = zip(map(chen2D_O, reverse(prod)));

* ``module`` -- returns a set() of fuel complexes. Fuel complexes are
  domain-level complexes with, ideally, constant concentration. `Nuskell` uses
  the keyword ``infty`` to translate a domain-level complex into a fuel
  complex. The example code below starts with the ``module main()`` which takes
  the input CRN as argument. We will discuss the properties of the ``crn``
  object in detail later (see :ref:`CRN_Object`).

  .. code-block:: none
  
    # The *required* module "main(crn)" translates a CRN into a set of fuel species.
    module main(crn) = sum(map(get_fuels, crn)) 
      where crn = irrev_reactions(crn);

    module get_fuels(reaction) = sum(map(infty, complexes))
      where 
        complexes = get_complexes(reaction.reactants, reaction.products) ;

* ``global`` -- returns a global variable, such as a global domain. 

  .. code-block:: none
  
    global toehold = short() ;


Note: `Nuskell` does *not* enforce the proper usage of ``class``, ``function``,
``module`` and ``macro``, they can be used interchangeably. Only, the
``global`` declarator is specific to the use of global variables.

.. _CRN_Object:

The *crn* object
----------------


-----------------------------------
Tutorial script 1 - Fist Steps
-----------------------------------

There are two required parts: (i) the ``formal`` class defines sequence and
structure of signal complexes, (ii) the ``main`` module produces a set of fuel
species from the input CRN. The following translation scheme translates CRNs in
binary format (two reactants, two products) into a DSD system, and *aborts* the 
translation if it encounters a non-binary reaction.
The classes define signal and fuel complexes, the modules contain instructions
to design fuel complexes specific for a reaction of formal species. 

.. code-block:: none

  # -----------------------------------------------------------------------------
  # Translate formal reactions with two reactants and two products.
  # Lakin et. al (2012) "Abstractions for DNA circuit design." [Figure 5]
  # -----------------------------------------------------------------------------
  #
  # Coded by Stefan Badelt (badelt@caltech.edu)
  
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

----------------------
Built-In Functions
----------------------

.. builtin base_level functions:
.. tail, complement, infty, unique, flip, rev_reactions, irrev_reactions
.. trailer:
.. apply, index, attribute

Built-in functions written in the `Nuskell` programming language:

  .. code-block:: none

    function range(x) = if x == 0 then [] else range(x - 1) + [x - 1] ;

    function sum(x) = if len(x) == 0 then empty elseif len(x) == 1 then x[0] else x[0] + sum(tail(x)) ;

    function len(x) = if x == [] then 0 else 1 + len(tail(x)) ;

    function reverse(x) = if x == [] then [] else reverse(tail(x)) + [x[0]] ;

    function rxn_degree(x, r) = if len(x) == 0 then [] elseif len(x[0].reactants) == r then [x[0]] + rxn_degree(tail(x), r) else rxn_degree(tail(x), r) ;

    function unirxn(x) = if len(x) == 0 then [] elseif len(x[0].reactants) == 1 then [x[0]] + unirxn(tail(x)) else unirxn(tail(x)) ;

    function birxn(x) = if len(x) == 0 then [] elseif len(x[0].reactants) == 2 then [x[0]] + birxn(tail(x)) else birxn(tail(x)) ;

    function map(f, x) = if len(x) == 0 then [] else [f(x[0])] + map(f, tail(x)) ;

    function map2(f, y, x) = if len(x) == 0 then [] else [f(y, x[0])] + map2(f, y, tail(x))

.. _repository: https://github.com/DNA-and-Natural-Algorithms-Group/nuskell
