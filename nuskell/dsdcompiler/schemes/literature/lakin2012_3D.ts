#
# Lakin, Youssef, Cardelli, Phillips "Abstractions for DNA circuit design."
# J. R. Soc. Interface 9 (68) (2012) 470-486 (2012)
#
# Note: * Implements Figure 5: {A + B -> C + D}
#
#       * no fuel species to reverse the release of second reactant
#
#       * generalized on the CRN level, which makes trimolecular reactions
#         incredibly inefficient.
#
# Coded by Stefan Badelt (badelt@caltech.edu)
#

# Define a global short toehold domain
global toehold = short();

# Write a class to define domains and structure of signal species
class formal(s) = "? t f" | ". . ."
  where { t = toehold ; f = long() };

class inter() = "? t f" | ". . ."
  where { t = toehold ; f = long() };

global gfuel = inter();

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

module prep(r, p) = narity
  where 
    narity = if len(r) == 2 and len(p) == 2 then
                bimol_fuels(r,p)
             else if len(r) > 2 then
                bimol_fuels([r[0], r[1]], [gfuel, i]) + bimol_fuels([gfuel, i], [r[0], r[1]]) + prep([i]+tail(tail(r)), p) + [gfuel]
                    where i = inter()
             else if len(p) > 2 then
                bimol_fuels(r, [j, p[0]]) + prep([j, gfuel], tail(p)) + [gfuel]
                  where j = inter()
             else
                abort('something went terribly wrong')
  ;

module reaction(rxn) = sum(map(infty, prep(r, p)) + extra)
  where {
    extra = if len(rxn.reactants) == 2 and len(rxn.products) == 2 then [] else [infty(gfuel)];
    r = if len(rxn.reactants) == 0 then 
          [gfuel, gfuel]
        else if len(rxn.reactants) == 1 then
          [rxn.reactants[0], gfuel]
        else
          rxn.reactants;
    p = if len(rxn.products) == 0 then 
          [gfuel, gfuel]
        else if len(rxn.products) == 1 then
          [gfuel, rxn.products[0]]
        else
          rxn.products
   };

# Write the module *main* that applies *rxn* to the crn.
module main(crn) = sum(map(reaction, crn)) 
  where crn = irrev_reactions(crn);
