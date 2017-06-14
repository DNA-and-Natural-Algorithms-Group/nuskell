#
# Lakin, Youssef, Cardelli, Phillips "Abstractions for DNA circuit design."
# J. R. Soc. Interface 9 (68) (2012) 470-486 (2012)
#
# Note: * Implements Figure 5: {A + B -> C + D}
#       * uses more fuel species (to reverse the release of all reactants)
#       * generalized on the DNA level, for all but {X->}
#       * generalized on the CRN level, for {X->} as {X->f}
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

macro pmac(s) = [ "i t j +" 
                | "( ( . +", 
                  "t* i*" 
                | ")  )",  
                  "t i"   
                | ". ." ]
  where {
    i = long();
    t = s.t;
    j = s.f 
  };

macro rmac(s) = 
  [ "a t i + " 
  | "( ( . + ", 
    "t* a*" 
  | ")  )",
    "a t i" | ". . ."]
  where {
    a = s.f ;
    t = s.t ;
    i = long() 
  };

class lakin2012_variant(r, p) 
  = [ "k + x + y l t* " 
    | "~ + ~ + ~ ~ .  ",
      "z t" 
    | "~ ."] + m
  where {
    t = toehold ; 
    [k, l, m] = flip(map(rmac, r), 3);
    l = reverse(l);
    [x, y, z] = flip(map(pmac, p), 3);
    y = reverse(y) 
  };

module reaction(rxn) = sum(map(infty, lakin2012_variant(r, p)) + extra)
  where {
    extra = if len(rxn.products) >= 1 then [] else [infty(gfuel)];
    r = rxn.reactants;
    p = if len(rxn.products) == 0 then 
          [gfuel]
        else
          rxn.products
   };

# Write the module *main* that applies *rxn* to the crn.
module main(crn) = sum(map(reaction, crn)) 
  where crn = irrev_reactions(crn);

