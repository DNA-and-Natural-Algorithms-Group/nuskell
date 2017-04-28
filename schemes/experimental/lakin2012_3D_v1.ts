#
# Lakin, Youssef, Cardelli, Phillips "Abstractions for DNA circuit design."
# J. R. Soc. Interface 9 (68) (2012) 470-486 (2012)
#
# Note: * Implements Figure 5: {A + B -> C + D}
#
#       * Not generalized.
#
# Coded by Stefan Badelt (badelt@caltech.edu)
#

# Define a global short toehold domain
global toehold = short();

# Write a class to define domains and structure of signal species
class formal(s) = "h t f" | ". . ."
  where { h = long() ; t = toehold ; f = long() };

macro pmac(s) = [ "i t j +" 
                | "( ( . +", 
                  "t* i*" 
                | ")  )",  
                  "t i"   
                | ". ." ]
  where {
    i = s.h;
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

class lakin_gate(r, p) 
  = [ "k + x + y l t* " 
    | "~ + ~ + ~ ~ .  ",
      "z t" 
    | "~ ."] + [m]
  where {
    t = toehold ; 
    [k, l, m] = flip(map(rmac, r), 3);
    l = reverse(l);
    [x, y, z] = flip(map(pmac, p), 3);
    y = reverse(y) 
  };

module rxn(r) = infty(gate) + infty(fuel)  + sum(map(infty,help))
  where {
    [gate, fuel, help] = lakin_gate(r.reactants, r.products) };

# Write the module *main* that applies *rxn* to the crn.
module main(crn) = sum(map(rxn, crn)) 
  where crn = irrev_reactions(crn);
