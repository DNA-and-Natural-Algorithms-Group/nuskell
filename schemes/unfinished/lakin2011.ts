# Lakin et. al (2011) Abstractions for DNA circuit design
# unbuffered gate (e.g. Figure 5)

# written by Stefan Badelt (badelt@caltech.edu)

global toehold = short();

class formal(s) = "h t X" 
                | ". . ."
  where {
    h = long() ;
    t = toehold ;
    X = long() 
  };

macro pmac(s) = [ "i t j +" 
                | "( ( . +", 
                  "t* i*" 
                | ")  )",  
                  "t i"   
                | ". ." ]
  where {
    i = s.h;
    t = s.t;
    j = s.X 
  };

macro rmac(s) = 
  [ "a t i + " 
  | "( ( . + ", 
    "t* a*" 
  | ")  )" ]
  where {
    a = s.X ;
    t = s.t ;
    i = long() 
  };

class lakin_gate(r, p) 
  = [ "k + x + y l t* " 
    | "~ + ~ + ~ ~ .  ",
      "z t" 
    | "~ ."]
  where {
    t = toehold ; 
    [k, l] = flip(map(rmac, r), 2);
    l = reverse(l);
    [x, y, z] = flip(map(pmac, p), 3);
    y = reverse(y) 
  };

module rxn(r) = infty(gate) + infty(fuel) 
  where 
    [gate, fuel] = lakin_gate(r.reactants, r.products) ;

module main(crn) = sum(map(rxn, crn)) 
  where 
    crn = irrev_reactions(crn)
