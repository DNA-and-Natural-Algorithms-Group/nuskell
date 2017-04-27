#
# Lakin, Stefanovic, Phillips "Modular verification of chemical reaction
# network encodings via serializability analysis" (2016), TCS
#
#   Note: * The scheme implements Figure 8 from the above publication. 
#           {x + y -> y + b}
#
#         * generalized on the CRN level for all but bimolecular reactions.
#
#         * requires autocatalytic format, otherwise incorrect, BUT:
#           as e.g. {A->B} is implemented as {A+f->f+B} it is correct for 
#           a wider class of CRNs!
#
#         * enumerate using --reject-remote !!!
#
# Coded by Stefan Badelt (badelt@caltech.edu)
#

global toehold = short();

class formal(s) = "t x"
                | ". ."
    where {
        t = toehold;
        x = long() };

class inter() = "t x"
               | ". ."
    where {
        t = toehold;
        x = long() };

global gfuel = inter();

class twooutput(r, p, c, a)
    = [ "x + t b + c + t y + t a + t* a* t* y* t* c* b* t* x* "
      | "( + ( ( + ( + ( ( + ( ( + .  )  )  )  )  )  )  )  )  ",
      # "y t" | ". .", # same as twoinput "y t" for A + Y <=> Y + B
        "b c t" | ". . ."]
    where {
      x = r[0].x;
      y = p[0].x;
      b = p[1].x;
      t = toehold };

class twoinput(r, c, a) =
    [ "x t + y t + c + a t + a + a* t* a* c* t* y* t* x* t*"
    | "( ( + ( ( + ( + ( ( + ( + )  )  )  )  )  )  )  )  . ",
      #"t a" | ". .",  # not in original Figure 15
      #"x t" | ". .",  # not in original Figure 15
      #"y t" | ". .",  # not in original Figure 15
      "t c a" | ". . ."]
    where {
      x = r[0].x;
      y = r[1].x;
      t = toehold };

class lakin3D_cat_rev(r,p) = lakin3D_cat(r,p) + lakin3D_cat(p,r) ;

class lakin3D_cat(r,p) = twoinput(r, c, a) + twooutput(r, p, c, a)
  where { 
    c = long(); a = long() };

module prep(r, p) = narity
  where 
    narity = if len(r) == 2 and len(p) == 2 then
                lakin3D_cat(r,p)
             else if len(r) > 2 then
                lakin3D_cat_rev([r[0], r[1]], [gfuel, i]) + prep([i]+tail(tail(r)), p) + [gfuel]
                  where i = inter()
             else if len(p) > 2 then
                lakin3D_cat(r, [p[0], j]) + prep([j, gfuel], tail(p)) + [gfuel]
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
    
module main(crn) = sum(map(reaction, crn))
  where crn = irrev_reactions(crn)
