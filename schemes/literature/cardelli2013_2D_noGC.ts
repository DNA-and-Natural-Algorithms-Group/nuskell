#
# Luca Cardelli "Two-Domain DNA Strand Displacement", Mathematical Structures
# in Computer Science. (2013)
#
#   Note: * The scheme implements Figures 4-10 (transducer, fork, 
#           catalyst, join, 3-way join and the implicit n-way join)
#           - Figure 3-5 (Transducer): {X -> Y}
#           - Figure 6 (2-way Fork):  {X -> Y + Z}
#           - Figure 7 (Catalyst):  {X + Y -> Y + Z}
#           - Figure 8,9 (2-way Join): {X + Y -> Z}
#           - Figure 10 (3-way Join): {W + X + Y -> Z}
#
#         * All reactios are implemented *without* garbage collection. This 
#           includes the cooperative binding complexes and additional domains
#           on the produce complexes.
#
#         * Intermediate ('garbage') species on react complexes are fuels 
#           in this implementation.
#
#         * Generalized on the CRN level for { -> X}
#
#         * Generalized on DNA level for higher order reactions.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
#          Stefan Badelt (badelt@caltech.edu)
#

global toehold = short();

class formal(s) = "t x"
                | ". ."
    where {
        t = toehold;
        x = long() };

class fuel() = "t x"
               | ". ."
    where {
        t = toehold;
        x = long() };

macro output(s)
    = ["t x +"
     | "( ( +",
       "x* t*"
     | ")  )",
       "x t"
     | ". ."]
    where {
        t = toehold;
        x = s.x };

class noutput(r, p, a)
    = [ "x + out1 t a + t* a* t* out2 x*"
      | "( +  ~   ( ( + .  )  )   ~   )"] + l
    where {
        x = (r[0]).x;
        [out1, out2r, l] = flip(map(output, p), 3);
        out2 = reverse(out2r);
        t = toehold };

macro imac(r) 
  = [ "x t + " | "( ( + ", 
      "t* x*" | ") )",
      "x t" | ". ."]
    where { t = r.t; x = r.x };

class gate2D_noGC(r,p)
    = [ "i + a t + a + a* t* a* j t*"
      | "~ + ( ( + ( + )  )  )  ~ . ",
        "t a"
      | ". ."] + noutput(r, reverse(p), a) + ifw
    where {
      [i, jr, ifw] = flip(map(imac,r),3);
      j = reverse(jr);
      a = long();
      t = toehold };

module reaction(r) =
    if len(r.reactants) == 0 then
      sum(map(infty,gate2D_noGC([i], r.products)+[i]))
        where { i = fuel() }
    else
      sum(map(infty, gate2D_noGC(r.reactants, r.products)));

module main(crn) = sum(map(reaction, crn))
  where crn = irrev_reactions(crn)
