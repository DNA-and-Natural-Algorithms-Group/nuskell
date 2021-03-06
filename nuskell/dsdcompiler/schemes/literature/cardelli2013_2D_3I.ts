#
# Luca Cardelli "Two-Domain DNA Strand Displacement", Mathematical Structures
# in Computer Science. (2013)
#
#   Note: * The scheme implements Figures 4-10 (transducer, fork, 
#           catalyst, join, 3-way join and the implicit n-way join)
#           but it uses a 3-domain strand to introduce an irreversilbe
#           step after all reactants have been consumed (Figure 15).
#           - Figure 3-5 (Transducer): {X -> Y}
#           - Figure 6 (2-way Fork):  {X -> Y + Z}
#           - Figure 7 (Catalyst):  {X + Y -> Y + Z}
#           - Figure 8,9 (2-way Join): {X + Y -> Z}
#           - Figure 10 (3-way Join): {W + X + Y -> Z}
#           - Figure 15 (3-domain transducer modification)
#             + same 3-domain mechanism for every other reation.
#
#         * Generalized on the CRN level for { -> X}
#
#         * Uses garbage collection with cooperative binding, which cannot be 
#           enumerated by peppercorn -> pathway equivalence fails.
#
#         * Generalized on DNA level for higher order reactions, but supports
#           only one catalyst optimization. e.g. X + Y + Z -> Z + Y + W will 
#           only optimize for Z, but not for Y.
#
#         * Note, the Catalyst case does not apply to {C -> C + A}, as it is
#           not obvious from the paper.
#
#         * Needs to be enumerated using --reject-remote
#
# Coded by Stefan Badelt (badelt@caltech.edu)
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

macro gcgates(r)
  = [ "t b + " | " ( ( + ",
      "b* t* " | " ) ) ",
      "b t" | ". .",
      " b y + t* y* b* t*"
    | " ( ( + .  )  )  . " ]
  where {
    b = long();
    y = r.x;
    t = toehold };

class noutput(r, p, a)
    = [ "x + gc1 out1 t a + t* a* t* out2 gc2 x*"
      | "( +  ~   ~   ( ( + .  )  )   ~    ~  )"] + l + gl + gcg
    where {
        x = (r[0]).x;
        [out1, out2r, l] = flip(map(output, p), 3);
        [gc1, gc2r, gl, gcg] = flip(map(gcgates, tail(r)),4);
        out2 = reverse(out2r);
        gc2 = reverse(gc2r);
        t = toehold };

class noutput_cat(r, p, a)
    = [ "x + gc1 out1 t a + t* a* t* out2 gc2 x*"
      | "( +  ~   ~   ( ( + .  )  )   ~    ~  )"] + l + gl + gcg
    where {
        x = (r[0]).x;
        [out1, out2r, l] = flip(map(output, p), 3);
        l = reverse(tail(reverse(l)));
        r = reverse(tail(reverse(r)));
        [gc1, gc2r, gl, gcg] = flip(map(gcgates, tail(r)),4);
        out2 = reverse(out2r);
        gc2 = reverse(gc2r);
        t = toehold };

macro imac(r) 
  = [ "x t + " | "( ( + ", 
      "t* x*" | ") )",
      "x t" | ". ."]
    where { t = r.t; x = r.x };

class gate2D_GC(r,p)
    = [ "i + b + a t + a + a* t* a* b* j t*"
      | "~ + ( + ( ( + ( + )  )  )  )  ~ . ",
        "t b a"
      | ". . ."] + nout
    where {
      a = long();
      nout = if len(r) > 1 and len(p) > 0 and r[-1] == p[0] then noutput_cat(r, reverse(p), a) else noutput(r, reverse(p), a) ;
      [i, jr, ifw] = flip(map(imac,r),3);
      j = reverse(jr);
      b = long();
      t = toehold };

module reaction(r) =
    if len(r.reactants) == 0 then
      sum(map(infty,gate2D_GC([i], r.products)+[i]))
        where { i = fuel() }
    else
      sum(map(infty, gate2D_GC(r.reactants, r.products)));

module main(crn) = sum(map(reaction, crn))
  where crn = irrev_reactions(crn)
