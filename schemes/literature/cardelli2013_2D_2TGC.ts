#
# Luca Cardelli "Two-Domain DNA Strand Displacement", Mathematical Structures
# in Computer Science. (2013)
#
#   Note: * The scheme implements Figures 4-6,8-10 (transducer, fork, 
#           join, 3-way join and the implicit n-way join) but it uses 
#           garbage collection as described in Figures 11, 12 and 14, 
#           i.e. avoids cooperative binding by using a second toehold.
#           - Figure 3-5 (Transducer): {X -> Y}
#           - Figure 6 (2-way Fork):  {X -> Y + Z}
#           - Figure 11,12 (2-toehold-2-way Join): {X + Y -> Z}
#           - Figure 14 (3-input join collector)
#           - Generalized n-way join collector as described in the text
#
#         * Does not implement Figure 7 (Catalyst):  {X + Y -> Y + Z}
#
#         * Generalized on the CRN level for { -> X}
#
#         * Generalized on the DNA level for higher order reactions.
#
# Coded by Stefan Badelt (badelt@caltech.edu)
#

global toehold = short();
global toehold2 = short();

class formal(s) = "t x"
                | ". ."
    where {
        t = toehold;
        x = long() };

class signal() = "t x"
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

macro gciter(r)
  = [ "x + t b + u" | "( + ( ( + (",
      "u* b* t* x*" | ")  )  )  ) ",
      "x t" | ". .",
      "b u" | ". .",
      "b + b* t*" | "( + ) ."]
  where {
    b = long();
    x = r.x ;
    u = toehold2;
    t = toehold
  };

macro gcgates(r)
  = [ "x y + u* y* b" 
    | "~ ( + .  )  ~" ] + multi + hx + hb + cb
  where {
    y = r[len(r)-1].x ;
    u = toehold2;
    multi = if len(r) > 1 then ["y + y* u*" | "( + ) ."] else [];
    [x, br, hx, hb, cb] = flip(map(gciter, reverse(tail(reverse(r)))), 5);
    b = reverse(br) 
  };

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

# shall we keep "t a" | ". ." as a fuel???
class gate2D_2T_GC(r,p)
    = [ "i + x u + a t + a + a* t* a* u* x* j t*"
      | "~ + ( ( + ( ( + ( + )  )  )  )  )  ~ . ",
        "u a" | ". .", "t a" | ". ."] + noutput(r, reverse(p), a) + ifw + gcg
    where {
      [i, jr, ifw] = flip(map(imac,reverse(tail(reverse(r)))),3);
      gcg = if len(r) > 1 then gcgates(tail(r)) else [];
      j = reverse(jr);
      x = r[len(r)-1].x;
      u = toehold2;
      a = long();
      t = toehold };

module reaction(r) =
    if len(r.reactants) == 0 then
      sum(map(infty,gate2D_noGC([i], r.products)+[i]))
        where { i = signal() }
    else if len(r.reactants) == 1 then
      sum(map(infty, gate2D_noGC(r.reactants, r.products)))
    else 
      sum(map(infty, gate2D_2T_GC(r.reactants, r.products)));

module main(crn) = sum(map(reaction, crn))
  where crn = irrev_reactions(crn)
