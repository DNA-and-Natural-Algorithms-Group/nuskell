#
# Luca Cardelli's translation scheme from "Two-Domain DNA Strand
# Displacement", Mathematical Structures in Computer Science.
# Can be found at http://lucacardelli.name/
#
#   Note: * The scheme implements Figures 4-10 (transducer, fork, 
#           catalyst, join, 3-way join and the implicit n-way join)
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
# generalized by Stefan Badelt (badelt@caltech.edu)
#

global toehold = short();

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

macro imac(r) 
  = [ "x t + " | "( ( + ", 
      "t* x*" | ") )" ]
    where { t = r.t; x = r.x };

class gate2D_GC(r,p)
    = [ "i + a t + a + a* t* a* j t*"
      | "~ + ( ( + ( + )  )  )  ~ . ",
        "t a"
      | ". ."] + noutput(r, reverse(p), a)
    where {
      [i, jr] = flip(map(imac,r),2);
      j = reverse(jr);
      a = long();
      t = toehold };

module reaction(r) =
    if len(r.reactants) == 0 then
      sum(map(infty,gate2D_GC([i], r.products)+[i]))
        where { i = signal() }
    else
      sum(map(infty, gate2D_GC(r.reactants, r.products)));

module main(crn) = sum(map(reaction, crn))
