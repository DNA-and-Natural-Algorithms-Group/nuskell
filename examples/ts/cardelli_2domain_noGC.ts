#
# Luca Cardelli's translation scheme from "Two-Domain DNA Strand
# Displacement", Mathematical Structures in Computer Science.
# Can be found at http://lucacardelli.name/
#
#   NOTE : The garbage collection module was omitted in this implementation.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
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

class outgate(r, p, a)
    = [ "x + out1 t a + t* a* t* out2 x*"
      | "( +  ~   ( ( + .  )  )   ~   )"] + l
    where {
        x = (r[0]).x;
        [out1, out2r, l] = flip(map(output, p), 3);
        out2 = reverse(out2r);
        t = toehold };

class oneinput(r, p)
    = [ "x t + a t + a + a* t* a* t* x* t*"
      | "( ( + ( ( + ( + )  )  )  )  )  . ",
        "x t"
      | ". .",
        "t a"
      | ". ."] + outgate(r, p, a)
    where {
        a = long();
        x = (r[0]).x;
        t = toehold };

class twoinput(r, p)
    = [ "x t + y t + a t + a + a* t* a* t* y* t* x* t*"
      | "( ( + ( ( + ( ( + ( + )  )  )  )  )  )  )  . ",
        "t a"
      | ". .",
        "x t"
      | ". .",
        "y t"
      | ". ."] + outgate(r, p, a)
    where {
        x = (r[0]).x;
        y = (r[1]).x;
        t = toehold;
        a = long() };

module reaction(r) =
    if len(r.reactants) == 1 then
        sum(map(infty,oneinput(r.reactants, r.products)))
    elseif len(r.reactants) == 2 then
        sum(map(infty,twoinput(r.reactants, r.products)))
    else 
      abort('reaction-type not implemented for this scheme');

module main(crn) = sum(map(reaction, crn))
