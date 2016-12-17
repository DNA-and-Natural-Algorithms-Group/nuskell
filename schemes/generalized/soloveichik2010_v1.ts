#
# David Soloveichik's translation scheme from "DNA as a universal substrate
# for chemical kinetics", Proceedings of the National Academy of Sciences, 
# 107: 5393-5398, 2010.
#
#   Note: * Implements the scheme as shown in Figures 2 and 3 and generalized 
#           to arbitrary reactions using the following rules:
#           - ' -> A' : 'f -> A' where f is a fuel species
#           - 'A -> ' : 'A -> f' where f is a fuel species
#           - reactions with arity > 2 are implemented by elongating fuel 
#               complexes and fuel (helper) strands
#
#   Note: * The choice for generalization of 'A-> ' can easily be optimized by 
#           omitting f, see comments in the "rxn" function
#
#   Note: * Higher arity reactions may also be implemented as: 
#             A + B + C -> X  becomes  A + B <=> i; i + C -> X
#           This would increase the number of complexes and reaction steps, but 
#             reduce the maximal length of strands.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com)
# modified by Stefan Badelt (badelt@caltech.edu)
#

class formal(s) = "? a b c"
                | ". . . ."
    where {
        a = short() ; 
        b = long() ;
        c = short() };

class signal() = "? a b c"
                | ". . . ."
    where {
        a = short() ; 
        b = long() ;
        c = short() };

macro gmac(s)
     = ["d a"
      | ". .",
        "d a b c +"
      | "( ( . . +",
        "a* d*"
      | ")  )"]
    where {
        d = long() ;
        a = s.a ; 
        b = s.b ;
        c = s.c };

macro inputmac(s, i)
     = ["b c a +"
      | "( ( ( +",
        "a* c* b*"
      | ")  )  )",
        "b c a"
      | ". . ."]
    where {
        b = s[i].b ;
        c = s[i].c ;
        a = s[i+1].a };

class maingate(s, l)
    = ["b x y g + y* x* c a*"                    # L_i
     | "~ ( ( ~ + )  )  ~ .",
       "h + i y*"                                # T_i
     | "~ + ~ .", bi]
    where {
        x = s[len(s)-1].b;
        y = s[len(s)-1].c;
        a = s[0].a ;
        [b, d, bi] = flip(map2(inputmac, s, range(len(s)-1)), 3) ;
        c = reverse(d) ;
        [g, h, j] = flip(map(gmac, l), 3) ;
        i = reverse(j) };

module rxn(r) =
  if len(r.reactants) == 0 then
    infty(i) + infty(l) + infty(t) + sum(map(infty,b))
      where {
        i = signal();
        [l, t, b] = maingate([i], r.products)}

  # Comment out the following elseif statement for optimization
  elseif len(r.products) == 0 then
    infty(l) + infty(t) + sum(map(infty,b))
      where {
        i = signal();
        [l, t, b] = maingate(r.reactants, [i])}
 
  else
    infty(l) + infty(t) + sum(map(infty,b))
      where 
        [l, t, b] = maingate(r.reactants, r.products);

module main(crn) = sum(map(rxn, crn))
    where
        crn = irrev_reactions(crn)

