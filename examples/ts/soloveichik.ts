#
# David Soloveichik's translation scheme from "DNA as a universal substrate
# for chemical kinetics", Proceedings of the National Academy of Sciences, 
# 107: 5393-5398, 2010.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
#

class formal(s) = "? a b c"
                | ". . . ."
    where {
        a = short() ; 
        b = long() ;
        c = short() };

class ugate(s, l)
     = ["b c d + c* b* a*"        # G_i
      | "( ( ~ + )  )  .",
        "e + f c*"                # T_i
      | "~ + ~ ."]
    where {
        a = s.a ;
        b = s.b ;
        c = s.c ;
        [d, e, g] = flip(map(gmac, l), 3) ;
        f = reverse(g) };

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

module unimolecular(r) = infty(g) + infty(t)
    where
        [g, t] = ugate(r.reactants[0], r.products);

class bgate(s1, s2, l)
    = ["b c d + e f g + f* e* d* c* b* a*"       # L_i
     | "( ( ( + ( ( ~ + )  )  )  )  )  .",
       "h + i f*"                                # T_i
     | "~ + ~ .",
       "b c d"                                   # B_i
     | ". . ."]
    where {
        a = s1.a ;
        b = s1.b ;
        c = s1.c ; 
        d = s2.a ;
        e = s2.b ;
        f = s2.c ;
        [g, h, j] = flip(map(gmac, l), 3) ;
        i = reverse(j) };

module bimolecular(r) = infty(l) + infty(t) + infty(b)
    where
        [l, t, b] = bgate(r.reactants[0], r.reactants[1], r.products);

module rxn(r) =
    if len(r.reactants) == 1 then
        unimolecular(r)
    elseif len(r.reactants) == 2 then
        bimolecular(r)
    else r[0];

module main(crn) = sum(map(rxn, crn))
    where
        crn = irrev_reactions(crn)
