#
# David Soloveichik's translation scheme from "DNA as a universal substrate
# for chemical kinetics", Proceedings of the National Academy of Sciences, 
# 107: 5393-5398, 2010.
#
# Generalized to implement reactions of arbitrary arity.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
#



class formal(s) = "? a b c"
                | "? . . ."
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

class spongate(l)
    = ["y g"
     | ". ~",
       "h + i y*"                                # T_i
     | "~ + ~ ."]
    where {
        y = short() ;
        [g, h, j] = flip(map(gmac, l), 3) ;
        i = reverse(j) };

module normalrxn(r) = infty(l) + infty(t) + sum(map(infty, b))
    where
        [l, t, b] = maingate(r.reactants, r.products);

module sponrxn(r) = infty(o) + infty(t)
    where
        [o, t] = spongate(r.products);

module rxn(r) = if len(r.reactants) == 0 then sponrxn(r) else normalrxn(r);

module main(crn) = sum(map(rxn, crn))
    where
        crn = irrev_reactions(crn)
