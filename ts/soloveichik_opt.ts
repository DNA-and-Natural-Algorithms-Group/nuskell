#
# David Soloveichik's translation scheme from "DNA as a universal substrate
# for chemical kinetics", Proceedings of the National Academy of Sciences, 
# 107: 5393-5398, 2010.
#
# Generalized to implement reactions of arbitrary arity.
# Uses delayed choice to optimize number of species
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
     | "~ + ~ .", bi, h, j]
    where {
        x = s[len(s)-1].b;
        y = s[len(s)-1].c;
        a = s[0].a ;
        [b, d, bi] = flip(map2(inputmac, s, range(len(s)-1)), 3) ;
        c = reverse(d) ;
        [g, h, j] = flip(map(gmac, l), 3) ;
        i = reverse(j) };

class maingateopt(s, l, h, j)
    = ["h + i y*"                                # T_i
     | "~ + ~ ."]
    where {
        y = s[len(s)-1].c;
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

module normalrxn(r) = [infty(l) + infty(t) + sum(map(infty, b)), h, j]
    where
        [l, t, b, h, j] = maingate(r.reactants, r.products);

module normalrxnopt(r, h, j) = [infty(t), h, j]
    where
        [t] = maingateopt(r.reactants, r.products, h, j);

module sponrxn(r) = [infty(o) + infty(t), 0, 0]
    where
        [o, t] = spongate(r.products);

function prefix(l1, l2) = if l1 == [] then 1 elseif l2 == [] or l1[0] != l2[0] then 0 else prefix(tail(l1), tail(l2));
function cutprefix(l, n) = if n == 0 then [] else [l[0]] + cutprefix(tail(l), n - 1);

# if there are two reactions R_1 and R_2 such that they have the same reactants
# and the products of R_1 is a prefix of the products of R_2
# then the gate L_i can be shared 
function optimize(x, i) =
    if len(crn[i].reactants) > 0 and
       prefix(r.products, crn[i].products) and
       crn[i].reactants == r.reactants then [cutprefix(previous[i][0], len(r.products)), cutprefix(previous[i][1], len(r.products))] else 0
    where {
        crn = x[0];
        previous = x[1];
        r = x[2] };

function trim(l) = if l == [] then [] elseif l[0] == 0 then trim(tail(l)) else [l[0]] + trim(tail(l));

module helper(crn, i, previous) =
    if i >= len(crn) then
        empty
    else
        result + helper(crn, i+1, previous + [[h, j]])
    where {
        r = crn[i];
        [result, h, j] = 
            if len(r.reactants) == 0 then
                sponrxn(r)
            else
                (if len(x) > 0 and len(r.products) > 0 then
                    normalrxnopt(r, x[0][0], x[0][1])
                else
                    normalrxn(r))
                where 
                    x = trim(map2(optimize, [crn, previous, crn[i]], range(i)))};

module main(crn) = helper(crn, 0, [])
    where
        crn = irrev_reactions(crn)
