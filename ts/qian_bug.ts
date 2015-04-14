#
# Lulu Qian's translation scheme from "Efficient Turing-universal
# computation with DNA polymers", Lecture Notes in Computer Science,
# 6518: 123-140, 2011.
# 
#   Note : this was a scheme that the authors of the above paper discovered
#          while they were working on that paper, but it was later shown to
#          have a bug and was not included in the published version of the
#          paper. It was reproduced here for demonstrative purposes.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
#

global toehold = short();

class formal(s) = "hist th recog"
                | ". . ."
    where {
        hist = long();
        th = toehold;
        recog = long()};

macro reactant(s)
    = ["a b +" | "( ( +",
       "b* a*" | ") )",
       "a b" | ". ."]
    where {
        a = s.recog;
        b = s.th };

macro product(s)
    = ["a b c +" | "( ( . +",
       "b* a*" | ") )"]
    where {
        a = s.hist;
        b = s.th;
        c = s.recog };

class maingate(r, p)
    = [["a b + c d e*"
     | "~ ~ +  ~ ~ ."], extra]
    where {
        [a, temp, extra] = flip(map(reactant, r), 3);
        d = reverse(temp);
        [b, temp2] = flip(map(product, p), 2);
        c = reverse(temp2);
        e = toehold;
        f = long() };

class suppgate(s)
    = "a b" | ". ."
    where {
        a = s.th;
        b = s.hist };

module rxn(r) = infty(g) + sum(map(infty, map(suppgate, reverse(tail(reverse(p)))))) + infty("t b t" | ". . .") + sum(map(infty, gates))
    where {
        p = if len(r.products) == 0 then [formal(0)] else r.products;
        [[g], gates] = maingate(r.reactants, p);
        t = toehold;
        b = (p[len(p) - 1]).hist };

module main(crn) = sum(map(rxn, irrev_reactions(crn)))
