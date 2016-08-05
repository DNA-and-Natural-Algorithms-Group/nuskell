#
# Lulu Qian's translation scheme from "Efficient Turing-universal
# computation with DNA polymers", Lecture Notes in Computer Science,
# 6518: 123-140, 2011.
#
#   NOTE : This is the corrected version of the scheme in the above paper
#          which was provided to us by the authors of that paper.
#          This scheme is referred to as "Qian_Fixed" in Seung Woo Shin's
#          Master's Thesis.  It is conceptually very similar to 
#          the slightly simplified version that appears in the Journal
#          of the Royal Society Interface revision of Qian et al.
#   NOTE:  Files qian_fixed.ts and qian_mastersthesis.ts are identical
#          and always will be.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
#

global toehold = short();

class formal(s) = "hist th recog"
                | ". . ."
    where {
        hist = long();
        th = toehold;
        recog = long() };

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
    = [["a f + g e + b + f + f* c e* g* f* d e*"
     | "~ ( + ( ( + ~ + ( + ) ~ ) ) ) ~ .",
       "e f g"
     | ". . ."], extra]
    where {
        [a, temp, extra] = flip(map(reactant, r), 3);
        d = reverse(temp);
        [b, temp2] = flip(map(product, p), 2);
        c = reverse(temp2);
        e = toehold;
        f = long();
        g = long() };

class suppgate(s)
    = "a b" | ". ."
    where {
        a = s.th;
        b = s.hist };

module rxn(r) = infty(g) + infty(h) + sum(map(infty, map(suppgate, r.products))) + sum(map(infty, gates))
    where
        [[g, h], gates] = maingate(r.reactants, r.products);

module main(crn) = sum(map(rxn, crn))
