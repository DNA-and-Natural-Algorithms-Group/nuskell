#
# Qian, Soloveichik, Winfree: "Efficient Turing-universal computation with DNA
# polymers", DNA Computing and Molecular Programming 16, 2011.
#
# Note:   * implements Figure 1 (X + Y -> A + B) 
#         * implements Figure 2 (X + Y <=> A + B)
#
#         * generalized on the DNA level
#         * automatically combines corresponding irreversible reactions into 
#           one reversible reaction
#
#         * pathway and bisimulation incorrect for irreversible reactions.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com)
#          Stefan Badelt (badelt@caltech.edu)
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

macro reactant_rev(s)
    = ["a b +" | "( ( +",
       "b* a*" | ") )",
       "a b" | ". ."]
    where {
        b = s.hist;
        a = s.th };

macro product(s)
    = ["a b c +" | "( ( . +",
       "b* a*" | ") )"]
    where {
        a = s.hist;
        b = s.th;
        c = s.recog };

macro product_rev(s)
    = ["a b c +" | ". ( ( +",
       "c* b*" | ") )"]
    where {
        a = s.hist;
        b = s.th;
        c = s.recog };

class maingate(r, p)
    = [["a b f + f* c d e*"
     | "~ ~ ( + )  ~ ~ .",
       "e f"
     | ". ."], extra]
    where {
        [a, temp, extra] = flip(map(reactant, r), 3);
        d = reverse(temp);
        [b, temp2] = flip(map(product, p), 2);
        c = reverse(temp2);
        e = toehold;
        f = long() };

class maingate_rev(r, p)
    = ["a b + c d e*"
     | "~ ~ + ~ ~ .",

       "k i + e* j l"
     | "~ ~ + .  ~ ~",

       extra]
    where {
        [a, temp, extra] = flip(map(reactant, r), 3);
        d = reverse(temp);
        [b, temp2] = flip(map(product, p), 2);
        c = reverse(temp2);

        [k, temp3] = flip(map(product_rev, r), 2);
        l = reverse(temp3);
        [i, temp4, waste] = flip(map(reactant_rev, p), 3);
        j = reverse(temp4);

        e = toehold};

class suppgate(s)
    = "a b" | ". ."
    where {
        a = s.th;
        b = s.hist };

module rxn(r) =
    if r.reversible then
        infty(fw) + infty(bw) + sum(map(infty, map(suppgate, r.products))) + sum(map(infty, gates))
        where
            [fw, bw, gates] = maingate_rev(r.reactants, r.products)
    else
        infty(g) + infty(h) + sum(map(infty, map(suppgate, r.products))) + sum(map(infty, gates))
        where
            [[g, h], gates] = maingate(r.reactants, r.products);

module main(crn) = sum(map(rxn, rev_reactions(crn)))
