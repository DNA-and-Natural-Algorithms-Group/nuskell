#
# Qian, Soloveichik, Winfree: "Efficient Turing-universal computation with DNA
# polymers", DNA Computing and Molecular Programming 16, 2011.
#
# Note:   * modified Figure 1 (X + Y -> A + B) by introducing an irreversible 
#           step after reactants are consumed, but before products are released.
#         * implements Figure 2 (X + Y <=> A + B)
#
#         * generalized on the DNA level
#         * automatically combines corresponding irreversible reactions into 
#           one reversible reaction
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com)
#          Erik Winfree (winfree@caltech.edu)
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
       "b* a*" | ") )"]
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

macro revreactant(s)
     = ["a b c +" | ". ( ( +",
        "c* b*" | ") )"]
     where {
         a = s.hist;
         b = s.th;
         c = s.recog };

macro revproduct(s)
     = ["a b +" | "( ( +",
        "b* a*" | ") )"]
     where {
         a = s.th;
         b = s.hist };

class maingate_uni(r, p)
    = ["a f + g e + b f + f* c e* g* f* d e*"
     | "~ ( + ( ( + ~ ( + )  ~ )  )  )  ~ .",
       "e f g"
     | ". . ."]
    where {
        [a, temp1] = flip(map(reactant, r), 2);
        d = reverse(temp1);
        [b, temp2] = flip(map(product, p), 2);
        c = reverse(temp2);
        e = toehold;
        f = long();
        g = long() };

class maingate_rev(r, p)
    = [ "a b + c d e*"
      | "~ ~ + ~ ~ .",
        "w x + e* y z"
      | "~ ~ + .  ~ ~"]
    where {
        [a, temp1] = flip(map(reactant, r), 2);
        d = reverse(temp1);
        [b, temp2] = flip(map(product, p), 2);
        c = reverse(temp2);
        e = toehold;
        [w, temp3] = flip(map(revreactant, r), 2);
        z = reverse(temp3);
        [x, temp4] = flip(map(revproduct, p), 2);
        y = reverse(temp4) };

class reversefuelstrand(s)
    = "a b" | ". ."
    where {
        a = s.recog;
        b = s.th };

class forwardfuelstrand(s)
    = "a b" | ". ."
    where {
        a = s.th;
        b = s.hist };

module rxn(r) =
    if r.reversible then
        infty(g) + infty(h) + sum(map(infty, map(forwardfuelstrand, r.products))) + sum(map(infty, map(reversefuelstrand, r.reactants)))
        where 
            [g, h] = maingate_rev(r.reactants, r.products)
    else
        infty(g) + infty(h) + sum(map(infty, map(forwardfuelstrand, r.products))) + sum(map(infty, map(reversefuelstrand, r.reactants)))
        where
            [g, h] = maingate_uni(r.reactants, r.products);

module main(crn) = sum(map(rxn, rev_reactions(crn)))
