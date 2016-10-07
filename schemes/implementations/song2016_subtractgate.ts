# Song et al. (2016) "Analog Computation by DNA Strand Displacement Circuits"
# A case-by-case example for the subtraction gate (Figure 7)

# Coded by Stefan Badelt (badelt@caltech.edu)

# ----------------------------------
# Format version 1: 'implementation'
# ----------------------------------
# Gs + Is2 <=> Gsp + i
# i + Ds -> w1 + w2
# Gsp + Is1 -> w3 + w4
# formal = {Is1, Is2}
# constant = {Gs, Ds}
# ----------------------------------
# => verifier: compare enumerated CRN to implementation CRN
# ... write a class for every formal and constant species that 
#     uses globaly shared domains

global a  = short();
global v1 = long();
global v2 = long();

# Formal takes a Species Object as input
class formal(s) = 
  if s.name == "Is1" then Is1(s)
  else if s.name == "Is2" then Is2(s)
  else print(s.name);

class Is1(s) = "h a v1" | ". . ."
  where h = long();

class Is2(s) = "h a v2" | ". . ."
  where h = long();

# Constant assumes a Species Object as input
class constant(s) = strands
  where strands = 
    if s.name == "Gs" 
      then Gs(s)
    else if s.name == "Ds" 
      then Ds(s)
    else empty();

class Gs(s) = "v2 a + v1 + v1* a* v2* a*"
            | "(  ( + (  +  )  )   )  . ";

class Ds(s) = "v2 + a* v2*" | "( + . )";

module rxn(r) = sum(map(infty, r.reactants + r.products));
  #where void = print(r.reactants, r.products);         #-^- remove

module main(crn) = sum(map(rxn, crn));

#-----
# 
# function solution(r) = map(infty, const)
#   where const = map(isconst, r.reactants + r.products)
# 
# function isconst(s) = if s.
#   
# 
# module rxn(r) = infty(
#   where {
#     gates = map(constant, r.reactants + r.products);
#     void = print('aaaa', gates)};
# 
# module main(crn) = sum(map(rxn, crn)) 


# ----------------------------
# Format version 2: 'condense'
# ----------------------------
# B + X <=> Y + i
# Z + i -> 
# A + Y -> 
# formal { A, B }
# constant { X, Y, Z }
# ----------------------------
# => problem: how does Y know that it will gain domains from A & B?
# => there are two different pathways for int+const and formal+const

# ----------------------------
# Format version 3: 'minimal'
# ----------------------------
# B + X -> Y 
# A + Y -> 
# formal { A, B }
# constant { X, Y }
# ----------------------------
# ===> problem: how does Y know that it will gain domains from A & B?

# ----------------------------
# Format version 3: 'tiny'
# ----------------------------
# A -> i 
# B + i -> 
# formal { A, B }
# ----------------------------
# ===> problem: how does Y know that it will gain domains from A & B?

# --------------------------
# Format version 3: 'analog'
# --------------------------
# A = A - B
# --------------------------
# ===> problem: requires a different framework...

