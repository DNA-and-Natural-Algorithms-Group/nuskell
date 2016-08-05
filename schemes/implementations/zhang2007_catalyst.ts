# Zhang et al. (2007) "Engineering Entropy-Driven Reactions and Networks Catalyzed by DNA"

# A case-by-case example for the catalyst reaction (Figure 1A + 1D)
# Note: Domain 2 is actually contains two toeholds (2a, 2b), and 2c

# Coded by Stefan Badelt (badelt@caltech.edu)

# -----------------------------------------
# Input implementation CRN (Figure 1A + 1D)
# -----------------------------------------
# C + S + F -> C + SB + OB + W
# OB + OR -> RQ + ROX
# formal = {C, ROX}
# constant = {S, F, OR}
# -----------------------------------------

global d1  = long(len:10);
global d2a = short(len:6);
global d2b = short(len:6);
global d2c = long(len:12);
global d3  = short(len:4);
global d4  = long(len:16);
global d5  = short(len:6);
global d6  = long(len:16);

class formal(s) = 
  if s.name == 'C' then C(s)
  else if s.name == 'ROX' then ROX(s)
  else empty() ;
    #where void = print('ignoring formal species:', s.name);

class constant(s) = 
  if s.name == 'S' then S(s)
  else if s.name == 'F' then F(s)
  else if s.name == 'OR' then OR(s)
  else empty() ;
    #where void = print('ignoring constant species:', s.name);
 
class C(s) = "d4 d5" | ". .";

class F(s) = "d2a d2b d2c d3 d4" | ". . . . .";

class S(s) = 
    " d1 d2a d2b d2c + d6 d3 d4 + d5* d4* d3* d2c* d2b* d2a* " 
  | " .   (   (   (  + .  (  (  +  .   )   )   )    )    )   ";

class OR(s) = 
    " d1 d2a + d2b* d2a* d1* "
  | " (   (  +  .    )    )  ";

class ROX(s) = " d1 d2a " | " . . ";

module rxn(r) = sum(map(infty, r.reactants + r.products));

module main(crn) = sum(map(rxn, crn));


