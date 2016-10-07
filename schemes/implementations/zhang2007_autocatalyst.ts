# Zhang et al. (2007) "Engineering Entropy-Driven Reactions and Networks Catalyzed by DNA"

# A DSD implementation of an autocatalyst reaction (Figure 4)

# Coded by Stefan Badelt (badelt@caltech.edu)

# -----------------------------------------
# Input implementation CRN (Figure 4)
# -----------------------------------------
# A + S -> I + SB
# I + F -> A + A + W
# SB + SR -> FQ + TET
# formal = {A, SB, TET}
# constant = {S, F, SR}
# -----------------------------------------

# These are the same as for the "Zhang 2007 catalyst" system
global d2b = short(len:6);
global d2c = long(len:12);
global d3  = short(len:4);

# This might look confusing, but it is including the information 
# from the Erratum (d4t is a subsequence of the *original* d4a:
global d4t  = long(len:7);
global d4a  = long(len:3);
global d4b  = long(len:6);

global d6  = long(len:16);

class formal(s) = 
  if s.name == 'A' then A(s)
  else if s.name == 'SB' then SB(s)
  else if s.name == 'TET' then TET(s)
  else empty() ;
    #where void = print('ignoring formal species:', s.name);

class constant(s) = 
  if s.name == 'S' then S(s)
  else if s.name == 'F' then F(s)
  else if s.name == 'SR' then SR(s)
  else empty() ;
    #where void = print('ignoring constant species:', s.name);

class A(s) = "d4t d4a d4b d2b d2c" | ". . . . .";

class SB(s) = "d6 d3 d4t d4a d4b" | ". . . . .";

class TET(s) = "d6 d3" | ". .";

class S(s) = 
    " d4t d4a d4b d2b d2c + d6 d3 d4t d4a d4b + d2b* d4b* d4a* d4t* d3* d2c* d2b* d4b* "
  | "  .   .   (   (   (  + .  (   (   (   (  +  .    )    )    )    )   )    )    )   ";

class F(s) = "d4b d2b d2c d3 d4t d4a d4b" | ". . . . . . .";

class SR(s) = 
    " d6 d3 + d4t* d3* d6* "
  | " (  (  +  .    )   )  ";

 
module rxn(r) = sum(map(infty, r.reactants + r.products));

module main(crn) = sum(map(rxn, crn));

