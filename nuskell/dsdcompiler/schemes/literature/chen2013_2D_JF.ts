#
# Chen, Dalchau, Srinivas, Phillips, Cardelli, Soloveichik, Seelig
# "Programmable chemical controllers made from DNA", Nature Nanotech., 2013.  
#
#   Note: * Implements supplementary Figure S7.6 fork and join modules.
#           {A->R; R->B; A+B->R; R->B+C; A+B+C->R; R->B+C+D}
#         * generalized on the DNA level higher order reactions and for {X->}
#         * generalized on the CRN level for {->X}
#         * includes the delayed choice optimization suggested by the authors
#
# Coded by Stefan Badelt (badelt@caltech.edu)
#

class formal(s) = "t x" | ". ."
    where { t = short(); x = long() };

class fuel() = "t x" | ". ."
    where { t = short(); x = long() };

class flux_domains() = "tr dr tq di" | ". . . ."
  where {
    tr = short(); 
    dr = long(); 
    tq = short();
    di = long() };

macro chen2D_I(s) = 
  [ "tb + b" | "( + (", 
    "b* tb*" | ") )" ]
    where { tb = s.t; b = s.x };

macro chen2D_O(s) = 
  [ "tb b + " | "( ( +", 
    "b* tb*" | ") )"]
    where { tb = s.t; b = s.x };

macro chen2D_fo(data, c) = "d t"|". ."
  where {
    i = data[0];
    tr = data[1];
    pr = data[2];
    d = if c == 0 then i else pr[(c-1)].x;
    t = if c == len(pr) then tr else pr[c].t };
module chen2D_JF_nm(react, prod, flux) = 
  [ " a inp tr + r tq + tq* r* tr* rinp a* ta*" 
  | " (  ~  (  + ( (  +  )  )   )   ~   )   . " ] + 
  [ " tr r " | ". ." ] + 
  [ " i + out + tr r + tq* r* tr* rout i* "
  | " ( +  ~  + (  ( +  .  )   )   ~   )  " ] + ofuels
    where {
      a = react[0].x;
      ta = react[0].t;
      [inp, rin] = flip(map(chen2D_I, tail(react)), 2);
      rinp = reverse(rin);
      [out, rou] = flip(map(chen2D_O, reverse(prod)), 2);
      rout = reverse(rou);
      r = flux.dr;
      tr = flux.tr;
      tq = flux.tq;
      i = flux.di;
      ofuels = map2(chen2D_fo, [i, tr, reverse(prod)], range(len(prod)+1))
    };

# find previous reactions with the same reactants in same order, return shared domains
function optimize(data, i) =
  if len(crn[i].reactants) > 0 and crn[i].reactants == r.reactants then previous[i] else 0 
  where {
    crn = data[0];
    previous = data[1];
    r = data[2] };

# remove every entry of l that does not contain data.
function trim(l) = 
  if l == [] then [] elseif l[0] == 0 then trim(tail(l)) else [l[0]] + trim(tail(l));

module delayedChoice(crn, i, previous) =
  if i >= len(crn) then
    empty
  elseif len(crn[i].reactants) == 0 then
    sum(map(infty, gates)) + delayedChoice(crn, i+1, previous + [newfl])
      where {
        f = fuel();
        newfl = [];
        flux = flux_domains();
        gates = [f] + chen2D_JF_nm([f], crn[i].products, flux) }
  else
    sum(map(infty, gates)) + delayedChoice(crn, i+1, previous + [newfl])
      where {
        r = crn[i];
        seen = trim(map2(optimize, [crn, previous, r], range(i)));
        flux = if len(seen) > 0 then seen[0] else flux_domains();
        newfl = if len(seen) > 0 then [] else flux;
        gates = chen2D_JF_nm(r.reactants, r.products, flux)
      };

module main(crn) = delayedChoice(crn, 0, [])
    where crn = irrev_reactions(crn);

