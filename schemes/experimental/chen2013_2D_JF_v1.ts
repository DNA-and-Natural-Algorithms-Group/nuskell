#
# Yuan-Jyue Chen's translation scheme from "Programmable chemical controllers
#   made from DNA", Nature, 2013.  #
#
#   Note: * Implements supplementary Figure S7.6
#         * includes the delayed choice optimization suggested by the authors
#
# Coded by Stefan Badelt (badelt@caltech.edu)
#

# Switch delayedChoice optimization ON or OFF
global ON = 1;
global OFF = 0;
global do_opti = ON;

global void = if do_opti then print("Optimization: ON") else print ("Optimization: OFF");

class formal(s) = "t x" | ". ."
    where { t = short(); x = long() };

#class fuel() = "t x" | ". ."
#    where { t = short(); x = long() };

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

macro chen2D_fi(data, c) = "d t"|". ."
  where {
    a = data[0];
    tr = data[1];
    pr = data[2];
    d = if c == 0 then a else pr[(c-1)].x;
    t = if c == len(pr) then tr else pr[c].t };

macro chen2D_fo(data, c) = "d t"|". ."
  where {
    i = data[0];
    tr = data[1];
    pr = data[2];
    d = if c == 0 then i else pr[(c-1)].x;
    t = if c == len(pr) then tr else pr[c].t };

module chen2D_JF_nm_del(react, prod, flux) = 
  [ " a inp tr + i r tq + tq* r* tr* rinp a* ta*" 
  | " (  ~  (  + . ( (  +  )  )   )   ~   )   . " ] + 
  [ " tr r " | ". ." ] + ifuels +
  [ " i + out + tr i + r + tq* r* i* tr* rout i* "
  | " ( +  ~  + (  ( + ( +  .  )  )   )   ~   )  " ] +
    ofuels
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
      ifuels = map2(chen2D_fi, [a, tr, react], range(len(react)+1));
      ofuels = map2(chen2D_fo, [i, tr, reverse(prod)], range(len(prod)+1))
    };

module chen2D_JF_nm(react, prod, flux) = 
  [ " a inp tr + i + r tq + tq* r* i* tr* rinp a* ta*" 
  | " (  ~  (  + ( + ( (  +  )  )  )   )   ~   )   . " ] + 
  [ " tr i r " | ". . ." ] + ifuels +
  [ " i + out + tr r + tq* r* tr* rout i* "
  | " ( +  ~  + (  ( +  .  )   )   ~   )  " ] +
    ofuels
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
      ifuels = map2(chen2D_fi, [a, tr, react], range(len(react)+1));
      ofuels = map2(chen2D_fo, [i, tr, reverse(prod)], range(len(prod)+1))
    };

# find previous reactions with the same reactants in same order, return shared domains
function optimize(data, i) =
  if len(crn[i].reactants) > 0 and crn[i].reactants == r.reactants then previous[i] else 0 
  where {
    #void = print('optimizing');
    crn = data[0];
    previous = data[1];
    r = data[2] };

# remove ever entry of l that does not contain data.
function trim(l) = 
  if l == [] then [] elseif l[0] == 0 then trim(tail(l)) else [l[0]] + trim(tail(l));

module delayedChoice(crn, i, previous) =
  if i >= len(crn) then
    empty
  else
    sum(map(infty, gates)) + delayedChoice(crn, i+1, previous + [newfl])
      where {
        r = crn[i];
        seen = if do_opti then trim(map2(optimize, [crn, previous, r], range(i))) else [];
        #void = print('s', seen);
        flux = if len(seen) > 0 then seen[0] else flux_domains();
        #void = print('f', flux);
        newfl = if len(seen) > 0 then [] else flux;
        #void = print('nf', newfl);
        gates = chen2D_JF_nm(r.reactants, r.products, flux)
        #void = print('g', gates)
      };

module main(crn) = delayedChoice(crn, 0, [])
    where crn = irrev_reactions(crn);

