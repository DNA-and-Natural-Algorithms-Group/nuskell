#
# Niranjan Srinivas PhD Thesis: "Programming chemical kinetics: engineering
# dynamic reaction networks with DNA strand displacement.", Caltech (2015)
#
# Note: * Implements Figures from Chapter 3: Adventures in programming and
#         debugging molecular reaction networks:
#         - Figure 3.4: {B + A -> X + Y}
#         - Figure 3.5: {B -> X}
#         - Figure 3.6: {B -> }
#         - Figure 3.7: { -> X}
#
#       * Generalized on the DNA level for higher order reactions.
#
# Coded by Stefan Badelt (badelt@caltech.edu)
#

class formal(s) = "? f m s" | ". . . ."
  where {
    f = short();
    m = long();
    s = short() };

macro prodg(s) =
  [ "hy fy my sy + "
  | " (  (  .  . + ",
    "fy* hy*" | ") )",
    "hy fy" | ". ."]
  where {
   hy = long();
   fy = s.f;
   my = s.m;
   sy = s.s };

module srinivas_pgate(r, p, hflux) = 
  [ "hx fx mx sx + n + m fx* hx* sa*"
  | "(  (  .  .  + ~ + ~  )   )   . ",
    "fx help" | ".  ~"] # helper
  where {
    sa = r[len(r)-1].s;
    hx = hflux; 
    fx = p[0].f; 
    mx = p[0].m; 
    sx = p[0].s;
    [n, m, help] = flip(map(prodg, tail(p)), 3);
    m = reverse(m)
  };

module srinivas_pgate_p1(r, p, hflux) = 
  [ "hx fx mx sx + fx* hx* sa*"
  | "(  (  .  .  +  )   )   . "] 
  where {
    sa = r[len(r)-1].s;
    hx = hflux;
    fx = p[0].f; 
    mx = p[0].m; 
    sx = p[0].s 
  };

module flux(r,p) = if len(p) == 0 then [[]] 
  else if len(p) == 1 then [["h f" | ". ."]]
  where {h = long(); f = p[0].f }
  else [["h" | "."]]
  where {h = long(); f = p[0].f };

module srinivas_rgate(r, p) = 
  # can deal with all len(p)
  # can deal with len(r) = 0 and len(r) > 1
  [ "n + mr sr fl + sr* mr* m fr*"
  | "~ + (  (  ~  +  )   )  ~  . ", "l" | "~" ]
  where {
    fr = r[0].f ;
    mr = r[len(r)-1].m ;
    sr = r[len(r)-1].s ;
    [fl] = flux(r,p);
    [n,m,l] = flip(map2(reactg, r, range(len(r)-1)),3);
    #void = print(r,p);
    #void = print([n,m,l]);
    m = reverse(m) };

module srinivas_rgate_r0(r, p) = 
  [ "mr sr fl + sr* mr* fr*"
  | "(  (  ~  +  )   )   . ", "f m s " | ". . ." ]
  where {
    f = short();
    m = long();
    s = short();
    fr = f;
    mr = m;
    sr = s;
    [fl] = flux(r,p) };

macro reactg(r, i) = 
  [ "mb sb fa + " 
  | "(  (  (  + ", 
    "fa* sb* mb*" 
  | " )   )   ) ",
    "mb sb fa"
  | ".  .  . " ]
  where {
    mb = r[i].m ;
    sb = r[i].s ;
    fa = r[i+1].f };

# At some point do the genralized version
# module rxn(r) = infty(react) + infty(produce) 
#   where
#     [react, produce] = srinivas_gates(r.reactants, r.products) ;

module rxn(r) = sum(map(infty, react + produce))
  where {
    react = 
      if len(r.reactants) == 0 then 
        #-- needs a specific reactant-fuel
        srinivas_rgate_r0(r.reactants, r.products)
      else 
        #-- adjusts back and fuel to number of products
        srinivas_rgate(r.reactants, r.products);

    produce = 
      if len(r.products) == 0 then 
        []
      else if len(r.products) == 1 then 
        if len(r.reactants) == 0 then 
          srinivas_pgate_p1(react, r.products, hist) where hist = react[0].fl[0].h
        else 
          srinivas_pgate_p1(r.reactants, r.products, hist) where hist = react[0].fl[0].h
      else 
        if len(r.reactants) == 0 then 
          srinivas_pgate(react, r.products, hist) where hist = react[0].fl[0].h
        else 
          srinivas_pgate(r.reactants, r.products, hist) where hist = react[0].fl[0].h 
  };

module main(crn) = sum(map(rxn, crn)) 
  where 
    crn = irrev_reactions(crn)

# If you have only one product, you need to make the flux different!

