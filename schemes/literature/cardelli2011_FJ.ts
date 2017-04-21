#
# Luca Cardelli's translation scheme from "Strand Algebras for DNA
# Computing", Natural Computing, 10: 407-428, 2011. Can be found at
# http://lucacardelli.name/
#
#   Note: * This implements the `fork' and `join' gates from the paper,
#           which are pointed out to be problematic in the paper. 
#           Schemes are shown in Figs 4-8.
#
#         * In addition to the reactions describe in the paper, the 
#           spontaneous reaction ' -> A' is implemented as 'f -> A' 
#           where f is a fuel strand.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com) 
# modified by Stefan Badelt (badelt@caltech.edu)
#

class formal(s) = "? t b"
                | ". . ."
    where {
        t = short();
        b = long() };

class signal() = "? t b"
               | ". . ."
    where {
        t = short();
        b = long() };

# Figure 6: 2-way Fork
class forkgate(x, y, z)
    = [ "xb yt yb + a zt zb + zt* a* yt* xb* xt*"
      | " (  (  . + (  (  . + )   )  )   )   .",
        "yt a zt"
      | ".  .  ."]
    where {
        xt = x.t;
        xb = x.b;
        yt = y.t;
        yb = y.b;
        zt = z.t;
        zb = z.b;
        a = long() };

module fork2(x, y, z) = infty(gb) + infty(gt)
    where
        [gb, gt] = forkgate(x, y, z);

# r = single reactant, p = list
module fork(r, p) =
    if len(p) == 2 then
        fork2(r, p[0], p[1])
    else
        fork2(i, p[0], p[1]) + 
        fork(r, tail(tail(p)) + [i])
        where
            i = signal();

# Figure 7: 2-way Join
class joingate(x, y, z)
    = [ "xb yt + yb a + b zt zb + zt* b* a* yb* yt* xb* xt*"
      | "(  (  + (  ( + ( (  .  +  )  )  )   )   )   )   .",
        "a b zt"
      | ". .  .",
        "xb yt"       # r1
      | ".  .",
        "yb a"        # r2
      | ".  ."]
    where {
        xt = x.t;
        xb = x.b;
        yt = y.t;
        yb = y.b;
        zt = z.t;
        zb = z.b;
        a = short();
        b = long() };

# Figure 8: join gate - garbage collection
class joincleanup(r1, r2)
    = [ "xb + yt c + d yb + a* yb* d* c* yt* xb* " 
      | " ( +  ( ( + ( (  + .   )  )  )   )   )  ",
        " c d " 
      | " . . ",
        " yb + yb* d* " 
      | " (  +  )  . ",
        " c + c* yt* "
      | " ( +  )  .  "]
    where {
        xb = r1.xb;
        yt = r1.yt;
        c  = long();
        d = short();
        yb = r2.yb;
        a = r2.a };

module join2(x, y, z) = sum(map(infty, [gb, gt, r1, r2] + joincleanup(r1, r2)))
  where 
    [gb, gt, r1, r2] = joingate(x, y, z);

# r = list, p = single product
module join(r, p) =
    if len(r) == 2 then
        join2(r[0], r[1], p)
    else
        # need a reversible reaction!
        fork2(i, r[0], r[1]) + join2(r[0], r[1], i) +
        join(tail(tail(r)) + [i], p)
        where
            i = signal();

class transgate(x, y)
    = [ "xb yt yb + a + a* yt* xb* xt*"
      | " (  (  . + ( + )  )   )   .",
        "yt a"
      | " . ."]
    where {
        xt = x.t;
        xb = x.b;
        yt = y.t;
        yb = y.b;
        a = long() };

# Figure 5: one input -> one output
module transducer(x, y) = infty(gb) + infty(gt)
    where [gb, gt] = transgate(x, y);

class anhlgate(x)
    = "xb + xb* xt*"
    | " ( +  )   ."
    where {
        xt = x.t;
        xb = x.b };

# Figure 4: one input -> no output
module annihilator(x) = infty(g)
    where
        g = anhlgate(x);

module joinfork(r, p) =
    if len(p) == 0 then
      if len(r) == 1 then
        annihilator(r[0])
      else
        join(r, i) + annihilator(i)
        where i = signal()

    elseif len(r) == 0 then
      # abort('reaction type not implemented')          
      if len(p) == 1 then
        transducer(i, p[0]) + infty(i)
        where i = signal()
      else
        fork(i, p) + infty(i)
        where i = signal()

    elseif len(r) == 1 and len(p) == 1 then
        transducer(r[0], p[0])
    elseif len(r) == 1 then
        fork(r[0], p)
    elseif len(p) == 1 then
        join(r, p[0])
    else
        join(r, i) + fork(i, p)
        where
            i = signal();

module reaction(r) = joinfork(r.reactants, r.products);

module main(crn) = sum(map(reaction, crn))
  where crn = irrev_reactions(crn);
