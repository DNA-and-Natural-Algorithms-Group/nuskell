#
# Luca Cardelli's translation scheme from "Strand Algebras for DNA
# Computing", Natural Computing, 10: 407-428, 2011. Can be found at
# http://lucacardelli.name/
#
#   Note : this implements the `fork' and `join' gates from the paper,
#          which are pointed out to be problematic in the paper. This
#          implements the schemes shown in Figs 4-7.
#
# Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
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

class joingate(x, y, z)
    = [ "xb yt + yb a + b zt zb + zt* b* a* yb* yt* xb* xt*"
      | "(  (  + (  ( + ( (  .  +  )  )  )   )   )   )   .",
        "a b zt"
      | ". .  .",
        "xb yt"
      | ".  .",
        "yb a"
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

module join2(x, y, z) = infty(gb) + infty(gt) + infty(r1) + infty(r2)
    where
        [gb, gt, r1, r2] = joingate(x, y, z);

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

module transducer(x, y) = infty(gb) + infty(gt)
    where
        [gb, gt] = transgate(x, y);

class anhlgate(x)
    = "xb + xb* xt*"
    | " ( +  )   ."
    where {
        xt = x.t;
        xb = x.b };

module annihilator(x) = infty(g)
    where
        g = anhlgate(x);

module fork(r, p) =
    if len(p) == 2 then
        fork2(r, p[0], p[1])
    else
        fork2(i, p[0], p[1]) + 
        fork(r, tail(tail(p)) + [i])
        where
            i = signal();

module join(r, p) =
    if len(r) == 2 then
        join2(r[0], r[1], p)
    else
        fork2(i, r[0], r[1]) + join2(r[0], r[1], i) +
        join(tail(tail(r)) + [i], p)
        where
            i = signal();

module joinfork(r, p) =
    if len(p) == 0 then
        if len(r) == 1 then
            annihilator(r[0])
        else
            fork(r, i) + annihilator(i)
            where
                i = signal()
    elseif len(r) == 0 then
      abort('reaction type not implemented')          

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
