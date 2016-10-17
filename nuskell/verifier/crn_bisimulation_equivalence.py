# crn_bisimulation_equivalence.py written by Qing Dong <dong.qing.nju@gmail.com>

#from sets import Set
import itertools, time, basis_finder, string, copy
# python's collections.Counter is a data type for multisets,
#  such as states of a CRN or interpretations
from collections import Counter

permissive_depth = 8
midsearchDepth = 3
calcMidDepth = True
intr = {}
f = True
max_depth = 0
permissive_failure = [[[],[]],[]]
printing = True
debug = False

def printRxn(rxn):
    first = True
    for x in rxn[0]:
        if x[0] not in string.letters:
            xname = "i" + x
        else:
            xname = x
        if not first:
            print "+",
        else:
            first = False
        if rxn[0][x] > 1:
            print rxn[0][x],
        print xname,
    print "->",
    first = True
    for x in rxn[1]:
        if x[0] not in string.letters:
            xname = "i" + x
        else:
            xname = x
        if not first:
            print "+",
        else:
            first = False
        if rxn[1][x] > 1:
            print rxn[1][x],
        print xname,
    print

def output(intrp):
    for sp in intrp:
        print "   ",
        k = sp
        if k[0:4] == 'impl':
            k = k[4:]
        printRxn([{k: 1}, intrp[sp]])
    print

def solve(a):
# Find a non-negative and non-trivial integer solution x of the equation ax=0.
# Return [] when there is no such solution.
# (This is the algorithm from Contejean & Devie 1994.
    q = len(a[0])
    
    def multi(x, y):
        s = 0
        for i in range(len(x)):
            s = s + x[i] * y[i]
        return s
    
    def sub(x):
        s = []
        for i in range(len(a)):
            s.append(multi(a[i],x))
        return s
    
    def Min(b,t):
        if not b:
            return True
        else:
            for i in range(len(b)):
                r = True;
                for j in range(q):
                    r = r and (b[i][j] <= t[j])
                if r:
                    return False
            return True
        
    e = []
    for i in range(q):
        e.append([])
        for j in range(len(a)):
            e[i].append(a[j][i])     
    p = [[0 for i in range(q)]]
    frozen = [[False for i in range(q)]]
    zero = [0 for i in range(len(a))]
    zero1 = [0 for i in range(q)]
    b = []
    while p:
        t = p.pop()
        if sub(t) == zero and t != zero1:
            if t[q-1] == 1:
                return t    # just get the first solution, not all solutions (unlike C&D 1994).
            b.append(list(t))
            frozen.pop()
        else:
            f = frozen.pop()
            for i in range(q):
                if not f[i] and (multi(sub(t), e[i]) < 0) or t == zero1:
                    tmp = list(t)
                    tmp[i] += 1
                    if Min(b, tmp):
                        if i == q-1:
                            f[i] = True
                        p.append(tmp)
                        frozen.append(list(f))
                        f[i] = True
    return []

# multiset difference is collections.Counter's '-' operator
# note: NOT collections.Counter's 'subtract' function
# symmetric difference is two copies of that

def msleq(x,y):
    # True if (multisets) x <= y (vector comparison)
    for k in x:
        if x[k] > y[k]:
            return False

    return True

def mstimes(s, l):
    # return l*s for integer l, multiset s
    c = Counter()
    for k in s:
        c[k] = l * s[k]
    return c

def msdiv(s, l):
    # return s/l for integer l, multiset s if l divides s, otherwise False
    # l divides s if l divides s[k] for every key k
    c = Counter()
    for k in s:
        c[k] = s[k]/l
        if c[k] * l != s[k]:
            return False
    return c

def subsets(x):
# generates all subsets of multiset x
# Python's 'itertools.product' method works almost perfectly, and should
#  work without duplicates
    ks = x.keys()
    vs = map(lambda v: xrange(v+1), x.values())
    for prod in itertools.product(*vs):
        # calls to keys and values with no dictionary modifications in between
        #  should produce the keys and values in the same order,
        #  and product should respect that order (I think)
        yield Counter(dict(zip(ks, prod))) + Counter()

def enum(n, s, weights=None):
# partition multiset s into n ordered (possibly empty) parts.  
# e.g. enum(2, [a, b]) = [ [[],[a,b]], [[a],[b]], [[b],[a]], [[a,b],[]] ]
# if weights are given, enumerates all lists l of n multisets
#  such that s = sum(weights[i] * l[i])
    if weights is None:
        weights = [1] * n
    if n == 0:
        yield []
        return
    if len(weights) < n or weights[0] < 0:
        raise Exception
    elif weights[0] == 0:
        for j in enum(n-1, s, weights[1:]):
            yield [Counter()] + j
        return
    if n == 1:
        sdivw = msdiv(s, weights[0])
        if sdivw is not False: yield [sdivw]
        return
    x = subsets(s)
    for i in x:
        ss = mstimes(i, weights[0])
        if not msleq(ss, s):
            continue
        for j in enum(n-1, s-ss):
            yield [i] + j

def interpret(s, intrp):
    # return interpretation of s according to intrp
    ss = s.copy()
    ks = ss.keys()
    for k in ks:
        if k in intrp:
            v = ss.pop(k)
            for i in range(v):
                ss += intrp[k]

    return ss

def interleq(x, y, intrp):
    # True if m(x) <= m(y) with m given by interpretation intrp
    return msleq(interpret(x,intrp),interpret(y,intrp))
       
def subst(crn, intrp):
# Substitute implementation species for formal species in CRN according to interpretation.
    return [[interpret(j,intrp) for j in rxn] for rxn in crn]

def checkT(T):
# checks the table to see if there is a whole row or a whole (non-trivial) column that's all false, so we have to roll back.
# Returns False if we have to roll back.  (i.e. our partial interpretation fails the delimiting condition)
    for i in range(len(T)):
        if True not in T[i]:
            return False
    for i in range(len(T[0])-1):
        r = False
        for j in range(len(T)):
            r = r or T[j][i]
        if not r:
            return False
    return True

def update(fcrn, icrn, fs):
# the um, er, um, completely recalculates the table from scratch.
# assumes subst has already been applied to implementation icrn
# This should be logically equivalent to the UpdateTable in the MS thesis for compiled DNA reactions (we think).
# WARNING:
# If an implementation CRN has directly catalytic species, the code below may fail (though thesis psuedocode is correct).
# E.g.  i3 + i7 --> i12 + i7 + i8
    m = len(icrn)
    r = []
    for i in range(len(icrn)):
        rr = []
        for j in range(len(fcrn)):
            t1 = fcrn[j][0] - icrn[i][0]
            t2 = icrn[i][0] - fcrn[j][0]
            t3 = fcrn[j][1] - icrn[i][1]
            t4 = icrn[i][1] - fcrn[j][1]
            if set(t2).isdisjoint(set(fs)) and set(t4).isdisjoint(set(fs)):
                if t1.keys() == []:
                    if t3.keys() == []:
                        rr.append(True)
                    else:
                        if t4.keys() == []:
                            rr.append(False)
                        else:
                            rr.append(True)
                else:
                    if t2.keys() == []:
                        rr.append(False)
                    else:
                        if t3.keys() == []:
                            rr.append(True)
                        else:
                            if t4.keys() == []:
                                rr.append(False)
                            else:
                                rr.append(True)
            else:
                rr.append(False)
        t1 = icrn[i][0] - icrn[i][1]
        t2 = icrn[i][1] - icrn[i][0]
        if (set(t1).isdisjoint(set(fs)) or not set(t2).issubset(set(fs))) and (set(t2).isdisjoint(set(fs)) or not set(t1).issubset(set(fs))):
            rr.append(True)
        else:
            rr.append(False)
        r.append(rr)
    return r

def perm(fcrn, icrn, fs, intrp, permcheck):
# check the permissive condition, using the global original formal crn and original implementation crn (without substitutions).
    global intr, f, permissive_depth, max_depth, permissive_failure
    global midsearchDepth, printing, calcMidDepth
    tr = []
    fr = []
    hasht = set([])
    f = True
    crn2 = subst(icrn, intrp)
    T = update(fcrn, crn2, fs)
    nulls = [k for k in intrp if not intrp[k].keys()] # null species
    now = time.clock()

    if printing:
        print "Testing permissive condition"
        print "Formal CRN:", fcrn
        print "Implementation CRN:", icrn
        print "Interpretation:", intrp
        print

    def cnstr(s):
        # given formal state s, generate minimal set of impl. states which
        #  interpret to a superset of s
        # = concatenation of any x such that s[0] in m(x) with
        #  all minimal implementations of s - m(x)
        if s.keys() == []:
            yield Counter()
        else:
            s1 = s.keys()[0]
            for k in intrp:
                if s1 in intrp[k]:
                    if (s - intrp[k])[s1] >= s[s1]:
                        print s, s1, k
                        print intrp
                        print intrp[k]
                        assert False
                    for out in cnstr(s - intrp[k]):
                        yield Counter({k:1}) + out

    def search(s, d):
        # s is the implementation state, d is the depth.
        # try to find (via trivial reactions) a state in which a reaction in fr can fire.
        # fr is a particular formal reaction along with all implementation reactions that interpret to it. 
        # tr is the list of all trivial reactions in the implementation.
        # fail if depth d of trivial reaction steps is exceeded without firing a reaction in fr.
        global f, printing
        if d > permissive_depth:
            f = False
            return False
        if list(s.elements()).sort() in hasht:
            return False
        else:
            hasht.add(list(s.elements()).sort()) 
        for i in fr:
            if (i[0] - s).keys() == []:
                return True
        for i in tr:
            if (i[0] - s).keys() == []:
                t = (s - i[0]) + i[1]
                if search(t, d+1):
                    return True
        return False

    def midsearch(start, goal, pickup, ignore, formal, k):
        global printing
        # search for a path from start to goal (states of non-null species)
        #  which produces at least pickup null species
        #  assuming it already has infinite copies of everything in ignore
        #  of length at most 2^k
        # if goal is None, the goal is any reaction in fr
        if goal is not None:
            if ignore.issuperset(goal - start):
                return True
            if not interleq(goal, start, intrp):
                return False

        if k == 0:
            if goal is None:
                for rx in fr:
                    if ignore.issuperset(rx[0] - start):
                        return True
                return False
            
            for rx in tr:
                if ignore.issuperset(rx[0] - start):
                    # every element of rx[0] - start (multiset)
                    #  is also an element of ignore (set)
                    # e.g. true on rx[0] - start = {|a,a,b|}, ignore = {a,b,c}
                    if msleq(goal, (start - rx[0]) + rx[1]):
                        return True

        else:
            if midsearch(start,goal,pickup,ignore,formal,k-1):
                return True
            for part in enum(2,pickup):
                for mid in cnstr(formal):
                    if midsearch(start,mid,part[0],ignore,formal,k-1) \
                       and ((not interleq(start, mid, intrp)) or
                            midsearch(mid,goal,part[1],ignore,formal,k-1)):
                        return True

        return False

    def loopsearch(start, formal):
        global f, printing
        if printing:
            print "Loop-searching from", start
            print " to any of", fr
        # search for a path from start to any reaction in fr
        rounds = midsearchDepth
        if calcMidDepth:
            nequiv = 0
            rounds = 0
            roundequiv = 1
            for point in cnstr(formal):
                nequiv += 1
                if nequiv > roundequiv:
                    rounds += 1
                    roundequiv *= 2

        if printing:
            print " Will search", nequiv, "states for", rounds, "rounds"
        for part in map(Counter.keys,enum(len(nulls) + 1,nulls)):
            if any([part[i] != [] and part[i+1] == [] \
                    for i in range(len(part) - 2)]):
                continue # avoid redundancy

            pickups = filter(None,part[:-1])
            if printing:
                print " Using null species partition:", pickups
            check1 = True
            place = start
            ignore = []
            for pickup in pickups:
                check2 = False
                for base in cnstr(formal):
                    if midsearch(place,base,[],ignore,formal,rounds):
                        if not interleq(place,base,intrp):
                            return True
                        elif midsearch(base,base,pickup,ignore,formal,rounds):
                            check2 = True
                            place = base
                            break

                if not check2:
                    check1 = False
                    break

                ignore.extend(pickup)

            if check1 and midsearch(place,None,[],ignore,formal,rounds):
                if printing:
                    print " Search success."
                return True

        f = False
        if printing:
            print " Search failure."
        return False

    def allsearch(formal):
        # check at once whether every implementation of state "formal"
        #  can implement the given formal reaction
        #  by storing for each state the list of states it can reach w/o
        #  null species (except those producible by loops)
        points = map(lambda x: [x,set([]),[]], cnstr(formal))
        # points[i][0] is the ith implementation state
        # points[i][1] is all null species loopable from that state
        # points[i][2][j] is True if points[j][0] is known to be reachable
        #  (from state i without initial null species)
        # exception: points[i][2] is True if a formal reaction is known to
        #  be reachable from state i
        l = len(points) # let's not recalculate this every time...
        rngl = range(l) # same here...
        for i in rngl:
            points[i][2] = l*[False]

        changed = True
        while changed:
            changed = False
            for i in rngl:
                if points[i][2] is not True:
                    for rx in fr:
                        if points[i][1].issuperset(rx[0] - points[i][0]):
                            points[i][2] = True
                            changed = True
                            break

                    if points[i][2] is True:
                        continue

                    for rx in tr:
                        if points[i][1].issuperset(rx[0] - points[i][0]):
                            left = points[i][0] - rx[0]
                            after = left + rx[1]
                            for j in rngl:
                                if msleq(points[j][0],after):
                                    if points[j][2] is True:
                                        points[i][2] = True
                                        changed = True
                                        break

                                    if points[j][2][i]:
                                        s = points[j][1].union(
                                            [x for x in left if x in nulls])
                                        if not s <= points[i][1]:
                                            points[i][1] |= s
                                            changed = True

                                    if not points[i][2][j]:
                                        points[i][2][j] = True
                                        changed = True

                                    for k in rngl:
                                        if (not points[i][2][k]) \
                                           and points[j][2][k]:
                                            points[i][2][k] = True
                                            changed = True

                        if points[i][2] is True:
                            break

        return all([p[2] is True for p in points])
    
    n = 0
    # build tr just once
    for i in T:
        if i.pop():
            tr.append(icrn[n])
        n += 1
    if printing: print "Reactions to check:", len(fcrn)
    for i in range(len(fcrn)):
        # build fr for this formal reaction
        if printing: print "Checking reaction", fcrn[i]
        fr = []
        for j in range(len(icrn)):
            if T[j][i]:
                fr.append(icrn[j])

        if printing:
            lstates = list(cnstr(fcrn[i][0]))
            print "Number of states:", len(lstates)
            print lstates
        ist = cnstr(fcrn[i][0])

        if permcheck == "whole":
            if not allsearch(fcrn[i][0]):
                if printing:
                    watch = time.clock()
                    print "Permissive test failed in time", watch - now
                if max_depth == -2:
                    return False
                if f == True:   # global variable saying whether permissive has been tried or not
                    max_depth = -1
                else:
                    max_depth = -2
                intr = intrp.copy()
                permissive_failure[0] = fcrn[i]
                permissive_failure[1] = ["Somewhere"]
                return False

            continue
        
        # At this point, "ist" contains an exhaustive and sufficient list of possible initial implementation states 
        # in which the current formal reaction #i must be able to fire (via a series of trivial reactions),
        # Note that we will only want to test states in "ist" that interpret to a state in which #i can fire.
        
        tested = []  # to avoid testing a superset of some state that's already been tested
        spaceefficient = True # ... but only if we have space to store them
        for j in ist:
            tmp = interpret(j,intrp)
            if msleq(fcrn[i][0], tmp):  # only test if reactants j interpret to allow #i to fire
                t = False
                for k in tested:
                    if msleq(k, j):
                        t = True
                        if printing: print "State", j, "is a superset of", k
                        break
                if t:
                    continue
                hasht = set([])
                found = ((permcheck=="pspace") and loopsearch(j,fcrn[i][0])) \
                        or ((not permcheck) and search(j,0))
                if not found: # didn't work, but keep track of stuff for user output
                    if printing:
                        watch = time.clock()
                        print "Permissive test failed in time", watch - now
                    printing = False
                    if max_depth == -2:
                        return False
                    if f == True:   # global variable saying whether permissive has been tried or not
                        max_depth = -1
                    else:
                        max_depth = -2
                    intr = intrp.copy()
                    permissive_failure[0] = fcrn[i]
                    permissive_failure[1] = j
                    return False
                elif not spaceefficient:
                    tested.append(j)
            elif printing:
                print "Not testing", j, "with interpretation", tmp
    intr = intrp.copy()
    if printing:
        watch = time.clock()
        print "Permissive test succeeded in time", watch - now
    return True

def equations(fcrn, icrn, fs, intrp, permcheck):
    # All unknown implementation reactions (i.e. those with some unknown species) must be trivial.
    # Build the matrix for the the "solve" function, to see whether the interpretation can be completed as required.
    # Also note that "unknow" is not used; the unknown species are recalculated here
    #  because "unknow" contains only those implementation reactions that have not been solved by the row search,
    #  but other rows (i.e. implementation reactions) may be implicitly solved by the partial interpretation.
    unknown = []
    ust = set([])
    sicrn = subst(icrn, intrp)
    for i in range(len(sicrn)):
        tmp = set(sicrn[i][0])-set(fs)
        tmp |= set(sicrn[i][1])-set(fs)
        ust |= tmp
        if tmp != set([]):
            unknown.append(i)  # "unknown" is the list of implementation reactions with an unknown species
    us = list(ust)  # all species that remain unknown in current (partial) interpretation
    
    # check atomic condition
    atoms = set()
    for k in intrp:
        sps = intrp[k].keys()
        if len(sps) == 1 and intrp[k][sps[0]] == 1:
            atoms.add(sps[0])

    atomsleft = list(set(fs) - atoms)
    l = len(atomsleft)
    for assign in itertools.permutations(us, l): # works even if l == 0
        # each assign is a tuple of implementation species to be interpreted
        #  as exactly one formal species, matching the order of atomsleft
        # if l == 0 then it.permutations(us, l) == [()], a list with
        #  element which is the zero-length tuple,
        #  so this loop will run exactly once with no change to itmp
        itmp = dict(intrp)
        for i in range(l):
            assert assign[i] not in itmp
            itmp[assign[i]] = Counter({atomsleft[i]: 1})

        sicrntmp = subst(icrn, itmp)
        T = update(fcrn, sicrntmp, fs)
        if (not checkT(T)) or any([not T[i][-1] for i in unknown]):
            # either the table is bad, or a reaction we are assuming
            #  is trivial can no longer be trivial
            continue

        ustmp = [sp for sp in us if sp not in assign]
        if not ustmp:
            out = perm(fcrn, icrn, fs, itmp, permcheck)
            if out:
                return out
            else:
                continue

        n = 0
        a = []
        for i in unknown:
            a.append([])
            for j in ustmp:
                a[n].append(sicrntmp[i][0][j]-sicrntmp[i][1][j])
                a[n].append(0)

            n += 1

        for isp in ustmp:
            itmp[isp] = Counter()

        check = True
        for fsp in fs:
            n = 0
            for j in unknown:
                a[n].pop()
                a[n].append(sicrn[j][0][fsp]-sicrn[j][1][fsp])
                n += 1

            s = solve(a)
            if s == []:
                check = False
                break
            else:
                for isp in us:
                    itmp[isp][fsp] = s[j]

        if check:
            for isp in ustmp:
                itmp[isp] = itmp[isp] + Counter()
            out = perm(fcrn, icrn, fs, itmp, permcheck)
            if out:
                return out
            else:
                continue

    return False

def searchr(fcrn, icrn, fs, unknown, intrp, d, permcheck, nontriv=False):
# Row search.  I.e. make sure every implementation reaction can interpret to a formal reaction (or be trivial).
    global max_depth, intr
    if unknown == []:
        return perm(fcrn, icrn, fs, intrp, permcheck)
    sicrn = subst(icrn, intrp)
    T = update(fcrn, sicrn, fs)
    if not checkT(T):
        return False
    if max_depth >= 0 and d > max_depth:
        intr = copy.deepcopy(intrp)
        max_depth = d
    min = len(fcrn)+1
    k = -1  # next row reaction we will solve
    nt = 0  # number of possibly trivial reactions according to table
    for i in unknown:
        tmp = T[i].count(True)
        if T[i][-1] == True:
            nt += 1
            tmp -= 1
        if tmp < min and tmp > 0:
            min = tmp
            k = i
    if nt == len(unknown) and not nontriv:  # all unsearched rows could be trivial -- try it!
        if equations(fcrn, icrn, fs, intrp, permcheck):
            return True
        else:
            # if we just tried equations and it didn't work, then we
            #  shouldn't try again unless something changes
            nontriv = True
    if k < 0:
        return False
    untmp = list(unknown)
    untmp.remove(k)
    if T[k][-1] == True:  # if implementation reaction #k can be trivial, leave it that way and try more rows
        if searchr(fcrn, icrn, fs, untmp, intrp, d, permcheck, nontriv):
            return True
    n = 0
    for c in range(len(fcrn)): # try to match implementation reaction #k with some formal reaction
        if T[k][c]:        
            ul = sicrn[k][0] - fcrn[c][0]
            kl = ul.keys()
            nl = len(kl)
            sl = fcrn[c][0] - sicrn[k][0]
            tmpl = enum(nl, sl, ul.values())
            ur = sicrn[k][1] - fcrn[c][1]
            kr = ur.keys()
            nr = len(kr)
            sr = fcrn[c][1] - sicrn[k][1]
            tmpr = enum(nr, sr, ur.values())
            for i in tmpl:
                intrpleft = dict(zip(kl, i))
                for j in tmpr:
                    intrpright = dict(zip(kr, j))
                    checkCompatible = True
                    for key in intrpleft:
                        if key in intrpright and \
                           any([intrpright[key][fsp] != intrpleft[key][fsp]
                                for fsp in fs]):
                            checkCompatible = False
                            break

                    if not checkCompatible:
                        continue

                    itmp = copy.deepcopy(intrp)
                    itmp.update(intrpleft)
                    itmp.update(intrpright)
                    if searchr(fcrn, icrn, fs, untmp, itmp, d+1, permcheck):
                        return True
    return False

def searchc(fcrn, icrn, fs, unknown, intrp, d, permcheck):
    # Search column.  I.e. make sure every formal reaction can be implemented.
    global max_depth, intr
    sicrn = subst(icrn, intrp)
    T = update(fcrn, sicrn, fs)
    if not checkT(T):
        return False
    if max_depth >= 0 and d > max_depth:
        intr = copy.deepcopy(intrp)
        max_depth = d
    min = len(icrn)+1
    c = -1  # this will be the next column to solve, if possible
    for i in unknown:
        tmp = 0
        for j in range(len(icrn)):
            if T[j][i]:
                tmp += 1
        if tmp < min:
            min = tmp
            c = i
    if debug:
        print "Solving column", c
    if c < 0:  # done with column search.  transition to row search!
        untmp = []
        for i in range(len(icrn)):
            if not (set(sicrn[i][0])-set(fs) == set([]) and \
                    set(sicrn[i][1])-set(fs) == set([])):
                untmp.append(i)
        return searchr(fcrn, icrn, fs, untmp, intrp, d, permcheck)
    else:
        untmp = list(unknown)
        untmp.remove(c)
        n = 0
        for k in range(len(icrn)):
            if T[k][c]:
                ul = sicrn[k][0] - fcrn[c][0]
                kl = ul.keys()
                nl = len(kl)
                sl = fcrn[c][0] - sicrn[k][0]
                tmpl = enum(nl, sl, ul.values())
                ur = sicrn[k][1] - fcrn[c][1]
                kr = ur.keys()
                nr = len(kr)
                sr = fcrn[c][1] - sicrn[k][1]
                tmpr = enum(nr, sr, ur.values())
                for i in tmpl:
                    intrpleft = dict(zip(kl, i))
                    for j in tmpr:
                        intrpright = dict(zip(kr, j))

                        checkCompatible = True
                        for key in intrpleft:
                            if key in intrpright and \
                               any([intrpleft[key][fsp] != intrpright[key][fsp]
                                    for fsp in fs]):
                                checkCompatible = False
                                break

                        if not checkCompatible:
                            continue

                        itmp = copy.deepcopy(intrp)
                        itmp.update(intrpleft)
                        itmp.update(intrpright)
                        if debug:
                            print itmp
                        if searchc(fcrn, icrn, fs, untmp, itmp, d+1, permcheck):
                            return True
    return False

def test(c1, c2, verbose = True, intrp = None, permcheck=False):
    global printing, debug
    global intr, max_depth, permissive_failure
    printing = printing and (verbose or debug)
    (fcrn, _) = c1
    fcrn = [[Counter(part) for part in rxn] for rxn in fcrn]
    (icrn, fs) = c2
    icrn = [[Counter(part) for part in rxn] for rxn in icrn]
    icrn = copy.deepcopy(icrn)
    for rxn in icrn:
        for k in rxn[0]:
            if k in fs:
                v = rxn[0].pop(k)
                rxn[0]['impl'+k] = v
        for k in rxn[1]:
            if k in fs:
                v = rxn[1].pop(k)
                rxn[1]['impl'+k] = v
    # maybe list->Counter conversion is temporary?

    if intrp is None: # default 1: no interpretation information
        intrp = {}
    elif intrp is True: # default 2: each fsp has a canonical implementation
        intrp = {fsp: Counter({fsp: 1}) for fsp in fs
                 if any([fsp in rxn[0] or fsp in rxn[1] for rxn in icrn])}

    if verbose:
        print "Original CRN:"
        for rxn in fcrn:
            print "   ",
            printRxn(rxn)
            print
        if icrn == []:
            print "Compiled CRN is empty"
            print
            return fcrn == []
        print "Compiled CRN:"
        for rxn in icrn:
            print "   ",
            printRxn(rxn)
        print
    elif icrn == []:
        return fcrn == []
    unknown = [i for i in range(len(fcrn))]
    if searchc(fcrn, icrn, fs, unknown, intrp, 0, permcheck):
        if verbose:
            print "Valid interpretation :"
            output(intr)
        return True
    else:
        if verbose:
            if max_depth >= 0:
                print "Delimiting condition cannot be satisfied."
                if max_depth >= len(fcrn):
                    print "There is implementation reaction not in formal CRN."
                else:
                    print "There is formal reaction not implemented."
                print "Max search depth reached :", max_depth
                print "with interpretation :"
                output(intr)
            else:
                print "Fail in permissive test with interpretation :"
                output(intr)
                print "at formal reaction :",
                printRxn(permissive_failure[0])
                print "on implementation status :", permissive_failure[1]
                if max_depth == -2:
                    print "with max trivial reaction chain length", permissive_depth, "reached."
            print
        return False

if __name__ == "__main__":
    # The name of the program
    program_name = "test"
    fcrn = [[['a'],['b']]]
    fs1 = []
    cs1 = []
    icrn = [[['a1'],['b1']],[['x'],['a1']],[['x'],['b1']],[['y'],['b1']],[['y'],['a1']],[['x'],['a0']],[['a0'],['a1']]]
    fs2 = ['a','b']
    cs2 = []
    v = test((fcrn, fs1), (icrn, fs2))
    if v:
        print "verify: compilation was correct."
    else:
        print "verify: compilation was incorrect."
