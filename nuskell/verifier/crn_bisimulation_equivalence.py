# crn_bisimulation_equivalence.py written by Qing Dong <dong.qing.nju@gmail.com>

#from sets import Set
import itertools, time, string, copy
import basis_finder

permissive_depth = 8
fcrn = []
icrn = []
intr = [[],[]]
f = True
max_depth = 0
permissive_failure = [[[],[]],[]]

def printRxn(rxn):
    first = True
    for x in rxn[0]:
        if x[0] not in string.letters:
            x = "i" + x
        if not first:
            print "+",
        else:
            first = False
        print x,
    print "->",
    first = True
    for x in rxn[1]:
        if x[0] not in string.letters:
            x = "i" + x
        if not first:
            print "+",
        else:
            first = False
        print x,
    print

def output(inter):
    tmp = []
    for i in range(len(inter[0])):
        tmp.append([[inter[0][i]],inter[1][i]])
    tmp.sort(key = lambda r: -len(r[1]))
    for i in range(len(tmp)):
        print "   ",
        printRxn(tmp[i])
    print

def solve(a):
    """
    Find a non-negative and non-trivial integer solution x of the equation
    ax=0. Return [] when there is no such solution. This is the algorithm
    from Contejean & Devie 1994.

    :param a: A is a matrix used for ax=0

    :return: x, such that ax=0
    """
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

def diff(x,y):
# things in x but not y (for multisets)
    r = list(x)
    for i in y:
        if i in r:
            r.remove(i) 
    return r

def sdiff(mset1, mset2):
# partitioned symmetric difference for multisets
# ( things in mset1 but not mset2, things in mset2 but not mset1 )
    r1 = list(mset1)
    r2 = []
    for i in range(len(mset2)):
        if mset2[i] in r1:
            r1.remove(mset2[i])
        else:
            r2.append(mset2[i])
    return (r1, r2)

def subsets(x):
# generates all subsets of multiset x (will include duplicates if x has duplicates)
# (would be better, and OK with rest of code, to remove the duplicates)
    r = [[y for (j, y) in enumerate(x) if (i >> j) & 1] for i in range(2**len(x))]
    return r

def enum(n,s):
# partition multiset s into n ordered (possibly empty) parts.  
# (doesn't sort, and with multisets it will have duplicates.  would be better to remove duplicates, e.g. in subsets().)
# e.g. enum(2, [a, b]) = [ [[],[a,b]], [[a],[b]], [[b],[a]], [[a,b],[]] ]
    if n == 0:
        return[[]]
    if n == 1:
        return [[s]]
    x = subsets(s)
    r = []
    for i in x:
        for j in enum(n-1, diff(s,i)):
            j.append(i)
            if not j in r:
                r.append(j)
    return r

def subst(crn, uslist, fslist):
# Substitute implementation species for formal species in CRN according to interpretation.
# uslist is a list of unknown species (i.e. not formal species yet in CRN, which has already been partially substituted).
# fslist is a list of multisets of formal species.  fslist[i] is used to replace species uslist[i] in crn.
    r = []
    for i in crn:
        t1 = []
        for j in i:
            t2 = []
            for k in j:
                if k in uslist:
                    t2.extend(fslist[uslist.index(k)])
                else:
                    t2.append(k)
            t1.append(t2)
        r.append(t1)
    return r
       
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

def update(crn1, crn2, fs):
# the um, er, um, completely recalculates the table from scratch.
# This should be logically equivalent to the UpdateTable in the MS thesis for compiled DNA reactions (we think).
# WARNING:
# If an implementation CRN has directly catalytic species, the code below may fail (though thesis psuedocode is correct).
# E.g.  i3 + i7 --> i12 + i7 + i8
    m = len(crn2)
    r = []
    for i in range(len(crn2)):
        rr = []
        for j in range(len(crn1)):
            (t1, t2) = sdiff(crn1[j][0], crn2[i][0])
            (t3, t4) = sdiff(crn1[j][1], crn2[i][1])
            if set(t2).isdisjoint(set(fs)) and set(t4).isdisjoint(set(fs)):
                if t1 == []:
                    if t3 == []:
                        rr.append(True)
                    else:
                        if t4 == []:
                            rr.append(False)
                        else:
                            rr.append(True)
                else:
                    if t2 == []:
                        rr.append(False)
                    else:
                        if t3 == []:
                            rr.append(True)
                        else:
                            if t4 == []:
                                rr.append(False)
                            else:
                                rr.append(True)
            else:
                rr.append(False)
        (t1, t2) = sdiff(crn2[i][0], crn2[i][1])
        if (set(t1).isdisjoint(set(fs)) or not set(t2).issubset(set(fs))) and (set(t2).isdisjoint(set(fs)) or not set(t1).issubset(set(fs))):
            rr.append(True)
        else:
            rr.append(False)
        r.append(rr)
    return r

def perm(fs, inter):
# check the permissive condition, using the global original formal crn and original implementation crn (without substitutions).
    global intr, fcrn, icrn, f, permissive_depth, max_depth, permissive_failure
    ist = []
    tr = []
    fr = []
    hasht = set([])
    f = True
    crn2 = subst(icrn, inter[0], inter[1])
    T = update(fcrn, crn2, fs)   

    def cnstr(s, c):   
        # for n-tuple of formal reactants s, construct all implementation n-tuples of implementation species
        # such that impl-n-tuple[i] is interpreted to a superset of formal-n-tuple[i].
        # These implementation n-tuples are accumulated in "ist" (adding to its previous content).
        # Note that the variables "ist" and "c" must be cleared before calling cnstr().
        if s == []:
            ist.append(list(c))
        else:
            ss = list(s)
            cc = list(c)
            t = ss.pop()
            cnstr(ss,cc)
            cc.append(0)
            l = len(cc)-1
            for i in t:
                cc[l] = i
                cnstr(ss,cc)
        return []
    
    def search(s, d):
        # s is the implementation state, d is the depth.
        # try to find (via trivial reactions) a state in which a reaction in fr can fire.
        # fr is a particular formal reaction along with all implementation reactions that interpret to it. 
        # tr is the list of all trivial reactions in the implementation.
        # fail if depth d of trivial reaction steps is exceeded without firing a reaction in fr.
        global f
        if d > permissive_depth:
            f = False
            return False
        s.sort()
        tmp = ''
        for i in s:
            tmp += i
        if tmp in hasht:
            return False
        else:
            hasht.add(tmp) 
        for i in fr:
            if diff(i[0], s) == []:
                return True
        for i in tr:
            if diff(i[0], s) == []:
                t = diff(s, i[0])
                t.extend(i[1])
                if search(t, d+1):
                    return True
        return False
    
    n = 0
    # build tr just once
    for i in T:
        if i.pop():
            tr.append(icrn[n])
        n += 1
    for i in range(len(fcrn)):
        # build fr for this formal reaction
        fr = []
        for j in range(len(icrn)):
            if T[j][i]:
                fr.append(icrn[j])
        fr.append(fcrn[i])   

        # build a list of possible implementation species for each formal species.
        # i.e. s[k] is a list containing formal species k and each implementation species that interprets to a superset containing it.
        s = []
        n = 0
        for j in fcrn[i][0]:
            s.append([j])
            for k in range(len(inter[0])):
                if j in inter[1][k]:
                    s[n].append(inter[0][k])
            n += 1
        ist = []
        c = []
        cnstr(s, c)
        # At this point, "ist" contains an exhaustive and sufficient list of possible initial implementation states 
        # in which the current formal reaction #i must be able to fire (via a series of trivial reactions),
        # specifically states with one implementation species per formal species, along with their subsets 
        # (in case the a species interprets to multiple reactants).
        # Note that we will only want to test states in "ist" that interpret to a state in which #i can fire.
        
        tested = []  # to avoid testing a superset of some state that's already been tested
        for j in ist:
            if diff(fcrn[i][0], j) == []:
                continue
            tmp = subst([[j,[]]], inter[0], inter[1])
            if diff(fcrn[i][0], tmp[0][0]) == []:  # only test if reactants j interpret to allow #i to fire
                t = False
                for k in tested:
                    if diff(k, j) == []:
                        t = True
                        break
                if t:
                    continue
                hasht = set([])
                if not(search(j,0)): # didn't work, but keep track of stuff for user output
                    if max_depth == -2:
                        return False
                    if f == True:   # global variable saying whether permissive has been tried or not
                        max_depth = -1
                    else:
                        max_depth = -2
                    intr = list(inter)
                    permissive_failure[0] = fcrn[i]
                    permissive_failure[1] = j
                    return False
                else:
                    tested.append(j)
    intr = list(inter)
    return True

def equations(crn1, crn2, fs, unknow, inter):
    # All unknown implementation reactions (i.e. those with some unknown species) must be trivial.
    # Build the matrix for the the "solve" function, to see whether the interpretation can be completed as required.
    # Note that crn2 has been substituted already according to the (partial) interpretation "inter".
    # Also note that "unknow" is not used; the unknown species are recalculated here
    #  because "unknow" contains only those implementation reactions that have not been solved by the row search,
    #  but other rows (i.e. implementation reactions) may be implicitly solved by the partial interpretation.
    unknown = []
    ust = set([])
    for i in range(len(crn2)):
        tmp = set(crn2[i][0])-set(fs)
        tmp |= set(crn2[i][1])-set(fs)
        ust |= tmp
        if tmp != set([]):
            unknown.append(i)  # "unknown" is the list of implementation reactions with an unknown species
    us = list(ust)  # all species that remain unknown in current (partial) interpretation
    n = 0
    a = []
    for i in unknown:
        a.append([])
        for j in us:
            a[n].append(crn2[i][0].count(j)-crn2[i][1].count(j))
        a[n].append(0)
        n += 1
    itmp = copy.deepcopy(inter)
    l = len(inter[0])
    for i in us:
        itmp[0].append(i)
        itmp[1].append([]) 
    for i in fs:
        n = 0
        for j in unknown:
            a[n].pop()
            a[n].append(crn2[j][0].count(i)-crn2[j][1].count(i))
            n += 1
        s = solve(a)
        if s == []:
            return False
        else:
            for j in range(len(us)):
                itmp[1][j+l].extend([i for k in range(s[j])])
    return perm(fs, itmp)

def searchr(crn1, crn2, fs, unknown, inter, d):
# Row search.  I.e. make sure every implementation reaction can interpret to a formal reaction (or be trivial).
    global max_depth, intr
    if unknown == []:
        return perm(fs, inter)
    T = update(crn1, crn2, fs)
    if not checkT(T):
        return False
    if max_depth >= 0 and d > max_depth:
        intr = list(inter)
        max_depth = d
    min = len(crn1)+1
    k = -1  # next row reaction we will solve
    nt = 0  # number of possibly trivial reactions according to table
    for i in unknown:
        tmp = T[i].count(True)
        if T[i][len(crn1)] == True:
            nt += 1
            tmp -= 1
        if tmp < min and tmp > 0:
            min = tmp
            k = i
    if nt == len(unknown):  # all unsearched rows could be trivial -- try it!
        if equations(crn1, crn2, fs, unknown, inter):
            return True
    if k < 0:
        return False
    untmp = list(unknown)
    untmp.remove(k)
    if T[k][len(crn1)] == True:  # if implementation reaction #k can be trivial, leave it that way and try more rows
        if searchr(crn1, crn2, fs, untmp, inter, d):
            return True
    n = 0
    for c in range(len(crn1)): # try to match implementation reaction #k with some formal reaction
        if T[k][c]:        
            ul = diff(crn2[k][0], crn1[c][0])
            nl = len(ul)
            sl = diff(crn1[c][0], crn2[k][0])
            tmpl = enum(nl, sl)
            ur = diff(crn2[k][1], crn1[c][1])
            nr = len(ur)
            sr = diff(crn1[c][1], crn2[k][1])
            tmpr = enum(nr, sr)
            ul.extend(ur)
            nl += nr
            for i in tmpl:
                for j in tmpr:
                    tmpi = list(i)
                    tmpi.extend(j)
                    crntmp = subst(crn2, ul, tmpi)
                    itmp = [[],[]]
                    itmp[0] = list(inter[0])
                    itmp[0].extend(ul)
                    itmp[1] = list(inter[1])
                    itmp[1].extend(tmpi)                       
                    if searchr(crn1, crntmp, fs, untmp, itmp, d+1):
                        return True
    return False

def searchc(crn1, crn2, fs, unknown, inter, d):
    # Search column.  I.e. make sure every formal reaction can be implemented.
    global max_depth, intr
    T = update(crn1, crn2, fs)
    if not checkT(T):
        return False
    if max_depth >= 0 and d > max_depth:
        intr = list(inter)
        max_depth = d
    min = len(crn2)+1
    c = -1  # this will be the next column to solve, if possible
    for i in unknown:
        tmp = 0
        for j in range(len(crn2)):
            if T[j][i]:
                tmp += 1
        if tmp < min:
            min = tmp
            c = i
    if c < 0:  # done with column search.  transition to row search!
        untmp = []
        for i in range(len(crn2)):
            if not (set(crn2[i][0])-set(fs) == set([]) and set(crn2[i][1])-set(fs) == set([])):
                untmp.append(i)
        return searchr(crn1, crn2, fs, untmp, inter, d)
    else:
        untmp = list(unknown)
        untmp.remove(c)
        n = 0
        for k in range(len(crn2)):
            if T[k][c]:
                ul = diff(crn2[k][0], crn1[c][0])
                nl = len(ul)
                sl = diff(crn1[c][0], crn2[k][0])
                tmpl = enum(nl, sl)
                ur = diff(crn2[k][1], crn1[c][1])
                nr = len(ur)
                sr = diff(crn1[c][1], crn2[k][1])
                tmpr = enum(nr, sr)
                ul.extend(ur)
                nl += nr
                for i in tmpl:
                    for j in tmpr:
                        tmpi = list(i)
                        tmpi.extend(j)
                        crntmp = subst(crn2, ul, tmpi)
                        itmp = [[],[]]
                        itmp[0] = list(inter[0])
                        itmp[0].extend(ul)
                        itmp[1] = list(inter[1])
                        itmp[1].extend(tmpi)                                                
                        if searchc(crn1, crntmp, fs, untmp, itmp, d+1):
                            return True
    return False

def test(c1, c2, verbose = True, inter = [[],[]]):
    (crn1, fs1) = c1[:2]
    (crn2, fs2) = c2[:2]
    global fcrn, icrn, intr, max_depth, permissive_failure
    fcrn = crn1
    icrn = crn2
    print "Original CRN:"
    for rxn in crn1:
        print "   ",
        printRxn(rxn)
    print
    if crn2 == []:
        print "Compiled CRN is empty"
        print
        return crn1 == crn2
    print "Compiled CRN:"
    for rxn in crn2:
        print "   ",
        printRxn(rxn)
    print
    unknown = [i for i in range(len(crn1))]
    if searchc(crn1, crn2, fs2, unknown, inter, 0):
        print "Valid interpretation :"
        output(intr)
        return True
    else:
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
    crn1 = [[['a'],['b']]]
    fs1 = []
    cs1 = []
    crn2 = [[['a1'],['b1']],[['x'],['a1']],[['x'],['b1']],[['y'],['b1']],[['y'],['a1']],[['x'],['a0']],[['a0'],['a1']]]
    fs2 = ['a','b']
    cs2 = []
    v = test((crn1, fs1, cs1), (crn2, fs2, cs2))
    if v:
        print "verify: compilation was correct."
    else:
        print "verify: compilation was incorrect."
