#!/usr/bin/env python
#
#
# Copyright (c) 2009-2015 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Finds the formal basis of the given CRN.
#

import sys, crn_parser

def print_reaction(rxn):
    first = True
    for x in rxn[0]:
        if first: first = False
        else: print "+",
        print x,
    print "->",
    first = True
    for x in rxn[1]:
        if first: first = False
        else: print "+",
        print x,
    print

def remove_const(crn, const):
    for rxn in crn:
        for x in const:
            while x in rxn[0]:
                rxn[0].remove(x)
            while x in rxn[1]:
                rxn[1].remove(x)
    return crn

def formal(s, fs):
    return filter(lambda x: x in fs, s)

def intermediate(s, fs):
    return filter(lambda x: x not in fs, s)

# computes the next state, assuming rxn can occur in s
def next_state(s, rxn):
    s = s[:]
    for x in rxn[0]:
        s.remove(x)
    s += rxn[1]
    return sorted(s)

# minimal_initial_state
def minimal_initial_state(pathway):
    initial = []
    current = []
    for rxn in pathway:
        for r in rxn[0]:
            if r in current:
                current.remove(r) 
            else:
                initial.append(r)
        for r in rxn[1]:
            current.append(r)
    return sorted(initial)

def contained(a, b):
    def contained_helper(a, b):
        if len(a) == 0:
            return True
        if len(b) < len(a):
            return False
        if a[0] == b[0]:
            return contained_helper(a[1:], b[1:])
        else:
            return contained_helper(a, b[1:])
    a = sorted(a)
    b = sorted(b)
    return contained_helper(a, b)

def final_state(p, S):
    c = S[:]
    for rxn in p:
        for x in rxn[0]:
            if x not in c: return None
            c.remove(x)
        for x in rxn[1]:
            c.append(x)
    return sorted(c)

def setminus(a, b):
    a = a[:]
    b = b[:]
    T = []
    for x in a:
        if x in b:
            b.remove(x)
        else:
            T.append(x)
    return T

def remove_duplicates(l):
    r = []
    if len(l) == 0: return []
    l.sort()
    while len(l) > 1:
        if l[0] != l[1]:
            r.append(l[0])
        l = l[1:]
    r.append(l[0])
    return r

def formal_state(S, fs):
    return len(formal(S, fs)) == len(S)

def decompose(p, fs):
    def helper(p1, p2, remaining, fs):
        i1 = minimal_initial_state(p1)
        i2 = minimal_initial_state(p2)
        if not formal_state(i1, fs) or not formal_state(i2, fs):
            return []
        if len(remaining) == 0:
            if len(p1) > 0 and len(p2) > 0:
                # valid decomposition
                return [[final_state(p1, i1), final_state(p2, i2)]]
            else:
                return []
        else:
            return helper(p1 + [remaining[0]], p2, remaining[1:], fs) \
                 + helper(p1, p2 + [remaining[0]], remaining[1:], fs)
    l = remove_duplicates(sorted(helper([], [], p, fs)))
    return l

def formal_closure(p, fs):
    initial = minimal_initial_state(p)
    S = initial[:]
    T = formal(S, fs)
    for r in p:
        S = next_state(S, r)
        f = formal(S, fs)
        T = T + setminus(f, T)
    return sorted(T)

def regular_final_state(p, fs):
    initial = minimal_initial_state(p)
    S = initial[:]
    state_list = [S]
    flag = True
    T = []
    RFS = []
    j = 0
    for r in p:
        S = next_state(S, r)
        state_list.append(S)
        f = formal(S, fs)
        if not contained(f, initial):
            flag = False
        if flag: j += 1
    i = 0
    for r in reversed(p):
        i += 1
        (R,P) = r
        f = formal(state_list[-i], fs)
        T = T + setminus(f, T)
        if j >= len(p) - i and \
           len(formal(setminus(state_list[-(i+1)],R),fs)) == 0:
            T.sort()
            if T not in RFS:
                RFS.append(T)
    return sorted(RFS)

def width(p):
    initial = minimal_initial_state(p)
    S = initial[:]
    w = len(S)
    for r in p:
        S = next_state(S, r)
        if len(S) > w: w = len(S)
    return w

def enumerate(p, w_max, i_max, crn, fs):
    global ret, ccheck, ebasis, done, inter
    if done: return
    w_p = width(p)
    if w_p > w_max: return
    initial = minimal_initial_state(p)
    if not formal_state(initial, fs): return
    if len(initial) > i_max: return
    final = final_state(p, initial)
    DFS = decompose(p, fs)
    # strong decomposability
    for [p1, p2] in DFS:
        if formal_state(p1, fs) or formal_state(p2, fs):
            return
    fc = formal_closure(p, fs)
    RFS = regular_final_state(p, fs)
    sig = (initial, final, w_p, fc, DFS, RFS)

    if sig in ret: return
    ret.append(sig)

    #print p
    #print "   ",(initial, final, w_p, fc, len(DFS), len(RFS))
    #if len(initial)>2 and formal_state(final, fs) and len(DFS) == 0: print p
    #print "   ",final, len(DFS)
    #print "     ",DFS
    if len(DFS) == 0:
        if final not in ccheck:
            ccheck.append(final)
            if not tidy(final, crn, fs):
                print "The given system is not tidy:"
                print "from initial state ", initial
                for rxn in p:
                    print_reaction(rxn)
                print
                done = 1
                return
        if len(p) > 0 and formal_state(final, fs):
            ### following is test for relaxed regularity condition
            def collapse(l):
                l2 = []
                for x in l:
                    if x in inter.keys():
                        y = inter[x]
                    else:
                        y = [x]
                    l2 += y
                return l2
            p1 = map(lambda rxn: [collapse(rxn[0]), collapse(rxn[1])], p)
            initial1 = minimal_initial_state(p1)
            final1 = final_state(p1, initial1)
            RFS1 = regular_final_state(p1, fs)
            ebasis.append(p)
            if final1 not in RFS1:
            ### test code ends (if delete, revive the following lines)
            #ebasis.append(p)
            #if final not in RFS:
                print "The given system is not regular:"
                print "from initial state ", initial
                for rxn in p:
                    print_reaction(rxn)
                print
                done = 1
                return
    for r in crn:
        enumerate(p + [r], w_max, i_max, crn, fs)
    
def tidy(S, crn, fs):
    # BFS. May have to eventually rewrite to impose a width bound.
    queue = [intermediate(S, fs)]
    mem = [queue[0]]
    while len(queue) > 0:
        S = queue[0]
        queue = queue[1:]
        if len(S) == 0: return True
        for rxn in crn:
            if len(rxn[0]) > 0 and contained(rxn[0], S):
                nS = intermediate(next_state(S, rxn), fs)
                if nS not in mem:
                    mem.append(nS)
                    queue.append(nS)
    return False

def find_basis(crn, fs, inter2 = {}):
    global ret, ccheck, ebasis, done, inter
    done = 0
    inter = inter2

    # fuel preprocessing: anything that can be produced by a reaction that
    # consumes only fuel species is also considered a fuel species.
    if False:
        fuel = []
        while True:
            flag = False
            for [R, P] in crn:
                spontaneous = True
                for x in R:
                    if x not in fuel:
                        spontaneous = False
                        break
                if spontaneous:
                    for x in P:
                        if x not in fs and x not in fuel:
                            fuel.append(x)
                            flag = True
            if not flag:
                break
        remove_target = []
        for [R, P] in crn:
            for x in fuel:
                while x in R:
                    R.remove(x)
                while x in P:
                    P.remove(x)
            R.sort()
            P.sort()
            if R == P:
                remove_target.append([R, P])
        for r in remove_target:
            crn.remove(r)

    # bf and br are for optimizations explained in 5.6 of SWS's Master's
    # thesis.
    b = 0 # the branching factor
    bf = 0 # max over (R,P) of |intermediate(R)|
    br = [] # set of (|formal(R)|,|intermediate(R)|)
    for [r, p] in crn:
        if len(r) > b: b = len(r)
        if len(intermediate(r, fs)) > bf: bf = len(intermediate(r, fs))
        if len(p) > b: b = len(p)
        br.append([len(formal(r, fs)), len(intermediate(r, fs))])
    br.sort()
    br = remove_duplicates(br)

    w_max = 0
    i_max = 0
    ccheck = []
    while not done:
        print "Current bounds : w_max", w_max, "i_max", i_max
        ret = []
        ebasis = []
        enumerate([], w_max, i_max, crn, fs)
        signatures = ret
        current_w = 0
        current_i = 0
        for (i, f, w, fc, dfs, rfs) in signatures:
            if len(dfs) == 0:
                if w > current_w:
                    current_w = w
                if len(i) > current_i:
                    current_i = len(i)
        w_t = current_w * bf + b 
        i_t = 0
        for [x, y] in br:
            i_t = max(i_t, current_i * y + x)
        if w_t <= w_max and i_t <= i_max: break
        w_max = w_t
        i_max = i_t
    
    if done: return None

    fbasis = []
    fbasis_raw = []
    for p in ebasis:
        ### following is test for history domain/waste handling
        def collapse(l):
            l2 = []
            for x in l:
                if x in inter.keys():
                    y = inter[x]
                else:
                    y = [x]
                l2 += y
            return l2
        initial = minimal_initial_state(p)
        final = final_state(p, initial)
        fbasis_raw.append([sorted(initial), sorted(collapse(final))])
        p1 = map(lambda rxn: [collapse(rxn[0]), collapse(rxn[1])], p)
        initial = minimal_initial_state(p1)
        final = final_state(p1, initial)
        ### test code ends (if delete, revive the following lines)
        #initial = minimal_initial_state(p)
        #final = final_state(p, initial)
        fbasis.append([sorted(initial), sorted(final)])
    fbasis = remove_duplicates(fbasis)
    # permissive test
    interrev = {}
    for x in inter.keys():
        for y in inter[x]:
            if y not in interrev.keys():
                interrev[y] = [[x]]
            else:
                interrev[y].append([x])
    for rxn in fbasis:
        def cartesian_product(l):
            if len(l) == 0:
                return []
            if len(l) == 1:
                return l[0]
            r = []
            for i in l[0]:
                for j in l[1]:
                    r.append(i+j)
            return cartesian_product([r]+l[2:])
        initial_states = cartesian_product(map(lambda x: interrev[x], rxn[0]))
        for initial in initial_states:
            initial = sorted(initial)
            flag = False
            for r in fbasis_raw:
                if r[0] == initial and r[1] == rxn[1]:
                    flag = True
                    break
            if not flag:
                print "Invalid history domain:"
                print "  Cannot get from ",initial," to",rxn[1]
                return None
    # permissive test end
    return fbasis

if __name__ == "__main__":
    crn_file = sys.argv[1]

    # Add appropriate extensions if necessary.
    if len(crn_file) < 4 or crn_file[-4:] != ".crn": crn_file += ".crn"

    (crn, formal_species, const_species) = crn_parser.parse_file(crn_file)
    crn = crn_parser.split_reversible_reactions(crn)

    crn = remove_const(crn, const_species)
    basis = find_basis(crn, formal_species)
    print "Formal basis :"
    for (r, p) in basis:
        print_reaction([r, p])
