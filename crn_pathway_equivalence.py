import time, basis_finder, string, crn_bisimulation_equivalence, sys

def printRxn(rxn, inter = {}):
    first = True
    for x in rxn[0]:
        if x[0] not in string.letters:
            if x in inter.keys() and inter[x]==[]:
                x = "w" + x
            else:
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
            if x in inter.keys() and inter[x]==[]:
                x = "w" + x
            else:
                x = "i" + x
        if not first:
            print "+",
        else:
            first = False
        print x,
    print

def findWastes(crn, formal):
    species = set()
    for rxn in crn:
        for x in rxn[0]:
            species.add(x)
        for x in rxn[1]:
            species.add(x)
    nonwastes = set(formal)
    while True:
        flag = False
        for x in species:
            if x not in nonwastes:
                for rxn in crn:
                    if x in rxn[0] and \
                      (len(nonwastes & set(rxn[1])) > 0 \
                      or len(nonwastes & set(rxn[0])) > 0):
                        flag = True
                        nonwastes.add(x)
                        break
        if not flag: break

    return species - nonwastes

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

def test(c1, c2, inter, verbose = True, integrated = False, interactive = False):
    (crn1, fs1) = c1
    (crn2, fs2) = c2
    #for rxn in crn2:
    #    print rxn
    #print "------"
    crn1 = [[sorted(rxn[0]), sorted(rxn[1])] for rxn in crn1]
    crn2 = [[sorted(rxn[0]), sorted(rxn[1])] for rxn in crn2]
    crn2.sort()

    wastes = findWastes(crn2, fs2)
    fs2 = (set(fs2)).union(wastes)
    for x in wastes: inter[x] = []

    # remove trivial reactions
    remove_target = []
    for [R, P] in crn2:
        R.sort()
        P.sort()
        if R == P:
            remove_target.append([R, P])
    for r in remove_target:
        crn2.remove(r)

    crn1size = len(crn1)
    crn2size = len(crn2)
    t1 = time.time()
    print "CRN 1 size :", crn1size
    print "CRN 2 size :", crn2size
    species = set()
    for rxn in crn2: species = species.union(set(rxn[0])).union(set(rxn[1]))
    print "# of species :", len(species)
    
    print "Original CRN:"
    for rxn in crn1:
        print "   ",
        printRxn(rxn)
    print
    print "Compiled CRN:"
    print "formal species = ", fs2
    for rxn in crn2:
        print "   ",
        printRxn(rxn, inter)
    print

    if interactive:
        print "Enter species to be treated as formal species, along with its interpretation:"
        print "(e.g. i187 -> A + B)"
        print "When done, press ctrl + D."
        for line in sys.stdin:
            z = map(lambda x: x.strip(), line.split("->"))
            y1 = z[0]
            y2 = map(lambda x: x.strip(), z[1].split("+"))
            if y1[0] == "i": y1 = y1[1:]
            inter[y1] = y2
            fs2.add(y1)
        print

    basis = basis_finder.find_basis(crn2, fs2, True, inter if integrated else None)
    if basis == None: # irregular or nontidy
        return False
    if integrated: # integrated hybrid
        (fbasis_raw, fbasis) = basis 

        for i in range(len(fbasis_raw)):
            fbasis_raw[i][0].sort()
            fbasis_raw[i][1].sort()
        fbasis_raw = remove_duplicates(fbasis_raw)
        for i in range(len(fbasis)):
            fbasis[i][0].sort()
            fbasis[i][1].sort()
        fbasis = remove_duplicates(fbasis)

        # bisimulation test
        fbasis2 = []
        for [initial, final] in fbasis_raw:
            def collapse(l):
                l2 = []
                for x in l:
                    if x in inter.keys():
                        y = inter[x]
                    else:
                        y = [x]
                    l2 += y
                return l2
            r = [sorted(initial), sorted(collapse(final))]
            fbasis2.append(r)
    else: # compositional hybrid
        fbasis_raw = basis

        for i in range(len(fbasis_raw)):
            fbasis_raw[i][0].sort()
            fbasis_raw[i][1].sort()
        fbasis_raw = remove_duplicates(fbasis_raw)

        # bisimulation test
        fbasis2 = []
        fbasis = []
        for [initial, final] in fbasis_raw:
            def collapse(l):
                l2 = []
                for x in l:
                    if x in inter.keys():
                        y = inter[x]
                    else:
                        y = [x]
                    l2 += y
                return l2
            r = [sorted(initial), sorted(collapse(final))]
            fbasis2.append(r)
            r = [sorted(collapse(initial)), sorted(collapse(final))]
            fbasis.append(r)
        fbasis = remove_duplicates(fbasis)
    # TODO : the following is not strictly correct because it tests for
    #        strong bisimulation instead of weak bisimulation.
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
            for r in fbasis2:
                if r[0] == initial and r[1] == rxn[1]:
                    flag = True
                    break
            if not flag:
                print "Permissive test failed:"
                print "  Cannot get from ",initial," to",rxn[1]
                return None
    # permissive test end
    basis = fbasis

    print "Basis of the compiled CRN:"
    for rxn in basis:
        print "   ",
        printRxn(rxn)
    print

    # delimiting test
    flag = True
    for rxn in crn1:
        if rxn not in basis:
            reactants = {}
            for x in rxn[0]:
                if x in reactants.keys():
                    reactants[x] += 1
                else:
                    reactants[x] = 1
            products = {}
            for x in rxn[1]:
                if x in products.keys():
                    products[x] += 1
                else:
                    products[x] = 1
            if reactants != products:
                print "Error : The formal pathway"
                print "    ",
                printRxn(rxn)
                print " is in the input CRN but not in the compiled CRN."
                flag = False
    for rxn in basis:
        if rxn not in crn1:
            reactants = {}
            for x in rxn[0]:
                if x in reactants.keys():
                    reactants[x] += 1
                else:
                    reactants[x] = 1
            products = {}
            for x in rxn[1]:
                if x in products.keys():
                    products[x] += 1
                else:
                    products[x] = 1
            if reactants != products:
                print "Error : The formal pathway"
                print "    ",
                printRxn(rxn)
                print " is in the compiled CRN but not in the input CRN."
                flag = False


    t2 = time.time()
    print "Elapsed time :", t2-t1
    return flag
