import time, basis_finder, string, crn_bisimulation_equivalence

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

def test(c1, c2, verbose = True, inter = [[],[]]):
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

    basis = basis_finder.find_basis(crn2, fs2, inter)
    if basis == None: # irregular or nontidy
        return False

    for i in range(len(basis)):
        basis[i][0].sort()
        basis[i][1].sort()

    print "Basis of the compiled CRN:"
    for rxn in basis:
        print "   ",
        printRxn(rxn)
    print

    #print "Proposed interpretation:"
    #for x in inter.keys():
    #    print x,"=>",inter[x]

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
