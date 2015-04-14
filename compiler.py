#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Compiler module.
#

import ts_parser, crn_parser, interpreter

def compile(ts_file, crn_file):
    """Compiles the given .crn file using the given .ts file."""
    # read the input files
    ts = ts_parser.parse_file(ts_file)
    (crn, formal_species, const_species) = crn_parser.parse_file(crn_file)

    domains, strands, formal_species, constant_species = \
        interpreter.interpret(ts, crn, formal_species)

    outfile = crn_file[:-4] + ".pil"
    print_as_PIL(outfile, domains, strands, \
                 formal_species, constant_species)

    outfile = crn_file[:-4] + ".dom"
    print_as_DOM(outfile, domains, strands, \
                 formal_species, constant_species)

def print_as_PIL(outfile, domains, strands, \
                 formal_species, constant_species):
    F = open(outfile, 'w')
    for x in domains:
        if x.name != "?":
            print >> F, "sequence " + x.name + " = " + "N" * x.length + \
                        " : " + str(x.length)
    for x in strands:
        total_length = 0
        print >> F, "strand " + x.name + " =",
        for d in x.domain_list:
            total_length += d.length
            print >> F, d.name,
        print >> F, ": " + str(total_length)

    print >> F, "# Formal species"
    for x in formal_species:
        print >> F, "structure " + x.name + " =",
        first = True
        for s in x.strand_list:
            if first: first = False
            else: print >> F, "+",
            print >> F, s.name,
        print >> F, ": "  + x.structure

    print >> F, "# Constant species"
    for x in constant_species:
        print >> F, "structure " + x.name + " =",
        first = True
        for s in x.strand_list:
            if first: first = False
            else: print >> F, "+",
            print >> F, s.name,
        print >> F, ": "  + x.structure

    # TODO : figure out a way to deal with the wildcards.

def print_as_DOM(outfile, domains, strands, \
                 formal_species, constant_species):
    F = open(outfile, 'w')
    for x in domains:
        if x.name != "?":
            print >> F, "sequence " + x.name + " : " + str(x.length)

    print >> F, "# Formal species"
    for x in formal_species:
        print >> F, x.name + " :"
        first = True
        t = []
        for s in x.strand_list:
            if first: first = False
            else: print >> F, "+",
            for d in s.domain_list:
                print >> F, d.name,
                t.append(d)
        print >> F
        i = 0
        while i < len(x.structure):
            if x.structure[i] == "+":
                print >> F, "+",
                i += 1
            else:
                if len(t) == 0:
                    break
                elif t[0].length == 0:
                    t = t[1:]
                else:
                    print >> F, x.structure[i],
                    i += t[0].length
                    t = t[1:]
        print >> F
                
    print >> F, "# Constant species"
    for x in constant_species:
        print >> F, x.name + " :"
        first = True
        t = []
        for s in x.strand_list:
            if first: first = False
            else: print >> F, "+",
            for d in s.domain_list:
                print >> F, d.name,
                t.append(d)
        print >> F
        i = 0
        while i < len(x.structure):
            if x.structure[i] == "+":
                print >> F, "+",
                i += 1
            else:
                if len(t) == 0:
                    break
                elif t[0].length == 0:
                    t = t[1:]
                else:
                    print >> F, x.structure[i],
                    i += t[0].length
                    t = t[1:]
        print >> F
