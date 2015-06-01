#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Parser module for domain specification description files (*.dom).
#

from pyparsing import *

ParserElement.setDefaultWhitespaceChars(" \t")

W = Word
G = Group
S = Suppress
O = Optional
L = Literal

identifier = W(alphas, alphanums)
number = W(nums, nums)
domain = G(OneOrMore((G(identifier + O("*")) | "?" | "+")))
dotparen = G(OneOrMore(W("().+", max=1)))
structure = domain + OneOrMore(LineEnd().suppress()) + dotparen
sequence = G(S("sequence") + identifier + S(":") + number + OneOrMore(LineEnd().suppress()))
molecule = G(identifier + S(":") + OneOrMore(LineEnd().suppress()) + structure)

document = StringStart() + ZeroOrMore(LineEnd().suppress()) + G(ZeroOrMore(sequence)) + G(OneOrMore(molecule + OneOrMore(LineEnd().suppress()))) + StringEnd()
document.ignore(pythonStyleComment)

def parse(data):
    return document.parseFile(data).asList()
