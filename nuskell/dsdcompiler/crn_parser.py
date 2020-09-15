#
#  nuskell/dsdcompiler/crn_parser.py
#  NuskellCompilerProject
#
# Copyright (c) 2010-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
from collections import namedtuple
from pyparsing import (Word, Literal, Group, Suppress, Combine, Optional, ParseException,
                       alphas, nums, alphanums, delimitedList, StringStart, StringEnd, LineEnd,
                       ZeroOrMore, OneOrMore, pythonStyleComment, ParseElementEnhance)

Reaction = namedtuple('reaction', 'reactants products k_fwd k_rev')

class CRNParseError(Exception):
    pass

def crn_document_setup(modular=False):
    """Parse a formal chemical reaction network.

    Args:
      modular <optional:bool>: Adds an additional nesting for modules within a
        CRN. Use one line per module (';' separates reactions).

    Format:
      # A list of reactions, optionally with reaction rates:
      # <- this is a comment!
      B + B -> C    # [k = 1]
      C + A <=> D   # [kf = 1, kr = 1]
      <=> A  [kf = 15, kr = 6]

      # Note that you can write multiple reactions in one line:
      A + 2C -> E [k = 13.78]; E + F <=> 2A  [kf = 13, kr = 14]

    Returns:

    """
    # NOTE: If you want to add support for multiple modules per line, you can use
    # the '|' character.

    W = Word
    G = Group
    S = Suppress
    O = Optional
    C = Combine
    L = Literal

    def T(x, tag):
        """ Return a *Tag* to distinguish (ir)reversible reactions """
        def TPA(tag):
            return lambda s, l, t: [tag] + t.asList()
        return x.setParseAction(TPA(tag))

    crn_DWC = "".join(
        [x for x in ParseElementEnhance.DEFAULT_WHITE_CHARS if x != "\n"])
    ParseElementEnhance.setDefaultWhitespaceChars(crn_DWC)

    identifier = W(alphas, alphanums + "_")

    multiplier = W(nums)
    species = G(O(multiplier) + identifier)

    number = W(nums, nums)
    num_flt = C(number + O(L('.') + number))
    num_sci = C(number + O(L('.') + number) + L('e') + O(L('-') | L('+')) + W(nums))
    gorf = num_sci | num_flt

    # Make specification of forward, backward, reverse more flexible
    kf = S('kf') | S('fw')
    kr = S('kr') | S('bw') | S('rv')

    k = G(S('[') + O(S('k') + S('=')) + gorf + S(']'))
    rev_k = G(S('[') + kf + S('=') + gorf + S(',') + kr + S('=') + gorf + S(']')) | G(S('[') + gorf + S(',') + gorf + S(']'))

    concentration = T(species + S('@') + G(L("initial") | L("i") | L("constant") | L("c")) + G(gorf), 'concentration')

    reaction = T(G(O(delimitedList(species, "+"))) +
                 S("->") +
                 G(O(delimitedList(species, "+"))) + O(k), 'irreversible')

    rev_reaction = T(G(O(delimitedList(species, "+"))) +
                     S("<=>") +
                     G(O(delimitedList(species, "+"))) + O(rev_k), 'reversible')

    expr = G(reaction | rev_reaction | concentration)

    if modular:
        module = G(expr + ZeroOrMore(S(";") + expr))
    else:
        module = expr + ZeroOrMore(S(";") + expr)

    crn = OneOrMore(module + ZeroOrMore(S(LineEnd())))

    document = StringStart() + ZeroOrMore(S(LineEnd())) + crn + StringEnd()
    document.ignore(pythonStyleComment)
    return document

def post_process(crn, defaultrate = 1, defaultmode = None, defaultconc = None):
    """ Process a parsed CRN to return namedtuples. """
    def flint(inp):
        return int(inp) if float(inp) == int(float(inp)) else float(inp)

    def remove_multipliers(species):
        flat = []
        for s in species:
            if len(s) == 1:
                flat.append(s[0])
            elif len(s) == 2:
                ss = [s[1]] * int(s[0])
                flat.extend(ss)
        return flat

    new = []
    species = dict()
    for line in crn:
        if line[0] == 'concentration':
            spe = line[1][0]
            ini = 'initial' if line[2][0][0] == 'i' else 'constant'
            num = line[3][0]
            species[spe] = (ini, flint(num))
            continue
        elif len(line) == 3:
            # No rate specified
            t, r, p = line
            r = remove_multipliers(r)
            p = remove_multipliers(p)
            if t == 'reversible':
                new.append(Reaction(r, p, defaultrate, defaultrate))
            elif t == 'irreversible':
                new.append(Reaction(r, p, defaultrate, 0))
            else:
                raise CRNParseError('Wrong CRN format!')
        elif len(line) == 4:
            t, r, p, k = line
            r = remove_multipliers(r)
            p = remove_multipliers(p)
            if t == 'reversible':
                assert len(k) == 2
                new.append(Reaction(r, p, flint(k[0]), flint(k[1])))
            elif t == 'irreversible':
                assert len(k) == 1
                new.append(Reaction(r, p, flint(k[0]), 0))
            else:
                raise CRNParseError('Wrong CRN format!')
        else:
            raise CRNParseError('Wrong CRN format!')
        for s in r + p:
            if s not in species:
                species[s] = (defaultmode, defaultconc)
    return new, species

def parse_crn_file(filename, process = True, **kwargs):
    """Parses a CRN from a file.

    Args:
      filename (<str>): The name of the file for parsing.
      process (optional <bool>): Process the format of the returned CRN.

    Returns:
      crn (<lol>): List of list representation of a CRN
      species (<set()>): A set of all involved species (only when process=True)

    """
    crn_document = crn_document_setup()
    if process:
        return post_process(crn_document.parseFile(
            filename, parseAll=True).asList(), **kwargs)
    else:
        return crn_document.parseFile(filename, parseAll=True).asList()

def parse_crn_string(data, process = True, **kwargs):
    """Parses a CRN from a string.

    Args:
      filename (<str>): The name of the file for parsing.
      process (optional <bool>): Process the format of the returned CRN.

    Returns:
      crn (<lol>): List of list representation of a CRN
      species (<set()>): A set of all involved species (only when process=True)

    """
    crn_document = crn_document_setup()
    if process:
        return post_process(crn_document.parseString(data).asList(), **kwargs)
    else:
        return crn_document.parseString(data).asList()


