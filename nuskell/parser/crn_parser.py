#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Parser module for chemical reaction network description files (*.crn).
#

from pyparsing import *

def crn_document_setup():
  """ Parse a chemical reaction network
  CRN:
    A + B -> C + D
    -> A
    X ->
    A + C <=> X
  returns: 
  (
    [['irreversible', ['A', 'B'], ['C','D']], 
     ['irreversible', [], ['A']], 
     ['irreversible', ['X'], []], 
     ['reversible', ['A', 'C'], ['X']], 
    ['A', 'B', 'C', 'D', 'X'], 
    []
  )

  -- this can read other input as well, but I haven't figured that out yet.
  """
  
  W = Word
  G = Group
  S = Suppress
  O = Optional
  L = Literal
  
  def T(x, tag):
    """ Return a *Tag* to distinguish (ir)reversible reactions """
    def TPA(tag):
      return lambda s, l, t: [tag] + t.asList()
    return x.setParseAction(TPA(tag))
  
  crn_DWC = "".join(
      [x for x in ParseElementEnhance.DEFAULT_WHITE_CHARS if x != "\n"])
  ParseElementEnhance.setDefaultWhitespaceChars(crn_DWC)
  
  identifier = W(alphas, alphanums+"_")
  
  reaction = T(G(O(delimitedList(identifier, "+"))) + \
             S("->") + \
             G(O(delimitedList(identifier, "+"))), "irreversible") + \
             S(LineEnd())
  rev_reaction = T(G(O(delimitedList(identifier, "+"))) + \
                 S("<=>") + \
                 G(O(delimitedList(identifier, "+"))), "reversible") + \
                 S(LineEnd())
  
  formal = G(L("formal") + S(L("=") + \
             L("{")) + delimitedList(identifier) + S("}"))
  constant = G(L("constant") + S(L("=") + \
               L("{")) + delimitedList(identifier) + S("}"))
  
  crn = ZeroOrMore(G(reaction | rev_reaction) | S(LineEnd())) + \
        O(formal) + ZeroOrMore(S(LineEnd())) + \
        O(constant) + ZeroOrMore(S(LineEnd()))
  
  document = StringStart() + crn + StringEnd()
  document.ignore(pythonStyleComment)
  return document

crn_document = crn_document_setup()
def species(crn):
  """ Returns the list of all the distinct species in the given CRN """
  species = set()
  for rxn in crn:
    # skip 'formal' or 'constant' statement
    if rxn[0] != "irreversible" and rxn[0] != "reversible": continue
    species = species.union(rxn[1]).union(rxn[2])
  return list(species)

def _post_process(crn):
  formal_species = species(crn)
  constant_species = []
  for x in crn:
    if x[0] == "formal":
      formal_species = x[1:]
    elif x[0] == "constant":
      constant_species = x[1:]
      formal_species = formal_species - constant_species
  for i in range(len(crn)):
    if crn[i][0] != "reversible" and crn[i][0] != "irreversible":
      crn = crn[:i]
      break
  return (crn, formal_species, constant_species)

def split_reversible_reactions(crn):
    """ called by basis_finder """
    new_crn = []
    for rxn in crn:
        if rxn[0] == "irreversible":
            new_crn.append(rxn[1:])
        if rxn[0] == "reversible":
            r = rxn[1]
            p = rxn[2]
            new_crn.append([r, p])
            new_crn.append([p, r])
    return new_crn

def parse_file(data):
    """Parses the given .crn file and returns the result in the form of
       (crn, formal species, constant species)."""
    crn = crn_document.parseFile(data, parseAll = True).asList()
    return _post_process(crn)

def parse_string(data):
    """Parses the given string and returns the result in the form of
       (crn, formal species, constant species)."""
    return _post_process(crn_document.parseString(data).asList())
