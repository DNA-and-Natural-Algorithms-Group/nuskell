#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Parser module for chemical reaction network description files (*.crn).
#

from pyparsing import (Word, Literal, Group, Suppress, Optional, ZeroOrMore,
    alphas, alphanums, delimitedList, StringStart, StringEnd, LineEnd, srange,
    pythonStyleComment, ParseElementEnhance)

def crn_document_setup():
  """Parse a chemical reaction network. 
  
  The crn document parser reads one chemical reaction per line, as well as
  special lines that define 'formal' and 'constant' species of the crn.

  Returns:
    a **pyparsing** document formating

  .. Input example:
  .. A + B -> C + D
  .. -> A
  .. X ->
  .. A + C <=> X
  .. formal = {A, B, C, D}
  .. constant = {}

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
  
  #identifier = W(alphas, alphanums+"_")
  identifier = W(srange("[a-zA-Z_]"), srange("[a-zA-Z0-9_]"))
  
  reaction = T(G(O(delimitedList(identifier, "+"))) + \
             S("->") + \
             G(O(delimitedList(identifier, "+"))), "irreversible")
  rev_reaction = T(G(O(delimitedList(identifier, "+"))) + \
                 S("<=>") + \
                 G(O(delimitedList(identifier, "+"))), "reversible")
  
  expr = G(reaction | rev_reaction) 

  reactline = expr + ZeroOrMore(S(";") + expr) + S(LineEnd())

  formal = G(L("formal") + S(L("=") + \
             L("{")) + delimitedList(identifier) + S("}"))
  constant = G(L("constant") + S(L("=") + \
               L("{")) + delimitedList(identifier) + S("}"))
  
  crn = ZeroOrMore(reactline | S(LineEnd())) + \
        O(formal) + ZeroOrMore(S(LineEnd())) + \
        O(constant) + ZeroOrMore(S(LineEnd()))
  
  document = StringStart() + crn + StringEnd()
  document.ignore(pythonStyleComment)
  return document

def species(crn):
  """Returns a list of all distinct species in the CRN.

  :param crn: a chemical reaction network 

  .. [['irreversible', ['A', 'B'], ['C','D']], 
  ..  ['irreversible', [], ['A']], 
  ..  ['irreversible', ['X'], []], 
  ..  ['reversible', ['A', 'C'], ['X']], 
  ..  ['A', 'B', 'C', 'D', 'X'], []]

  :return: list of species.
  """
  species = set()
  for rxn in crn:
    # skip 'formal' or 'constant' statement
    if rxn[0] != "irreversible" and rxn[0] != "reversible": continue
    species = species.union(rxn[1]).union(rxn[2])
  return list(species)

def _post_process(crn):
  """Take a crn and return it together with a list of formal and
  constant species.
  """

  all_species = species(crn)
  formal_species = []
  constant_species = []
  for x in crn:
    # overwrite the species assignment
    if x[0] == "formal":
      formal_species = x[1:]
    elif x[0] == "constant":
      constant_species = x[1:]

  #check that formal/constant assignments make sense!
  asp = set(all_species)
  fsp = set(formal_species)
  csp = set(constant_species)

  if csp :
    if fsp :
      if asp != fsp | csp :
        pass
        #raise ValueError, "missing species in CRN input"
      elif fsp & csp :
        raise ValueError, "species declared as formal & constant"
    else :
      fsp = asp - csp
  elif fsp :
    csp = asp - fsp
    if asp != fsp | csp :
      pass
      #raise ValueError, "missing species in CRN input"
    elif fsp & csp :
      raise ValueError, "species declared as formal & constant"
  else :
    fsp = asp

  formal_species = list(fsp)
  constant_species = list(csp)

  for i in range(len(crn)):
    if crn[i][0] != "reversible" and crn[i][0] != "irreversible":
      crn = crn[:i]
      break
  return (crn, formal_species, constant_species)

def split_reversible_reactions(crn):
  """Replace every reversible reaction with the two corresponding irreversible
  reactions, remove the 'reversible' and 'irreversible' tags.
  """
  new_crn = []
  for [x, r, p] in crn:
    assert (x == 'irreversible' or x == 'reversible')
    new_crn.append([r,p])
    if x == "reversible":
      new_crn.append([p,r])
  return new_crn

def combine_reversible_reactions(crn) : 
  """Condesnse two irreversible reactions into the corresponding reversible
  reactions, add 'reversible' and 'irreversible' tags.
  """
  if type(crn) != list:
    raise RuntimeError("The argument of `rev_reactions' should be a list.")
  new_crn = []
  removed = []
  for r in crn:
    if r in removed:
      continue

    if len(r) == 3 : tag = r[2]
    else : tag = 'irreversible'

    for r2 in crn: 
      if sorted(r[0]) == sorted(r2[1]) and \
          sorted(r[1]) == sorted(r2[0]):
            tag = 'reversible' 
            removed.append(r2)
            break
    r = [r[0],r[1],tag]
    new_crn.append(r)
  return new_crn

def parse_crn_file(filename):
    """Parses the given file and returns the result in the form of
       (crn, formal species, constant species)."""
    crn_document = crn_document_setup()
    crn = crn_document.parseFile(filename, parseAll = True).asList()
    return _post_process(crn)

def parse_crn_string(data):
    """Parses a crn in string format and returns the result in the form of
    (crn, formal_species, constant_species).

    :returns: 
    .. ([['irreversible', ['A', 'B'], ['C','D']], 
    ..   ['irreversible', [], ['A']], 
    ..   ['irreversible', ['X'], []], 
    ..   ['reversible', ['A', 'C'], ['X']], 
    ..   ['A', 'B', 'C', 'D', 'X'], [] ])
    """
    crn_document = crn_document_setup()
    return _post_process(crn_document.parseString(data).asList())
