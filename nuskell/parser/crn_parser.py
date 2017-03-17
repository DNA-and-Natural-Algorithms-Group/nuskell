#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Parser module for chemical reaction network description files (*.crn).
#

from pyparsing import (Word, Literal, Group, Suppress, Optional, ZeroOrMore, Combine, nums,
    alphas, alphanums, delimitedList, StringStart, StringEnd, LineEnd, srange, OneOrMore,
    pythonStyleComment, ParseElementEnhance)

def crn_document_setup(modular = False):
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
      return lambda s, l, t: t.asList() + [tag]
    return x.setParseAction(TPA(tag))
  
  crn_DWC = "".join(
      [x for x in ParseElementEnhance.DEFAULT_WHITE_CHARS if x != "\n"])
  ParseElementEnhance.setDefaultWhitespaceChars(crn_DWC)
  
  identifier = W(alphas, alphanums+"_")

  multiplier = W(nums)
  species = G(O(multiplier) + identifier)

  rate = C(W(nums) + O((L('.') + W(nums)) | (L('e') + O('-') + W(nums))))

  k = G(S('[') + S('k') + S('=') + rate + S(']'))
  rev_k = G(S('[') + S('kf') + S('=') + rate + S(',') + \
                     S('kr') + S('=') + rate + S(']'))

  reaction = T(G(O(delimitedList(species, "+"))) + \
             S("->") + \
             G(O(delimitedList(species, "+"))) + O(k), 'irreversible')

  rev_reaction = T(G(O(delimitedList(species, "+"))) + \
                 S("<=>") + \
                 G(O(delimitedList(species, "+"))) + O(rev_k), 'reversible')
  

  expr = G(reaction | rev_reaction) 

  if modular :
    module = G(expr + ZeroOrMore(S(";") + expr))
  else :
    module = expr + ZeroOrMore(S(";") + expr)

  formal = G(O(S(";")) + L("formal") + S(L("=") + \
             L("{")) + O(delimitedList(identifier)) + S("}"))

  constant = G(O(S(";")) + L("constant") + S(L("=") + \
               L("{")) + O(delimitedList(identifier)) + S("}"))
  
  crn = OneOrMore(module + ZeroOrMore(S(LineEnd()))) + O(formal) + ZeroOrMore(S(LineEnd())) + O(constant) + ZeroOrMore(S(LineEnd()))

  document = StringStart() + ZeroOrMore(S(LineEnd())) + crn + StringEnd()
  document.ignore(pythonStyleComment)
  return document

def _post_process(crn):
  """Take a crn and return it together with a list of formal and
  constant species.
  """
  def remove_multipliers(species) :
    flat = []
    for s in species:
      if len(s) == 1:
        flat.append(s[0])
      elif len(s) == 2:
        ss = [s[1]] * int(s[0])
        flat.extend(ss)
    return flat

  new = []
  asp = set()
  fsp = set()
  csp = set()
  fsp_manual = False
  csp_manual = False
  for line in crn:
    if line[0] == "formal":
      fsp = fsp.union(line[1:])
      fsp_manual = True
    elif line[0] == "constant":
      csp = csp.union(line[1:])
      csp_manual = True
    elif len(line) == 3:
      # No rate specified
      r, p, t = line
      r = remove_multipliers(r)
      p = remove_multipliers(p)
      if t == 'reversible':
        new.append([r,p,[None, None]])
      elif t == 'irreversible':
        new.append([r,p,[None]])
      else :
        raise ValueError('Wrong CRN format!', line)
    elif len(line) == 4:
      r, p, k, t = line
      r = remove_multipliers(r)
      p = remove_multipliers(p)
      if t == 'reversible':
        assert len(k) == 2
        new.append([r,p,k])
      elif t == 'irreversible':
        assert len(k) == 1
        new.append([r,p,k])
      else :
        raise ValueError('Wrong CRN format!', line)
    else :
      raise ValueError('Wrong CRN format!', line)
    asp = asp.union(r).union(p)
  crn = new

  # Check that formal/constant assignments make sense!
  if fsp_manual and csp_manual :
    # NOTE: Also empty sets can be forced!
    pass
  elif csp_manual :
    fsp = asp - csp
  elif fsp_manual :
    csp = asp - fsp
  else :
    fsp = asp

  if fsp & csp :
    raise ValueError, "species declared as formal & constant"

  return crn, sorted(list(fsp)), sorted(list(csp))

def split_reversible_reactions(crn):
  """ 
  Replace every reversible reaction with the two corresponding irreversible
  reactions.
  """
  new = []
  for [r, p, k] in crn :
    #if None in k :
    #  print Warning('# Set missing rates to 1.')
    #  k[:] = [x if x != None else 1 for x in k]

    if len(k) == 2:
      new.append([r,p,[k[0]]])
      new.append([p,r,[k[1]]])
    else :
      new.append([r,p,k])
  return new

def combine_reversible_reactions(crn) : 
  """Condense two irreversible reactions into the corresponding reversible reactions. """
  new_crn = []
  removed = []
  for rxn in crn:
    if rxn in removed:
      continue
    [r, p, k] = rxn
    assert type(r) == list and type(p) == list and type(k) == list

    for rxn2 in crn: 
      #if rxn in removed:
      #  continue
      [r2, p2, k2] = rxn2
      if sorted(r) == sorted(p2) and sorted(p) == sorted(r2):
        if len(k) == 2 or len(k2) == 2 :
          raise ValueError('reaction specified twice!')
        else :
          removed.append(rxn2)
          k += k2
        break
    new_crn.append([r, p, k])
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
    .. ([[['A', 'B'], ['C','D'], 'irreversible'],
    ..   [[], ['A'], 'irreversible'], 
    ..   [['X'], [], 'irreversible'], 
    ..   [['A', 'C'], ['X'], 'reversible'], 
    ..   ['A', 'B', 'C', 'D', 'X'], [] ])
    """
    crn_document = crn_document_setup()
    return _post_process(crn_document.parseString(data).asList())
