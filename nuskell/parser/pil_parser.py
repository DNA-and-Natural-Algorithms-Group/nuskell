#
#
# Copyright (c) 2017 Caltech. All rights reserved.
# Written by Stefan Badelt (badelt@caltech.edu).
#
#
# Parser module for pil files in kernel notation.
#

from pyparsing import (Word, Literal, Group, Suppress, Optional, ZeroOrMore, Combine, White,
    OneOrMore, alphas, alphanums, nums, delimitedList, StringStart, StringEnd, Forward, 
    LineEnd, pythonStyleComment, ParseElementEnhance)
# from pyparsing import Group, Forward, Word, Combine, Literal, Optional, Suppress, ZeroOrMore, OneOrMore, StringEnd, delimitedList, nestedExpr, alphanums

def pil_document_setup():
  crn_DWC = "".join(
      [x for x in ParseElementEnhance.DEFAULT_WHITE_CHARS if x != "\n"])
  ParseElementEnhance.setDefaultWhitespaceChars(crn_DWC)

  # letter = (* any alphanumeric character *)
  # identifier = letter, {letter}, [*]
  # pattern = expression, {space, expression}
  # expression = domain | loop | wildcard | break
  # domain = identifier
  # loop = identifier, "(", [pattern] ,")"
  # wildcard = "?" | "_" | "~" (* ... etc.*)
  # break = "+"
 
  def T(x, tag):
    def TPA(tag):
      return lambda s, l, t: [tag] + t.asList()
    return x.setParseAction(TPA(tag))

  W = Word
  G = Group
  S = Suppress
  O = Optional
  C = Combine
  L = Literal
 
  # NOTE: Exchange comment to forbid/allow names starting with digits
  #identifier = W(alphanums + "_-")
  identifier = W(alphas, alphanums + "_-")
  number = W(nums, nums)
  gorf = C(W(nums) + O((L('.') + W(nums)) | (L('e') + O('-') + W(nums))))
  domain = G(T(S("length") + identifier + S("=") + number + OneOrMore(LineEnd().suppress()),'domain'))

  # NOTE: exchange the comment for asense if you want to allow input in form of "x( ... y)", 
  # but also double-check if that really works...
  sense  = Combine(identifier + O(L("^")) + O(L("*")))
  #asense = (Combine(sense + S(")")) | S(")"))
  asense = S(")")

  sbreak = L("+")

  pattern = Forward()
  # NOTE: Remove S(White()) for backward compatiblility: )) is not allowed anymore.
  #loop = (Combine(sense + S("(")) + O(pattern) + asense)
  innerloop = S(White()) + pattern + S(White()) | G(S(White()))
  loop = (Combine(sense + S("(")) + innerloop + asense)
  pattern << G(OneOrMore(loop | sbreak | sense))

  unit = L('M') | L('mM') | L('uM') | L('nM') | L('pM')
  conc = G(S('@') + L('initial') + gorf + unit) | G(S('@') + L('constant') + gorf + unit)

  cplx = G(T(identifier + S("=") + OneOrMore(pattern) + O(conc) + OneOrMore(LineEnd().suppress()),'complex'))

  stmt = domain | cplx

  document = StringStart() + ZeroOrMore(LineEnd().suppress()) + OneOrMore(stmt) + StringEnd()
  document.ignore(pythonStyleComment)

  return document

def parse_pil_file(data):
  document = pil_document_setup()
  return document.parseFile(data).asList()

def parse_pil_string(data):
  document = pil_document_setup()
  return document.parseString(data).asList()
