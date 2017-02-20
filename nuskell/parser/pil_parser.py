#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Parser module for domain specification description files (*.dom).
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
  L = Literal
 
  identifier = W(alphas, alphanums + "_-")
  number = W(nums, nums)
  domain = G(T(S("length") + identifier + S("=") + number + OneOrMore(LineEnd().suppress()),'domain'))

  # NOTE: exchange the comment for asense if you want to allow input in form of "x( ... y)", 
  # but also double-check if that really works...
  sense  = Combine(identifier + O(L("^")) + O(L("*")))
  #asense = (Combine(sense + S(")")) | S(")"))
  asense = S(")")

  sbreak = L("+")

  pattern = Forward()
  # NOTE: Remove S(White()) for backward compatiblility: )) is not allowed anymore.
  loop = (Combine(sense + S("(")) + S(White()) + O(pattern) + S(White()) + asense)

  expression = (loop | sbreak | sense)
  pattern << G(OneOrMore(expression))

  cplx = G(T(identifier + S("=") + OneOrMore(pattern) + OneOrMore(LineEnd().suppress()),'complex'))

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
