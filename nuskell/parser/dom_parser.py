#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Parser module for domain specification description files (*.dom).
#

from pyparsing import (Word, Literal, Group, Suppress, Optional, ZeroOrMore,
    OneOrMore, alphas, alphanums, nums, delimitedList, StringStart, StringEnd,
    LineEnd, pythonStyleComment, ParseElementEnhance)

def dom_document_setup():
  crn_DWC = "".join(
      [x for x in ParseElementEnhance.DEFAULT_WHITE_CHARS if x != "\n"])
  ParseElementEnhance.setDefaultWhitespaceChars(crn_DWC)

  W = Word
  G = Group
  S = Suppress
  O = Optional
  L = Literal
  
  identifier = W(alphas, alphanums + "_")
  number = W(nums, nums)
  domain = G(OneOrMore((G(identifier + O("*")) | "?" | "+")))
  dotparen = G(OneOrMore(W("().+", max=1)))
  structure = domain + OneOrMore(LineEnd().suppress()) + dotparen
  sequence = G(S("sequence") + identifier + S(":") + number + OneOrMore(LineEnd().suppress()))
  molecule = G(identifier + S(":") + OneOrMore(LineEnd().suppress()) + structure)
  
  document = StringStart() + ZeroOrMore(LineEnd().suppress()) + G(ZeroOrMore(sequence)) + G(OneOrMore(molecule + OneOrMore(LineEnd().suppress()))) + StringEnd()
  document.ignore(pythonStyleComment)
  return document

def parse_dom_file(data):
  document = dom_document_setup()
  return document.parseFile(data).asList()

def parse_dom_string(data):
  document = dom_document_setup()
  return document.parseString(data).asList()
