#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# Parser module for translation scheme description files (*.ts).
#

from pyparsing import (Word, Literal, Group, Suppress, Optional, Forward,
    OneOrMore, ZeroOrMore, nums, alphas, alphanums, delimitedList,
    operatorPrecedence, ParserElement, opAssoc, StringStart, StringEnd,
    pythonStyleComment, quotedString, ParseElementEnhance)

"""
grammar:

  document ::= statements | statement 

  statement ::= class | function | module | macro | global
  class     ::= 'class'    + identifier + '(' + identifiers + ')' + '=' + expr
  function  ::= 'function' + identifier + '(' + identifiers + ')' + '=' + expr
  module    ::= 'module'   + identifier + '(' + identifiers + ')' + '=' + expr
  macro     ::= 'macro'    + identifier + '(' + identifiers + ')' + '=' + expr
  global    ::= 'global'   + id_list + '=' + expr

  id_list ::= '[' + id_list + ']' | identifier

  expr ::= if_expr | where_expr | quote_expr

  if_expr ::= 'if' + expr + 'then' + expr + 
              'elsif' + expr + 'then' + expr + 'else' + expr

  where_expr ::= test + where_clause | test
  test ::= 
  where_clause ::= 'where' + '{' + asigns + '}' | 'where' + asgn

  args ::= 
  
"""

def ts_document_setup():
  """The gramar to parse a translation scheme."""
  ParserElement.enablePackrat()
  ParseElementEnhance.setDefaultWhitespaceChars(" \n\t\r")
  
  W = Word
  G = Group
  S = Suppress
  O = Optional
  L = Literal
  
  def TPAOp(s, l, t):
    """ """
    def helper(t):
      if len(t) == 1:
        return t[0]
      else:
        return [t[1], t[0], helper(t[2:])]
    return [helper((t.asList())[0])]
  
  def T(x, tag):
    def TPA(tag):
      return lambda s, l, t: [tag] + t.asList()
    return x.setParseAction(TPA(tag))
  
  expr = Forward()
  exprlist = delimitedList(expr) 

  paren = S("(") + expr + S(")")
  list  = G(T(S("[") + O(exprlist) + S("]"), "list"))
  identifier = G(T(W(alphas, alphanums + "_"), "id"))
  number = G(T(W(nums), "num"))
  domains = G(S('"') + \
              OneOrMore(G(identifier + O("*")) | "?" | "+") + \
              S('"'))
  dotparen = G(S('"') + OneOrMore(W("().~+", max=1)) + S('"'))
  dna = G(T(domains + S("|") + dotparen, "dna"))
  
  atom = paren | list | identifier | number | dna
  
  function_call = G(T(S("(") + O(exprlist) + S(")"), "apply"))
  index = G(T(S("[") + expr + S("]"), "index"))
  attribute = G(T(S(".") + identifier, "attribute"))
  trailer = function_call | index | attribute
  atom_trailers = G(T(atom + ZeroOrMore(trailer), "trailer"))
  
  factor = Forward()
  factor << (G(T(S("-") + factor, "uminus")) | atom_trailers)
  
  test = operatorPrecedence(factor,
      [(W("*/",max=1), 2, opAssoc.LEFT, TPAOp),
       (W("+-",max=1), 2, opAssoc.LEFT, TPAOp),
       ((L("==") | L(">=") | L("<=") | L(">") | L("<") | L("!=")), \
           2, opAssoc.LEFT, TPAOp),
       ((L("and") | L("or")), 2, opAssoc.LEFT, TPAOp)])

  
  id_list = Forward()
  id_list << (G(T(S("[") + delimitedList(id_list) + S("]"), "idlist")) \
           | identifier)
  
  asgn = G(id_list + S("=") + expr)
  where_clause = G(S("where") + S("{") + \
      delimitedList(asgn, ";") + S("}")) | G(S("where") + asgn)

  quote_expr = G(T(quotedString, 'quote'))
  dict_expr = G(T(OneOrMore(O(S(',')) + G(identifier + S(':') + (quote_expr| number))), 'dict')) 

  where_expr = G(T(test + O(where_clause), "where"))
  if_expr = G("if" + expr + S("then") + expr + \
      ZeroOrMore(S("elseif") + expr + S("then") + expr) + \
      S("else") + expr)
  
  expr << (dict_expr | if_expr | where_expr | quote_expr )
  
  class_def    = G("class" + identifier + \
      G(S("(") + O(delimitedList(identifier)) + S(")")) + S("=") + expr)
  function_def = G("function" + identifier + \
      G(S("(") + O(delimitedList(identifier)) + S(")")) + S("=") + expr)
  module_def   = G("module" + identifier + \
      G(S("(") + O(delimitedList(identifier)) + S(")")) + S("=") + expr)
  macro_def    = G("macro" + identifier + \
      G(S("(") + O(delimitedList(identifier)) + S(")")) + S("=") + expr)
  global_def   = G("global" + id_list + S("=") + expr)
  
  stmt = class_def | function_def | module_def | macro_def | global_def
  
  document = StringStart() + O(delimitedList(stmt, ";")) + StringEnd()
  document.ignore(pythonStyleComment)
  return document

def parse_ts_file(filename):
  """Parses the given .ts file and returns the result as a list."""
  ts_document = ts_document_setup()
  return ts_document.parseFile(filename, parseAll = True).asList()

def parse_ts_string(data):
    """Parses the given string and returns the result as a list."""
    ts_document = ts_document_setup()
    return ts_document.parseString(data).asList()
