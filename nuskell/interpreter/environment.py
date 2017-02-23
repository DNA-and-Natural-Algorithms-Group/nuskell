# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2016 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (badelt@caltech.edu)
#
# The interpreter environment for translation schemes.
#

"""The nuskell programming language environment.

Nuskell code is interpreted using public functions of the **Environment**
class.  All other classes and functions are needed internally to set up the
interpreter environment.  Built-in functions and expressions are encapuslated
in the respective namespace. 

.. When modifying this file, use 2-whitespace character indents, 
.. otherwise follow the Google Python Style Guide:
.. http://google.github.io/styleguide/pyguide.html
"""

import sys
from copy import copy

from nuskell.objects import Domain, Complex, TestTube, flatten

class RuntimeError(Exception):
  """Exception for environment.py module.

  Args:
    value (str): Human readable string describing the exception.
  """
  def __init__(self, value):
    self.value = "Runtime error : " + value

  def __str__(self):
    return self.value

class void(object):
  """An empty object returned by the print function"""
  pass

class Species(object):
  """The Species object.

  Args:
    name (str): The name of a species.
  """
  def __init__(self, name):
    self.name = name

class Reaction(object):
  """Format of a single reaction:

  Args:
    reactands (List[str]) : a list of educts, eg. ['A', 'B']
    products (List[str]) : a list of products, eg. ['C']
    reversible (bool) : True or False
  """
  def __init__(self, reactants, products, reversible):
    self.reactants = reactants
    self.products = products
    self.reversible = reversible

class Function(object):
  """A function object of the nuskell language."""
  def __init__(self, args, body):
    self.args = args
    self.body = body

class NusComplex(Complex):
  """Nucleic Acid Complex (Sequence, Structure) pair.

  This is an extension of nuskell.objects.Complex(), which stores an additional
  dictionary of attributes. This dictionary maps domain names as specified in
  the translation-scheme to the unique Domain() object.

  Args:
    attr (Dict[name]=Domain): A mapping between names in the translation scheme
                              and Domain objects 
    sequence (List[str]): Domain Objects
    structure (List[str]): A list of dot-parens characters ".(.+)."
  """
  def __init__(self, attr=dict(), **kwargs):
    super(NusComplex, self).__init__(**kwargs)
    self.attributes = attr

  @property
  def flatten_cplx(self):
    """Resolves nested complexes and their corresponding '~' structure wildcard.

    This function is very specific to Nuskell translation using macros, but it
    might be generalized to handle 'Complex of Complexes' behavior of some sort.
    """
    def resolve(s, c):
      # This is a function with is very specific to nuskell. Sometimes,
      # A complex can contain a complex structure in the sequence, and 
      # a wildcard in the structure.
      if isinstance(s, Domain):
        assert (c in ['(','.',')'])
        return s, c
      elif s == '+' and c == '+' :
        return s, c
      elif isinstance(s, Complex):
        assert (c == '~')
        return s.sequence, s.structure
      elif type(s) == list:
        assert (c == '~')
        # map returns a list of [[seq, str], ..]
        # zip transposes this into tuples ([seq1,seq2,...]),([str1,str2,...])
        return [list(a) for a in zip(*map(lambda i : resolve(s[i],'~'),
          range(len(s))))]
      elif s == '?' :
        assert (c in ['(','.',')'])
        return Domain(list('N' * 15), prefix='h'), c
        #return s, c
      else :
        raise NotImplementedError

    newseq = []
    newstr = []
    for i in range(len(self._sequence)):
      if self._sequence[i] == [] and self._structure[i] == '~' : 
        continue
      se, ss = resolve(self._sequence[i], self._structure[i])
      newseq.append(se)
      newstr.append(ss)
    
    newseq = flatten(newseq)
    newstr = flatten(newstr)

    start = True
    for i in reversed(range(len(newseq))):
      if newseq[i] == '+' and newstr[i]=='+' :
        if start or i == 0 :
          newseq.pop(i)
          newstr.pop(i)
        start = True
      else :
        start = False

    self._sequence = newseq
    self._structure = newstr
    return self._sequence, self._structure

class builtin_expressions(object):
  """ Builtin expressions of Nuskell translation schemes.

  """
  # Say something about special trailer functions here?

  def interpret_expr(self, expr):
    """Recursive interpretation the body of a global variable."""
    operators = { "*" : lambda x, y: x * y,
                  "/" : lambda x, y: x / y,
                  "+" : lambda x, y: x + y,
                  "-" : lambda x, y: x - y,
                  "==" : lambda x, y: x == y,
                  "!=" : lambda x, y: x != y,
                  ">" : lambda x, y: x > y,
                  "<" : lambda x, y: x < y,
                  ">=" : lambda x, y: x >= y,
                  "<=" : lambda x, y: x <= y }

    keywords = {
        'id' : lambda self, content: (self._env, self._ref_binding(content[0])),
        'if' : self._if,  # return self._env, value
        'or' : self._or,  # return self._env, operand1 or operand2
        'and': self._and, # return self._env, operand1 and operand2
        'num': lambda self, content: (self._env, int(content[0])), 
        'quote': lambda self, content: (self._env, content[0]), 
        'dict': self._dict,
        'dna': self._dna, # return self._env, NusComplex()
        'list': self._list, # return self._env, content
        'where': self._where, # return self._env, value
        'uminus': self._uminus, # return self._env, -value (integer only!)
        }

    tag = expr[0]
    content = expr[1:]

    if tag in operators.keys():
      self._env, operand1 = self.interpret_expr(content[0])
      self._env, operand2 = self.interpret_expr(content[1])
      return self._env, operators[tag](operand1, operand2)

    elif tag == 'trailer' :
      trailerkeys = {
          'apply' : self._apply,
          'index' : self._index,
          'attribute' : self._attribute
          }
      # function call, list attribute of an object
      self._env, head = self.interpret_expr(content[0])

      for x in content[1:]:
        key = x[0]
        args = x[1:]
        self._env, head = trailerkeys[key](self, head, args)
      return self._env, head

    elif tag in keywords.keys():
      return keywords[tag](self, content)

    else :
      raise RuntimeError("Unknown expression: `" + tag + "'")
      return self._env, None

  # Context-dependent
  @staticmethod
  def _if(theenv, content) :
    """ Evaluate a (nested) if clause
    
    :param theenv: the environment object
    :param content: an expression that wants to be interpreted
    
    """
    while len(content) > 1:
      theenv._env, test = theenv.interpret_expr(content[0])
      if test:
        theenv._env, value = theenv.interpret_expr(content[1])
        return theenv._env, value
      else:
        content = content[2:]
    theenv._env, value = theenv.interpret_expr(content[0])
    return theenv._env, value

  def _or(self, theenv, content):
    theenv._env, operand1 = theenv.interpret_expr(content[0])
    if operand1:
      return theenv._env, operand1
    theenv._env, operand2 = theenv.interpret_expr(content[1])
    return theenv._env, operand1 or operand2

  def _and(self, theenv, content):
    theenv._env, operand1 = theenv.interpret_expr(content[0])
    if not operand1:
      return theenv._env, False
    theenv._env, operand2 = theenv.interpret_expr(content[1])
    return theenv._env, operand1 and operand2

  def _dict(self, theenv, content):
    kwargs = {}
    for asign in content:
      if asign[1][0] == 'num' : asign[1][1] = int(asign[1][1])
      if asign[1][0] == 'quote' : asign[1][1] = asign[1][1][1:-1]
      kwargs[asign[0][1]] = asign[1][1]
    return theenv._env, kwargs

  def _dna(self, theenv, content):
    domains = content[0]
    dotparen = content[1]
    attributes = {}
    for i in range(len(domains)):
      # TODO history domains
      if domains[i] != "?" and domains[i] != "+":
        starred = (len(domains[i]) == 2)
        dom = domains[i][0]
        # Get the binding for e.g.: ['id', 'd13']
        theenv._env, dom_value = theenv.interpret_expr(dom)        
        attributes[dom[1]] = dom_value
        if starred:
          dom_value = ~dom_value
        domains[i] = dom_value
    return theenv._env, NusComplex(sequence = domains, structure = dotparen,
        attr=attributes)

  def _list(self, theenv, content) :
    for i in range(len(content)):
      theenv._env, content[i] = theenv.interpret_expr(content[i])
    return theenv._env, content

  def _where(self, theenv, content): 
    """ """
    theenv._create_level()
    if len(content) > 1:
      # if this is not a trivial where clause 
      for asgn in content[1]:
        id_list = asgn[0]
        value = asgn[1]

        theenv._env, value = theenv.interpret_expr(value)
        id_list = self._bfunc_Obj.remove_id_tags(id_list)

        for key, value in self._bfunc_Obj.asgn_pattern_match(id_list, value) :
          theenv._env = theenv._create_binding(key, value)

    theenv._env, value = theenv.interpret_expr(content[0])
    theenv._destroy_level()
    return theenv._env, value

  def _uminus(self, theenv, content): 
    theenv._env, value = theenv.interpret_expr(content[0])
    if type(value) != int:
      raise RuntimeError(
          "The unary minus operator can only be used with integers.")
    return theenv._env, -value

  # trailer functions
  def _apply(self, theenv, head, args):
    for i in range(len(args)):
      theenv._env, args[i] = theenv.interpret_expr(args[i])
    theenv._env, head = theenv._eval_func(head, args)
    return theenv._env, head

  def _index(self, theenv, head, args):
    theenv._env, subscript = theenv.interpret_expr(args[0])
    if type(head) != list:
      raise RuntimeError("Only lists can be indexed.")
    if type(subscript) != int:
      raise RuntimeError("Subscript should be an integer.")
    head = head[subscript]
    return theenv._env, head

  def _attribute(self, theenv, head, args):
    identifier = args[0][1] # strip the tag
    if isinstance(head, NusComplex):
      head = head.attributes[identifier]
    elif identifier in head.__dict__:
      head = head.__dict__[identifier]
    else:
      raise RuntimeError(
          "The attribute `"+identifier+"' could not be found.")
    return theenv._env, head

class builtin_functions(object):
  """Builtin functions of Nuskell translation schemes.

  This object is typically initialized within the environment after the
  respective functions have been bound.  Most methods do not require any class
  variables and are therefore declared as staticmethods.  As a consequence, the
  functions can be accessed without prior initialization of the Object. 
  """
  ################################
  ### Builtin-function section ###
  ################################

  def __init__(self, sdl=6, ldl=15) :
    self._short_domain_length = sdl
    self._long_domain_length = ldl

  def eval_builtin_functions(self, f, args):
    """Evaluate built-in functions. These functions do not alter the
    environment and they are contained in their own class for readability.
    """

    keywords = {
        'print': self._print,
        'abort': self._abort,
        'tail' : self.tail,
        'flip' : self.flip,
        'long' : self.long,
        'short' : self.short,
        'infty' : self.infty,
        #'unique' : self.unique,
        'complement' : self.complement,
        'rev_reactions' : self.rev_reactions,
        'irrev_reactions' : self.irrev_reactions,
        }

    if f in keywords.keys() :
      return keywords[f](args)
    else :
      raise RuntimeError("`" + f + "' could not be resolved.")
      return None

  def long(self, args):
    """A function that returns branch-migration domains. """
    if args :
      kwargs = args[0]
    else :
      kwargs = {'len' : self._long_domain_length}

    if 'len' not in kwargs and 'sense' not in kwargs:
      kwargs['len'] = self._long_domain_length

    # Nuskell Defaults: 
    tag = kwargs['tag'] if 'tag' in kwargs else 'd'
    uppercon = kwargs['sense'] if 'sense' in kwargs else 'H' * kwargs['len']
    lowercon = kwargs['antis'] if 'antis' in kwargs else 'D' * kwargs['len']

    upper = Domain(list(uppercon), prefix=tag)
    lower = upper.get_ComplementDomain(list(lowercon))
    return upper

  def short(self, args):
    """A function that returns toehold domains.
 
    # supported keywords: 'len', 'sense', 'antisense', 'tag'
    # len is short for standard 'sense & antisense'
    # tag supports toehold and branch-migration which is
    # currently equal to short() and long()
   
    """
    if args :
      kwargs = args[0]
    else :
      kwargs = {'len' : self._short_domain_length}

    if 'len' not in kwargs and 'sense' not in kwargs:
      kwargs['len'] = self._short_domain_length

    # Nuskell Defaults: 
    tag = kwargs['tag'] if 'tag' in kwargs else 't'
    uppercon = kwargs['sense'] if 'sense' in kwargs else 'H' * kwargs['len']
    lowercon = kwargs['antis'] if 'antis' in kwargs else 'D' * kwargs['len']

    upper = Domain(list(uppercon), prefix=tag)
    lower = upper.get_ComplementDomain(list(lowercon))
    return upper

  def _print(self, args):
    """Print statment, primarily to debug nuskell scripts"""
    for a in args:
      if isinstance(a, Complex) :
        print map(str, a.sequence)
      else :
        print str(a),
    print
    return void()

  def _abort(self, args):
    """Raise SystemExit and print the error message"""
    raise SystemExit('EXIT: ' + args[0])

  @staticmethod
  def tail(args) :
    """ Returns the list minus the first element.. """
    if type(args[0]) != list:
      raise RuntimeError("`tail' should have a list as its argument.")
    return args[0][1:]

  @staticmethod
  def flip(args) :
    """A matrix transpose function for lists of lists. 

    :param args: a list that contains both the lol and and integer to specify
    the dimenstion.

    for example:
    input: [[[x,y,z],[a,b,c],[n,l,k],[u,v,w]],3]
    returns: [[x,a,n,u],[y,b,l,v],[z,c,k,w]]
    """
    if type(args[0]) != list:
      raise RuntimeError("The first argument of `flip' should be a list.")
    if type(args[1]) != int:
      raise RuntimeError("The second argument of `flip' should be an integer.")

    l = args[0]
    n = args[1]
    res = [0] * n
    for i in range(n):
      res[i] = [0] * len(l)
    for i in range(len(l)):
      if len(l[i]) != n:
        raise RuntimeError(
            "The argument of `flip' does not have a required form.")
      for j in range(n):
        res[j][i] = l[i][j]
    return res

  @staticmethod
  def infty(args) :
    """
    'infty' is a built-in function that takes a Complex instance and puts
    that the species is supposed to have "infinite" concentration. The
    resulting object will be a TestTube instance
    """
    if not isinstance(args[0], Complex):
      raise RuntimeError("The argument of `infty' should be a Complex")
    args[0].flatten_cplx
    if args[0].sequence == [] and args[0].structure == [] :
      return TestTube()
    else :
      return TestTube(complexes={args[0].name: (args[0], float("inf"), True)})
      #return TestTube(complexes={args[0].name: (args[0], float("inf"))})

    # Translate NusCompex -> Complex, so that the environment returns a 
    # Standard TestTube Object.
    #standard = Complex(sequence=args[0].sequence, structure=args[0].structure)
    #return TestTube(complexes={standard.name: standard})

  @staticmethod
  def const_nM(args) :
    if not isinstance(args[0], Complex):
      raise RuntimeError("The first argument of `const_nM' should be a Complex")
    if type(args[1]) != float:
      raise RuntimeError("The second argument of `const_nM' should be a concentration <float>")
    args[0].flatten_cplx
    if args[0].sequence == [] and args[0].structure == [] :
      return TestTube()
    else :
      return TestTube(complexes={args[0].name: (args[0], args[1])})
 
  #@staticmethod
  #def unique(args) :
  #  if type(args[0]) != int:
  #    raise RuntimeError("The first argument of `unique' should be an integer.")
  #  return Domain(args[0])

  @staticmethod
  def complement(args) :
    """ Returns the complement of a given input. It can interpret Complex, 
    Domains and Dot-bracket lists, as well as lists of lists.

    :param args: an expression in form of an array, where the first element
    in the array contains the data of interest. ... yeah, a bit complicated...

    :return: complement of the input
    """
    x = args[0]
    if type(x) == list:
      # args[0] forces us to introduce additional lists ...
      for i in range(len(x)) : x[i] = [x[i]]
      return reversed(map(self._complement, x))
    elif isinstance(x, Domain):
      print 'Untested:', ~x
      return ~x
    elif isinstance(x, Complex):
      raise NotImplementedError
      return Complex(
          self._complement([x.sequence]),
          self._complement([x.structure]),
          dict(x.attributes))
    elif x == "(":
      return ")"
    elif x == ")":
      return "("
    else:
      return x
  
  @staticmethod
  def rev_reactions(args) : 
    if type(args[0]) != list:
      raise RuntimeError("The argument of `rev_reactions' should be a list.")
    crn = args[0]
    new_crn = []
    removed = []
    for r in crn:
      if r in removed:
        continue
      reversible = r.reversible
      for r2 in crn: 
        if sorted(r.reactants) == sorted(r2.products) and \
            sorted(r.products) == sorted(r2.reactants):
              reversible = True
              removed.append(r2)
              break
      r = Reaction(r.reactants, r.products, reversible)
      new_crn.append(r)
    return new_crn

  @staticmethod
  def irrev_reactions(args) : 
    """
    a function that divides every reversible reaction into a pair of
    irreversible reaction.
    """
    if type(args[0]) != list:
      raise RuntimeError("The argument of `irrev_reactions' should be a list.")
    crn = args[0]
    new_crn = []
    for r in crn:
      if r.reversible:
        new_crn.append(Reaction(r.reactants, r.products, False))
        new_crn.append(Reaction(r.products, r.reactants, False))
      else:
        new_crn.append(r) 
    return new_crn

  @staticmethod
  def asgn_pattern_match(id, value):
    """Does the pattern matching for list assignments.
       ex) [a, b, c] = [1, 2, 3]
    """
    if type(id) == list:
      if type(value) != list or len(id) != len(value):
        raise RuntimeError(
            "Pattern matching failed while assigning values to " + str(id) + ".")
      result = []
      for i in range(len(id)):
        result += builtin_functions.asgn_pattern_match(id[i], value[i])
      return result
    elif type(id) == str:
      return [(id, value)]
    else:
      raise RuntimeError("""The left hand side of an assignment should be 
      either an identifier or a list-tree consisting only of identifiers.""")

  @staticmethod
  def remove_id_tags(l):
    """Helper function to remove all tags from the given id_list."""
    kwd = l[0]
    if kwd == "idlist":
      return map(builtin_functions.remove_id_tags, l[1:])
    elif kwd == "id":
      return l[1]
    else:
      raise RuntimeError("""The left hand side of an assignment should be 
      either an identifier or a list-tree consisting only of identifiers.""")
  
class Environment(builtin_expressions):
  """The Nuskell environment to interpret translation schemes. 
  
  Environment inherits all builtin expressions. It provides the interface to
  interpret code in a translation scheme, and to execute the code by assigning
  formal species and applying them to the main() function.
  """

  def __init__(self, name='default_env', sdlen = 6, ldlen = 15):
    """ Initialize the environment. 

    Collects optional arguments for the setup of the Environment and
    initializes built-in functions: short(), long(), tail(), flip(), etc...
    
    Args:
      name (optional, str): The name of the Environment.
      sdlen (optional, int): Dlfault length of built-in short() domains
      ldlen (optional, int): Default length of built-in long() domains

    Note:
      short() and long() domain lengths can also be changed at each assignment,
      e.g. short(7) is a toehold of length 7, long(10) is a branch-migration
      domain of length 10.
    """
    self._name = name
    self._env = [{}]

    # Setup the builtin functions. 
    self._env, self._bfunc_Obj = self._init_builtin_functions(sdlen, ldlen)

  # Public functions #
  def interpret(self, code):
    """ Returns the Environment (the final namespace).
    
    Creates bindings for variables and functions defined in the body of the
    translation scheme. 

    Note: 
      At this point all keywords (class, function, macro, module) are treated
      exactly the same, only the **global** keyword is special, because the 
      global expressions are **interpreted first** and then bound.

    Returns:
      self.env: An updated Environment
    """
    for stmt in code:
      kwd = stmt[0]
      body = stmt[1:]

      if kwd == "global":
        # create binding to interpretation of the expression
        id_list = body[0]
        id_list = self._bfunc_Obj.remove_id_tags(id_list)

        value = body[1]
        self._env, value = self.interpret_expr(value)

        for key, value in self._bfunc_Obj.asgn_pattern_match(id_list, value) :
          self._create_binding(key, value)
      else:
        # create binding to the function, without interpretation
        # e.g. kwd = module; id = 'rxn'; args = ['r']; body_ = [where...]
        assert body[0][0] == "id"
        id = body[0][1]
        args = map(lambda x: x[1], body[1]) # remove 'id' tags from args
        body_ = body[2] # the ['where' [...]] part
        self._create_binding(id, Function(args, body_))
    return self._env

  def translate_formal_species(self, fs_list):
    """ Apply the formal() function to the formal species in the input CRN.

    First, the bindings for formal species are created, then the formal()
    function is applied to every formal species in the input CRN.

    Returns:
      self.formal_species_dict : A dictionary of {fs:NusComplex()}.
    """
    formal_species_objects = map(Species, fs_list)

    # compile the formal species
    self._create_binding("__formalspecies__", formal_species_objects)

    # map(formal, __formalspecies__)
    self._env, fs_result = self.interpret_expr(["trailer",
      ["id", "map"], ["apply", 
        ["id", "formal"], 
        ["id", "__formalspecies__"]]])

    self.formal_species_dict = {}
    for i in range(len(fs_list)):
      self.formal_species_dict[fs_list[i]] = fs_result[i]

    return self.formal_species_dict

  def translate_constant_species(self, cs_list, crn_parsed):
    """Implementation CRN translator function.

    This function builds a TestTube object from species that have been declared
    as 'constant' in the input CRN. This enables to specify an implementation 
    CRN directly as input and use the translation scheme merely to describe how
    formal and constant species should look like. 

    TODO: This function might be combined with the regular formalCRN-to-DSD
    approach, but we should think about when this makes sense.

    """
    constant_species_objects = map(Species, cs_list)

    # compile the formal species
    self._create_binding("__constantspecies__", constant_species_objects)

    self._env, cs_result = self.interpret_expr(["trailer",
      ["id", "map"], ["apply", 
        ["id", "constant"], 
        ["id", "__constantspecies__"]]])

    self.constant_species_dict = {}
    for i in range(len(cs_list)):
      self.constant_species_dict[cs_list[i]] = cs_result[i]

    # replace every instance of a constant species in the CRN with the
    # respective cs_result instance. Skip all non-constant species
    crn_remap = []
    for r in crn_parsed:
      react = []
      for e, field in enumerate(r):
        if e == 0 :
          react.append(field)
        else :
          spec = []
          for s in field :
            if s in self.constant_species_dict :
              spec.append(self.constant_species_dict[s])
          react.append(spec)
      crn_remap.append(react)

    crn_object = map(
        lambda x: Reaction(x[1], x[2], x[0] == "reversible"), crn_remap)

    self._create_binding("__crn__", crn_object)
    self._env, self.constant_species_solution = self.interpret_expr( 
        ["trailer", ["id", "main"], ["apply", ["id", "__crn__"]]])
 
    return self.constant_species_solution

  def translate_reactions(self, crn_parsed):
    """Execute the main() function of the translation scheme.
    
    The input CRN is replaced with previously initialized formal species objects.

    Args:
      crn_paresed (List[Lists]): A crn in crn_parser format.

    Returns:
      self.main_solution (Solution) : The Solution object with all constant
      species

    Raises:
      RuntimeError: If the compiled formal species cannot be found

    """
    if not self.formal_species_dict :
      raise RuntimeError('Could not find the compiled formal species!')

    # replace every fs (str) with fs(NusComplex())
    crn_remap = map(
        lambda x: [x[2]] + map( lambda y: map(
            lambda z: self.formal_species_dict[z], y), x[:2]), crn_parsed)

    crn_object = map(
        lambda x: Reaction(x[1], x[2], x[0] == "reversible"), crn_remap)

    # main(__crn__)
    self._create_binding("__crn__", crn_object)
    self._env, self.constant_species_solution = self.interpret_expr( 
        ["trailer", ["id", "main"], ["apply", ["id", "__crn__"]]])

    return self.constant_species_solution

  # Private Functions #
  def _init_builtin_functions(self, sdlen, ldlen):
    self._create_binding("print", Function(["s"], "print"))
    self._create_binding("abort", Function(["s"], "abort"))
    self._create_binding("tail", Function(["l"], "tail"))
    self._create_binding("flip", Function(["l", "n"], "flip"))
    self._create_binding("long", Function([], "long"))
    self._create_binding("short", Function([], "short"))
    self._create_binding("infty", Function(["species"], "infty"))
    #self._create_binding("unique", Function(["l"], "unique"))
    #self._create_binding("complement", Function(["l"], "complement"))
    self._create_binding("rev_reactions", Function(["crn"], "rev_reactions"))
    self._create_binding("irrev_reactions", Function(["crn"], "irrev_reactions"))
    self._create_binding("empty", TestTube())
    return self._env, builtin_functions(sdlen, ldlen)

  def _create_binding(self, name, value):
    """ Create binding of a Nuskell function.

    Args:
      name (str): Name of the function.
      value (...): A built-in data type (Function, TestTube, Species,
        Domain, Reaction, Complex, etc ...), either as single value or list

    Returns: 
      An updated environment inclding the function binding in the top-level.
    """

    bindings = (Function, TestTube, Species, Domain, Reaction, NusComplex, 
        void, int, list)
    if isinstance(value, list) :
      assert all(isinstance(s, bindings) for s in value)
    else :
      assert isinstance(value, bindings)

    self._env[-1][name] = value
    return self._env

  # Private environment modification functions #
  def _create_level(self): 
    # Create a new level for function bindings.
    # It is commonly triggered by a 'where' statement or during the evaluation
    # of a non-built-in function call.

    # interpret_expr -> _where
    # interpret_expr -> _trailer -> _apply -> _eval_func
    self._env.append({})
    return self._env
  
  def _destroy_level(self):
    # Revert to the old level for function bindings.
    # It is commonly triggered by a 'where' statement or after the evaluation
    # of a non-built-in function call.

    # interpret_expr -> _where
    # interpret_expr -> _trailer -> _apply -> _eval_func
    self._env.pop()
    return self._env

  def _ref_binding(self, name):
    # Search levels from last to first to find a reference to an exisiting
    # function binding.

    # Args: 
    #   name (str) : Name of a function

    # Return: 
    #   Function binding (self._env[?][name])
    for level in reversed(self._env):
      if name in level.keys():
        return level[name]
    raise RuntimeError("Cannot find a binding for `" + name + "'.")
    return None

  def _eval_func(self, f, args):
    if type(f.body) == str: # the function is a built-in function
      return self._env, self._bfunc_Obj.eval_builtin_functions(f.body, args)
  
    def hardcopy_list(l):
      if type(l) != list:
        return copy(l)
      return map(hardcopy_list, l)
  
    self._create_level()
    if not isinstance(f, Function):
      raise RuntimeError(str(f) + "is not a function.")
  
    if len(f.args) > len(args):
      raise RuntimeError("The function `" + f.name + "' requires at least " \
          + str(len(f.args)) + " arguments but only found " + str(args) +\
          " arguments.")
  
    for i in range(len(f.args)):
      arg_name = f.args[i]
      arg_value = args[i]
      self._env = self._create_binding(arg_name, arg_value)
  
    # hardcopy the function body so that it does not change during evaluation.
    self._env, value = self.interpret_expr(hardcopy_list(f.body))
    self._destroy_level()
    return self._env, value


