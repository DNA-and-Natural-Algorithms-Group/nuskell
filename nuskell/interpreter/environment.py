# -*- coding: utf-8 -*-
#
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
# Modified by Stefan Badelt (badelt@caltech.edu)
#
#

"""The DSD translation scheme language environment.

This module comprises all classes and functions that are used to set up the
nuskell compiler environment. The environment has multiple levels. Built-in
functions are defined in the base-level upon initialization, additional nuskell
code is embedded using public functions of the **Environment** class. At this
point, only the interpreter module communicates with the environment module in
order to set up and execute the nuskell script together with the CRN of
interest.

.. When modifying this file, use 2-whitespace character indents, 
.. otherwise follow the Google Python Style Guide:
.. http://google.github.io/styleguide/pyguide.html
"""

import sys
import nuskell.include.dnaobjects as dnaobjects
from copy import copy

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

class NusDomain(dnaobjects.Domain):
  """The Nuskell Domain class, based on DNAObjects
  
  A DNA domain is a sequence of consecutive nucleotides. Typically, a domain,
  (as well as its complementary domain) is present in multiple DNA strands in a
  DSD circuit.  The domain-ID is a unique descriptor refering to all DNA
  domains with that ID.

  Args: 
    kwargs from DNAObects.Domain
    tag (string): Metadata at the domain-level, e.g. 'toehold',
                  'branch-migration', or 'history'

  """
 
  def __init__(self, domaintag='', **kwargs):
    """ 
    Args: 
      domaintag: You can specify different types of domains for 
        your sequence designer: 'toehold', 'branch-migration', etc...
    """
    super(NusDomain, self).__init__(**kwargs)
    self.tag=domaintag

    if 'bottom' in kwargs:
      self.bottom = kwargs['bottom']
    else :
      self.bottom = 'D' * len(self.constraints)

    if len(self.bottom) != len(self.constraints) :
      print 'Warning: top and bottom constraints have unequal length' 
      raise SystemExit ('this will conflict with later secondary structure '+
          'constraints in the design process')

  def __str__(self):
    if self.tag == 'toehold':
      return "t" + str(self.id)
    else:
      return "d" + str(self.id)

  def __invert__(self):
    return NusCDomain(self)

  def __eq__(self, other):
    if not isinstance(other, NusDomain): return False
    return self.id == other.id
  def __ne__(self, other):
    return not (self == other)

class NusCDomain(dnaobjects.ComplementaryDomain):
  """ Inherits from DNAObjects

  It makes use of __init__() ...
  
  """
  #def __init__(self, bottom=None, **kwargs):
  #  super(NusCDomain, self).__init__(**kwargs)

  def constraint(self):
    return self._complement.bottom

  def __invert__(self):
    return self._complement

  ## (In)equality
  def __eq__(self, other):
    """ Returns True iff their complements are equal."""
    if not isinstance(other, NusCDomain): return False
    return self._complement.__eq__(other.complement)
  def __ne__(self, other):
    """ Returns True iff they are not equal."""
    return not self.__eq__(other)

  def __str__( self ):
    if self._complement.tag == 'toehold':
      return "t" + str(self._complement.id) + "*"
    else:
      return "d" + str(self._complement.id) + "*"

class Structure(object):
  """The Structure of a DNA complex.

  The structure is specified using a list of domains and a corresponding
  dot-parens notation. 

  Args:
    domains (List[str]): names of domains
    dotparens (List[str]): A list of dot-parens characters ".(.+)."
    attr (Dict[name]=Domain): a mapping between names and Domain objects
  """
  def __init__(self, domains, dotparens, attr):
    self.domains = domains
    self.dotparens = dotparens
    self.attributes = attr

class Solution(object):
  """A set of structure objects.

  The solution object contains all species that have to be present in
  ``infinite'' concentration in a DSD circuit.

  Args:
    initial (set[Solution]) : A set of Solution objects
  """
  def __init__(self, initial):
    assert isinstance(initial, set)
    all(isinstance(s, Structure) for s in initial)
    self.molecules = initial

  def __add__(self, other):
    """Returns the union of two solution objects."""
    assert isinstance(other, Solution)
    return Solution(self.molecules.union(other.molecules))
  

class builtin_expressions(object):
  """Supported expressions of the nuskell programming language.

  *) Environment dependent functions (prefix ed_) interpret their arguments in
  the context of the current environment, and they *often* actually update the
  environment with supplied content.

  *) trailer functions (prefix tf_) have an additional head argument that is
  ...
  """
  #def __init__(self, env):
  #  # Shall we initialize this every time with a new environment? -- NO
  #  self.env = env

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
        'dna': self._dna, # return self._env, Structure()
        'list': self._list, # return self._env, content
        'where': self._where, # return self._env, value
        'uminus': self._uminus, # return self._env, -value (integer only!)
        #'trailer': self._trailer, # return self._env, head
        }

    tag = expr[0]
    content = expr[1:]

    #print 't', tag, content
    if tag in operators.keys():
      self._env, operand1 = self.interpret_expr(content[0])
      self._env, operand2 = self.interpret_expr(content[1])
      return self._env, operators[tag](operand1, operand2)
    elif tag == 'trailer' :
      """ nested stuff ... 
      has two list arguments, the function name (x[0], args=x[1:]
      function call, list indexing, or accessing attribute of an object
      """
      trailerkeys = {
          'apply' : self._apply,
          'index' : self._index,
          'attribute' : self._attribute
          }

      #print tag, content

      # function call, list attribute of an object
      self._env, head = self.interpret_expr(content[0])

      #print tag, head

      for x in content[1:]:
        key = x[0]
        args = x[1:]
        #print 'k', key, 'a', args
        self._env, head = trailerkeys[key](self, head, args)

      return self._env, head
    elif tag in keywords.keys():
      return keywords[tag](self, content)
    else :
      raise RuntimeError("Unknown expression `" + tag + "' was found.")
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
    #print content
    domains = content[0]
    dotparen = content[1]
    attributes = {}
    for i in range(len(domains)):
      if domains[i] != "?" and domains[i] != "+":
        starred = (len(domains[i]) == 2)
        dom = domains[i][0]
        #print 'early1', dom, starred
        theenv._env, dom_value = theenv.interpret_expr(dom)        
        #print 'early', dom, ':', dom_value
        attributes[dom[1]] = dom_value
        #print dom[1], ':', dom_value
        if starred:
          dom_value = ~dom_value
        domains[i] = dom_value
        #print 'what', domains, ':', dom_value
    #print "r", domains, dotparen, attributes
    return theenv._env, Structure(domains, dotparen, attributes)

  def _list(self, theenv, content) :
    #print 'l', content
    for i in range(len(content)):
      theenv._env, content[i] = theenv.interpret_expr(content[i])
    #print 'l', content
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
        id_list = remove_id_tags(id_list)

        for key, value in asgn_pattern_match(id_list, value) :
          theenv._env = theenv._create_binding(key, value)
      # evaluate the final value

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
    #print "apply", head, args
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
    if isinstance(head, Structure):
      head = head.attributes[identifier]
    elif identifier in head.__dict__:
      head = head.__dict__[identifier]
    else:
      raise RuntimeError(
          "The attribute `"+identifier+"' could not be found.")
    return theenv._env, head

class builtin_functions(object):
  """Builtin functions of the nuskell programming language.

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

    if 'len' not in kwargs and 'top' not in kwargs:
      kwargs['len'] = self._long_domain_length

    if 'len' in kwargs :
      # length is a short-cut for top + bottom in standard 3-letter alphabet
      for x in ['constraints', 'top', 'bottom']: 
        assert x not in kwargs
      kwargs['constraints'] = 'H' * kwargs['len']
      kwargs['bottom'] = 'D' * kwargs['len']
      del kwargs['len']
    elif 'top' in kwargs :
      assert 'constraints' not in kwargs
      assert 'bottom' in kwargs
      kwargs['constraints'] = kwargs['top']
      del kwargs['top']

    if 'tag' in kwargs:
      tag = kwargs['tag']
    else :
      tag = 'branch-migration'
    return NusDomain(domaintag=tag, **kwargs)

  def short(self, args):
    """A function that returns toehold domains.
 
    # supported keywords: 'len', 'top', 'bottom', 'tag'
    # len is short for standard 'top & bottom'
    # tag supports toehold and branch-migration which is
    # currently equal to short() and long()
   
    """
    if args :
      kwargs = args[0]
    else :
      kwargs = {'len' : self._short_domain_length}

    if 'len' not in kwargs and 'top' not in kwargs:
      kwargs['len'] = self._short_domain_length

    if 'len' in kwargs :
      # length is a short-cut for top + bottom in standard 3-letter alphabet
      for x in ['constraints', 'top', 'bottom']: 
        assert x not in kwargs
      kwargs['constraints'] = 'H' * kwargs['len']
      kwargs['bottom'] = 'D' * kwargs['len']
      del kwargs['len']
    elif 'top' in kwargs :
      assert 'constraints' not in kwargs
      assert 'bottom' in kwargs
      kwargs['constraints'] = kwargs['top']
      del kwargs['top']

    if 'tag' in kwargs:
      tag = kwargs['tag']
    else :
      tag = 'toehold'

    return NusDomain(domaintag=tag, **kwargs)

  def _print(self, args):
    """Print statment, primarily to debug nuskell scripts"""
    for a in args:
      if isinstance(a, Structure) :
        print [str(d) for d in a.domains]
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
    'infty' is a built-in function that takes a structure instance and puts
    that the species is supposed to have "infinite" concentration. The
    resulting object will be a Solution instance
    """
    if not isinstance(args[0], Structure):
      raise RuntimeError("The argument of `infty' should be a structure")
    return Solution(set([args[0]]))

  #@staticmethod
  #def unique(args) :
  #  if type(args[0]) != int:
  #    raise RuntimeError("The first argument of `unique' should be an integer.")
  #  return Domain(args[0])

  @staticmethod
  def complement(args) :
    """ Returns the complement of a given input. It can interpret Structures, 
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
    elif isinstance(x, NusDomain):
      print 'Untested:', ~x
      return ~x
    elif isinstance(x, Structure):
      return Structure(
          self._complement([x.domains]),
          self._complement([x.dotparens]),
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
      reversible = False
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

class Environment(builtin_expressions):
  """The environment of the tls language. It collects built-in functions
  together with the functions specified in the translation scheme, in order to
  compile a DNA strand displacement circuit.
  """

  def __init__(self, name, sdlen = 6, ldlen = 15):
    """
    Initializes the environment. Takes a name which usually refers to the
    current translation scheme, and parses a sample code that comprises
    built-in utilities. It seems that the length of domains is not so important
    at this point, ... It influences the sequence length in the .pil file, 
    but choosing an appropriate lenght might be more important for sequence
    design itself, rather than for domain level design. 
    
    :param name: Name of the translation scheme
    :param sdlen: Length of inbuilt **short()** domains
    :param ldlen: Length of inbuilt **long()** domains

    """
    self._name = name
    self._env = [{}]

    # Setup the builtin functions. 
    self._env, self._bfunc_Obj = self._init_builtin_functions(sdlen, ldlen)
    #print self._bfunc_Obj._long_domain_length 

  def _init_builtin_functions(self, sdlen, ldlen):
    """Initialization of builtin functions of the nuskell environment

    **TODO** this is where we want to discribe what: tail, complement, infty,
    short, long, unique, flip, empty, rev_reactions, irrev_reactions, ... do!
    """
    self._create_binding("print", Function(["s"], "print"))
    self._create_binding("abort", Function(["s"], "abort"))
    self._create_binding("tail", Function(["l"], "tail"))
    self._create_binding("flip", Function(["l", "n"], "flip"))
    self._create_binding("long", Function([], "long"))
    self._create_binding("short", Function([], "short"))
    self._create_binding("infty", Function(["species"], "infty"))
    #self._create_binding("unique", Function(["l"], "unique"))
    self._create_binding("complement", Function(["l"], "complement"))
    self._create_binding("rev_reactions", Function(["crn"], "rev_reactions"))
    self._create_binding("irrev_reactions", Function(["crn"], "irrev_reactions"))
    self._create_binding("empty", Solution(set()))
    return self._env, builtin_functions(sdlen, ldlen)

  def _create_binding(self, name, value):
    """Create a binding for a new nuskell function in the top level

    :param name: Function name
    :param value: Some built-in data type to be adressed (Function, Solution,
    Species, Domain, Reaction, Structure), either as single value or list

    :return: updated environment
    """

    bindings = (Function, Solution, Species, NusDomain, Reaction, Structure, 
        void, int, list)
    #print 'n:', name, 'v:', type(value), value
    if isinstance(value, list) :
      assert all(isinstance(s, bindings) for s in value)
    else :
      assert isinstance(value, bindings)

    #print "create binding:", len(self._env), name, value

    self._env[-1][name] = value
    return self._env

  # Private environment modification functions #
  def _create_level(self): # _create_env_level
    """Private: Create a new level for function bindings.
    It is commonly triggered by a 'where' statement or during the evaluation of
    a non-built-in function call.

    interpret_expr -> _where
    interpret_expr -> _trailer -> _apply -> _eval_func
    """
    self._env.append({})
    return self._env
  
  def _destroy_level(self): # _destroy_env_level_
    """Private: Revert to the old level for function bindings.
    It is commonly triggered after a 'where' statement or after the evaluation
    of a non-built-in function call.

    interpret_expr -> _where
    interpret_expr -> _trailer -> _apply -> _eval_func
    """
    self._env.pop()
    return self._env

  def _ref_binding(self, name):
    """Get the reference to an exisiting function binding, by searching all
    levels.

    :param name: Function name

    :return: Function binding (self._env[?][name])
    """
    for level in reversed(self._env):
      if name in level.keys():
        return level[name]
    raise RuntimeError("Cannot find a binding for `" + name + "'.")
    return None

  ### documentation ends here ... ###
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

  def interpret(self, code):
    """Creates bindings for variables and functions defined in the body of the
    code. Returns the environment (the final namespace).

    Note: 
      This function creates bindings for every function in the nuskell code.
      At this point all keywords (class, function, macro, module) are treated
      exactly the same, only the **global** keyword is special, because the 
      global expressions are **interpreted first** and then bound.
    """
    for stmt in code:
      kwd = stmt[0]
      body = stmt[1:]

      if kwd == "global":
        # create binding to interpretation of the expression
        id_list = body[0]
        id_list = remove_id_tags(id_list)

        value = body[1]
        self._env, value = self.interpret_expr(value)

        for key, value in asgn_pattern_match(id_list, value) :
          self._create_binding(key, value)
      else:
        # create binding to the function, without interpretation
        # e.g. kwd = module; id = 'rxn'; args = ['r']; body_ = [where...]
        assert body[0][0] == "id"
        #print "k:", kwd, "b:", body
        id = body[0][1]
        args = map(lambda x: x[1], body[1]) # remove 'id' tags from args
        body_ = body[2] # the ['where' [...]] part
        #print "id:", id, 'a', args, 'b2', body_
        self._create_binding(id, Function(args, body_))
    return self._env

  def translate_formal_species(self, fs_list):
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

    return self._env, self.formal_species_dict

  def translate_reactions(self, crn_parsed):
    """Execute the main() function of the translation scheme.
    
    The crn is combined with the compiled formal species.

    Args:
      crn_paresed (List[Lists]): A crn in crn_parser format.

    Returns:
      self.env : The updated environment
      self.main_solution (Solution) : The final Solution object

    Raises:
      RuntimeError: If the compiled formal species cannot be found

    """
    if not self.formal_species_dict :
      raise RuntimeError('Could not find the compiled formal species!')
    crn_remap = map(
        lambda x: [x[0]] + map(
          lambda y: map(
            lambda z: self.formal_species_dict[z], y), x[1:]), crn_parsed)
    crn_object = map(
        lambda x: Reaction(x[1], x[2], x[0] == "reversible"), crn_remap)

    # main(__crn__)
    self._create_binding("__crn__", crn_object)
    self._env, self.constant_species_solution = self.interpret_expr( 
        ["trailer", ["id", "main"], ["apply", ["id", "__crn__"]]])

    return self._env, self.constant_species_solution

  # Public functions #
  @property
  def name(self) :
    return self._name

  @property
  def env(self) :
    return self._env

  # TODO: A method that prints the environment in some human readable way.
  # def __str__(): pass

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
      result += asgn_pattern_match(id[i], value[i])
    return result
  elif type(id) == str:
    return [(id, value)]
  else:
    raise RuntimeError("""The left hand side of an assignment should be 
    either an identifier or a list-tree consisting only of identifiers.""")

def remove_id_tags(l):
  """Helper function to remove all tags from the given id_list."""
  kwd = l[0]
  if kwd == "idlist":
    return map(remove_id_tags, l[1:])
  elif kwd == "id":
    return l[1]
  else:
    raise RuntimeError("""The left hand side of an assignment should be 
    either an identifier or a list-tree consisting only of identifiers.""")

