"""The tls language environment.

Recommended Input: 
  *) tls code parsed by the tls_parser module
  *) crns parsed by the crn_parser module

"""

import sys
from copy import copy

# exceptions
class void(object):
  """An empty object returned by the print function"""
  pass

class RuntimeError(Exception):
  def __init__(self, value):
    self.value = "Runtime error : " + value
  def __str__(self):
    return self.value

class Species(object):
  def __init__(self, name):
    self.name = name

class Reaction(object):
  """Format of a single reaction:

  :param reactands: a list of educts ['A', 'B']
  :param products: a list of products ['C']
  :param reversible: True or False
  """
  def __init__(self, reactants, products, reversible):
    self.reactants = reactants
    self.products = products
    self.reversible = reversible

class Function(object):
  """ A function of the nuskell language """
  def __init__(self, args, body):
    self.args = args
    self.body = body

class Solution(object):
  """A set of Structure Objects found *in Solution* """
  def __init__(self, initial):
    assert isinstance(initial, set)
    all(isinstance(s, Structure) for s in initial)
    self.molecules = initial

  def __add__(self, other):
    assert isinstance(other, Solution)
    return Solution(self.molecules.union(other.molecules))
  
class Domain(object):
  """ This is a weird class, it is possible to define domains with 
  the same id but different length, the comparison was messed up and
  the global domain_length dictionary is never used ...

  TODO: This needs a rewrite in my opinion, 
    *) remove domain_length 
    *) rewrite __init__ and __eq__
  """
  domain_id = 0
  #domain_length = {}

  def __init__(self, length, id = None):
    assert type(length) is int
    assert (type(id) is int) or (id is None)

    if id == None:
      Domain.domain_id += 1
      id = Domain.domain_id
      #Domain.domain_length[id] = length
    self.id = id
    self.length = length

  def __str__(self):
    if self.id < 0: # complement of something
      return "d" + str(-self.id) + "*"
    else:
      return "d" + str(self.id)

  def __invert__(self):
    return Domain(self.length, -self.id)

  def __eq__(self, other):
    if not isinstance(other, Domain): return False
    return self.id == other.id

  def __ne__(self, other):
    return not (self == other)

class Structure(object):
  def __init__(self, domains, dotparens, attr):
    self.domains = domains
    self.dotparens = dotparens
    self.attributes = attr

class _builtin(object):
  """A collection of built-in nuskell function types. Three flavours of
  functions are combined in this class. 
  
  *) Base-level functions (prefix bl_) are independent of the current
  environment, they take a list of arguments and return a processed version of
  these arguments.

  *) Environment dependent functions (prefix ed_) interpret their arguments in
  the context of the current environment, and they *often* actually update the
  environment with supplied content.

  *) trailer functions (prefix tf_) have an additional head argument that is
  ...
  """

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
    sys.exit('or has not been tested')
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

  def _dna(self, theenv, content):
    #print content
    domains = content[0]
    dotparen = content[1]
    attributes = {}
    for i in range(len(domains)):
      if domains[i] != "?" and domains[i] != "+":
        starred = (len(domains[i]) == 2)
        dom = domains[i][0]
        #print 'early1', dom
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

  # eval_builtin_functions
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
    raise SystemExit, 'EXIT: ' + args[0]

  def _tail(self, args) :
    """ Returns the last element.. """
    if type(args[0]) != list:
      raise RuntimeError("`tail' should have a list as its argument.")
    return args[0][1:]

  # # eval_builtin_functions
  # def _complain(self, args) :
  #   """ Write an Error message and quit() """
  #   if type(args[0]) != list:
  #     raise RuntimeError("`complain' with a string!")
  #   #print args[0]
  #   return args[1:]

  def _complement(self, args) :
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
    elif isinstance(x, Domain):
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
  
  def _infty(self, args) :
    """
    'infty' is a built-in function that takes a structure instance and puts
    that the species is supposed to have "infinite" concentration. The
    resulting object will be a Solution instance
    """
    if not isinstance(args[0], Structure):
      raise RuntimeError("The argument of `infty' should be a structure")
    return Solution(set([args[0]]))

  def _unique(self, args) :
    if type(args[0]) != int:
      raise RuntimeError("The first argument of `unique' should be an integer.")
    return Domain(args[0])

  def _flip(self, args) :
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

  def _rev_reactions(self, args) : 
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

  def _irrev_reactions(self, args) : 
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

class Environment(_builtin):
  """ The environment of the tls language. It collects built-in functions
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
    self._env = self.base_level_functions()
    self._short_domain_length = sdlen
    self._long_domain_length = ldlen

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

  def _create_binding(self, name, value):
    """Create a binding for a new nuskell function in the top level

    :param name: Function name
    :param value: Some built-in data type to be adressed (Function, Solution,
    Species, Domain, Reaction, Structure), either as single value or list

    :return: updated environment
    """

    bindings = (Function, Solution, Species, Domain, Reaction, Structure, void)
    #print 'n:', name, 'v:', type(value), value
    if isinstance(value, list) :
      assert all(isinstance(s, bindings) for s in value)
    else :
      assert isinstance(value, bindings)

    #print "create binding:", len(self._env), name, value

    self._env[-1][name] = value
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
      return self._env, self._eval_builtin_func(f.body, args)
  
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

  def _eval_builtin_func(self, f, args):
    """Evaluate built-in functions. These functions do not alter the
    environment and they are contained in their own class for readability.
    """

    keywords = {
        'print': self._print,
        'abort': self._abort,
        'tail' : self._tail,
        'flip' : self._flip,
        'long' : lambda x: Domain(self._long_domain_length),
        'short' : lambda x: Domain(self._short_domain_length),
        'infty' : self._infty,
        'unique' : self._unique,
        #'complain' : self._complain,
        'complement' : self._complement,
        'rev_reactions' : self._rev_reactions,
        'irrev_reactions' : self._irrev_reactions,
        }

    if f in keywords.keys() :
      return keywords[f](args)
    else :
      raise RuntimeError("`" + f + "' could not be resolved.")
      return None

  # Public functions #
  @property
  def name(self) :
    return self._name

  @property
  def env(self) :
    return self._env

  def print_environment(self) :
    """ Print a snapshot of the current environment """
    for dic in self._env:
      for func in dic :
        if isinstance(dic[func], Function):
          print "f", func, dic[func].args, dic[func].body
        elif isinstance(dic[func], Solution):
          print "s", func, dic[func].molecules
        else :
          print "unknown thingy"

  def base_level_functions(self):
    """ Initialization of builtin functions of the nuskell environment

    **TODO** this is where we want to discribe what: tail, complement, infty,
    short, long, unique, flip, empty, rev_reactions, irrev_reactions, ... do!
    """
    self._create_binding("print", Function(["s"], "print"))
    self._create_binding("abort", Function(["s"], "abort"))
    self._create_binding("tail", Function(["l"], "tail"))
    self._create_binding("complement", Function(["l"], "complement"))
    self._create_binding("infty", Function(["species"], "infty"))
    self._create_binding("short", Function([], "short"))
    self._create_binding("long", Function([], "long"))
    self._create_binding("unique", Function(["l"], "unique"))
    self._create_binding("flip", Function(["l", "n"], "flip"))
    self._create_binding("empty", Solution(set()))
    self._create_binding("rev_reactions", Function(["crn"], "rev_reactions"))
    self._create_binding("irrev_reactions", Function(["crn"], "irrev_reactions"))
    return self._env

  def interpret(self, code):
    """ 
    Creates bindings for variables and functions defined in the body of the
    code. Returns the environment (the final namespace).
    """
    for stmt in code:
      kwd = stmt[0]
      body = stmt[1:]
      if kwd == "global":
        #print "GLOBAL", kwd, 'b:', body
        id_list = body[0]
        id_list = remove_id_tags(id_list)

        value = body[1]
        #print 'id', id_list, 'v', value
        self._env, value = self.interpret_expr(value)

        for key, value in asgn_pattern_match(id_list, value) :
          self._create_binding(key, value)
      else:
        # function declaration
        assert body[0][0] == "id"
        #print "k:", kwd, "b:", body
        id = body[0][1]
        args = map(lambda x: x[1], body[1])
        body_ = body[2]
        #print id, body_, args
        self._create_binding(id, Function(args, body_))
    return self._env


  def interpret_expr(self, v):
    """ Recursive interpretation the body of a global variable """
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
        'dna': self._dna, # return self._env, Structure(domains, dotparen, attributes)
        'list': self._list, # return self._env, content
        'where': self._where, # return self._env, value
        'uminus': self._uminus, # return self._env, -value (integer only!)
        #'trailer': self._trailer, # return self._env, head
        }

    tag = v[0]
    content = v[1:]

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


      for x in content[1:]:
        key = x[0]
        args = x[1:]
        self._env, head = trailerkeys[key](self, head, args)

      return self._env, head
    elif tag in keywords.keys():
      return keywords[tag](self, content)
    else :
      raise RuntimeError("Unknown expression `" + tag + "' was found.")
      return self._env, None


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
  """Helper function to remove all tags from the given id_list. """
  kwd = l[0]
  if kwd == "idlist":
    return map(remove_id_tags, l[1:])
  elif kwd == "id":
    return l[1]
  else:
    raise RuntimeError("""The left hand side of an assignment should be 
    either an identifier or a list-tree consisting only of identifiers.""")

