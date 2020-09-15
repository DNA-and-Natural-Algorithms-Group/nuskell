#
#  nuskell/dsdcompiler/interpreter.py
#  NuskellCompilerProject
#
# Copyright (c) 2010-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
"""The nuskell programming language interpreter.

Nuskell code is interpreted using public functions of the **Environment**
class.  All other classes and functions are needed internally to set up the
interpreter environment.  Built-in functions and expressions are encapuslated
in the respective namespace.

.. When modifying this file, follow the Google Python Style Guide:
.. http://google.github.io/styleguide/pyguide.html
"""
from copy import copy

from .objects import NuskellDomain, NuskellComplex, DSD_Complex, DSDDuplicationError

class NuskellEnvError(Exception):
    """Nuskell Environment Error.

    Args:
      msg (str): Error description.
    """

    def __init__(self, msg):
        self.message = msg
        super(NuskellEnvError, self).__init__(self.message)

class NuskellExit(SystemExit):
    """Nuskell Scheme Exit """

    def __init__(self, msg):
        super(NuskellExit, self).__init__(msg)


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


class NusFunction(object):
    """
    A function object of the Nuskell language.

    Args:
      args (): The arguments of a function.
      body (): The function body.

    """

    def __init__(self, args, body):
        self.args = args
        self.body = body

def flatten(l):
    if l == []:
        return l
    if isinstance(l[0], list):
        return flatten(l[0]) + flatten(l[1:])
    return l[:1] + flatten(l[1:])

class NusComplex(DSD_Complex):
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


    def __new__(cls, attr=dict(), *kargs, **kwargs):
        # The new method returns the present instance of an object, if it exists
        self = DSD_Complex.__new__(cls)
        try:
            super(NusComplex, self).__init__(*kargs, **kwargs)
            self.attributes = attr
        except DSDDuplicationError as e:
            return e.existing
        return self

    def __init__(self, attr=dict(), *kargs, **kwargs):
        # Remove default initialziation to get __new__ to work
        pass

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
            if isinstance(s, NuskellDomain):
                assert (c in ['(', '.', ')'])
                return s, c
            elif s == '+' and c == '+':
                return s, c
            elif isinstance(s, NusComplex):
                assert (c == '~')
                return s.sequence, s.structure
            elif isinstance(s, list):
                assert (c == '~')
                # map returns a list of [[seq, str], ..]
                # zip transposes this into tuples
                # ([seq1,seq2,...]),([str1,str2,...])
                return [list(a) for a in zip(*map(lambda i: resolve(s[i], '~'),
                                                  range(len(s))))]
            elif s == '?':
                assert (c in ['(', '.', ')'])
                return NuskellDomain(prefix='h', dtype='long'), c
                # return s, c
            else:
                raise NotImplementedError

        newseq = []
        newstr = []
        for i in range(len(self._sequence)):
            if self._sequence[i] == [] and self._structure[i] == '~':
                continue
            se, ss = resolve(self._sequence[i], self._structure[i])
            newseq.append(se)
            newstr.append(ss)

        newseq = flatten(newseq)
        newstr = flatten(newstr)

        start = True
        for i in reversed(range(len(newseq))):
            if newseq[i] == '+' and newstr[i] == '+':
                if start or i == 0:
                    newseq.pop(i)
                    newstr.pop(i)
                start = True
            else:
                start = False

        self._sequence = newseq
        self._structure = newstr
        return self._sequence, self._structure


class NuskellExpressions(object):
    """
    Builtin expressions for Nuskell translation schemes.
    """
    # Say something about special trailer functions here?

    def interpret_expr(self, expr):
        """Recursive interpretation the body of a global variable."""
        operators = {"*": lambda x, y: x * y,
                     "/": lambda x, y: x / y,
                     "+": lambda x, y: x + y,
                     "-": lambda x, y: x - y,
                     "==": lambda x, y: x == y,
                     "!=": lambda x, y: x != y,
                     ">": lambda x, y: x > y,
                     "<": lambda x, y: x < y,
                     ">=": lambda x, y: x >= y,
                     "<=": lambda x, y: x <= y}

        keywords = {
            'id': lambda self, content: (self._env, self._ref_binding(content[0])),
            'if': self._if,  # return self._env, value
            'or': self._or,  # return self._env, operand1 or operand2
            'and': self._and,  # return self._env, operand1 and operand2
            'num': lambda self, content: (self._env, int(content[0])),
            'quote': lambda self, content: (self._env, content[0]),
            'dict': self._dict,
            'dna': self._dna,  # return self._env, NusComplex()
            'list': self._list,  # return self._env, content
            'where': self._where,  # return self._env, value
            'uminus': self._uminus,  # return self._env, -value (integer only!)
        }

        tag = expr[0]
        content = expr[1:]

        if tag in operators.keys():
            self._env, operand1 = self.interpret_expr(content[0])
            self._env, operand2 = self.interpret_expr(content[1])
            return self._env, operators[tag](operand1, operand2)

        elif tag == 'trailer':
            trailerkeys = {
                'apply': self._apply,
                'index': self._index,
                'attribute': self._attribute
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

        else:
            raise NuskellEnvError("Unknown expression: `" + tag + "'")
            return self._env, None

    # Context-dependent
    @staticmethod
    def _if(theenv, content):
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
            if asign[1][0] == 'num':
                asign[1][1] = int(asign[1][1])
            if asign[1][0] == 'quote':
                asign[1][1] = asign[1][1][1:-1]
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
        return theenv._env, NusComplex(sequence=domains, structure=dotparen,
                                       attr=attributes, memorycheck=False)

    def _list(self, theenv, content):
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

                for key, value in self._bfunc_Obj.asgn_pattern_match(
                        id_list, value):
                    theenv._create_binding(key, value)

        theenv._env, value = theenv.interpret_expr(content[0])
        theenv._destroy_level()
        return theenv._env, value

    def _uminus(self, theenv, content):
        theenv._env, value = theenv.interpret_expr(content[0])
        if not isinstance(value, int):
            raise NuskellEnvError(
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
        if not isinstance(head, list):
            raise NuskellEnvError("Only lists can be indexed.")
        if not isinstance(subscript, int):
            raise NuskellEnvError("Subscript should be an integer.")
        try:
            head = head[subscript]
        except IndexError:
            raise NuskellEnvError(
                "Bug in translation scheme, expected element in empty list.")
        return theenv._env, head

    def _attribute(self, theenv, head, args):
        identifier = args[0][1]  # strip the tag
        if isinstance(head, NusComplex):
            head = head.attributes[identifier]
        elif identifier in head.__dict__:
            head = head.__dict__[identifier]
        else:
            raise NuskellEnvError(
                "The attribute `" + identifier + "' could not be found.")
        return theenv._env, head


class NuskellFunctions(object):
    """Builtin functions of Nuskell translation schemes.

    This object is initialized within the :obj:`NuskellEnvironment()` after the
    respective functions have been bound.  Most methods do not require any class
    variables and are therefore declared as staticmethods.  As a consequence, the
    functions can be accessed without prior initialization of the Object.

    """
    ################################
    ### Builtin-function section ###
    ################################

    def __init__(self, *kargs, **kwargs):
        super(NuskellFunctions, self).__init__(*kargs, **kwargs)

    def eval_builtin_functions(self, f, args):
        """Evaluate built-in functions.

        These functions are independent of the environment.

        Args:
          f (str) : The built-in function name.
          args (list): The arguments for the function call.

        Raises:
          NuskellEnvError: "Function could not be found."

        Returns:
          The results form f(args)
        """

        keywords = {
            'print': self._print,
            'abort': self._abort,
            'tail': self.tail,
            'flip': self.flip,
            'long': self.long,
            'short': self.short,
            'infty': self.infty,
            'unique' : self.unique,
            'complement': self.complement,
            'rev_reactions': self.rev_reactions,
            'irrev_reactions': self.irrev_reactions,
        }

        if f in keywords.keys():
            return keywords[f](args)
        else:
            raise NuskellEnvError("`" + f + "' could not be found.")

    def unique(self, args):
        """A function that returns branch-migration domains.
        """
        if type(args[0]) != int:
            raise RuntimeError("The first argument of `unique' should be an integer.")

        dom = NuskellDomain(length=args[0], prefix='u')
        cdm = ~dom
        return dom


    def long(self, args):
        """A function that returns branch-migration domains.
        """
        if args:
            raise args

        dom = NuskellDomain(dtype='long', prefix='d')
        cdm = ~dom
        return dom

    def short(self, args):
        """A function that returns toehold domains.
        """
        dom = NuskellDomain(dtype='short', prefix='t')
        cdm = ~dom
        return dom

    def _print(self, args):
        """Print statment, primarily to debug nuskell scripts"""
        for a in args:
            if isinstance(a, NusComplex):
                print(list(map(str, a.sequence)))
            elif isinstance(a, NuskellComplex):
                print(list(map(str, a.sequence)))
            else:
                print(str(a), end='')
        print()
        return void()

    def _abort(self, args):
        """Raise NuskellExit and print the error message"""
        raise NuskellExit(args[0])

    @staticmethod
    def tail(args):
        """ Returns the list minus the first element.. """
        if not isinstance(args[0], list):
            raise NuskellEnvError("`tail' should have a list as its argument.")
        return args[0][1:]

    @staticmethod
    def flip(args):
        """A matrix transpose function for lists of lists.

        Args:
          args (list): Contains both the lol and and integer to specify the
            dimenstion.

        Example:
          input: [[[x,y,z],[a,b,c],[n,l,k],[u,v,w]],3]
          returns: [[x,a,n,u],[y,b,l,v],[z,c,k,w]]
        """
        if not isinstance(args[0], list):
            raise NuskellEnvError(
                "The first argument of `flip' should be a list.")
        if not isinstance(args[1], int):
            raise NuskellEnvError(
                "The second argument of `flip' should be an integer.")

        l = args[0]
        n = args[1]
        res = [0] * n
        for i in range(n):
            res[i] = [0] * len(l)
        for i in range(len(l)):
            if len(l[i]) != n:
                raise NuskellEnvError(
                    "The argument of `flip' does not have a required form.")
            for j in range(n):
                res[j][i] = l[i][j]
        return res

    @staticmethod
    def infty(args):
        """
        'infty' is a built-in function that takes a Complex instance and puts
        that the species is supposed to have "infinite" concentration. The
        resulting object will be a dictionary
        """
        if not isinstance(args[0], NusComplex):
            raise NuskellEnvError("The argument of `infty' should be a Complex")
        args[0].flatten_cplx
        if args[0].sequence == [] and args[0].structure == []:
            return []
        else:
            try : 
                final = NuskellComplex(args[0].sequence, 
                                       args[0].structure, 
                                       prefix = 'final')
                final.concentration = ('constant', float('inf'), 'M')
                return [final]
            except DSDDuplicationError as e:
                final = e.existing
                return []

    @staticmethod
    def complement(args):
        """ Returns the complement of a given input. It can interpret Complex,
        Domains and Dot-bracket lists, as well as lists of lists.

        :param args: an expression in form of an array, where the first element
        in the array contains the data of interest. ... yeah, a bit complicated...

        :return: complement of the input
        """
        x = args[0]
        if isinstance(x, list):
            # args[0] forces us to introduce additional lists ...
            for i in range(len(x)):
                x[i] = [x[i]]
            return list(reversed(map(self._complement, x)))
        elif isinstance(x, Domain):
            print('Untested:', ~x)
            return ~x
        elif isinstance(x, NusComplex):
            raise NotImplementedError
        elif x == "(":
            return ")"
        elif x == ")":
            return "("
        else:
            return x

    @staticmethod
    def rev_reactions(args):
        if not isinstance(args[0], list):
            raise NuskellEnvError(
                "The argument of `rev_reactions' should be a list.")
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
    def irrev_reactions(args):
        """
        a function that divides every reversible reaction into a pair of
        irreversible reaction.
        """
        if not isinstance(args[0], list):
            raise NuskellEnvError(
                "The argument of `irrev_reactions' should be a list.")
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
        if isinstance(id, list):
            if not isinstance(value, list) or len(id) != len(value):
                raise NuskellEnvError(
                    "Pattern matching failed while assigning values to " + str(id) + ".")
            result = []
            for i in range(len(id)):
                result += NuskellFunctions.asgn_pattern_match(id[i], value[i])
            return result
        elif isinstance(id, str):
            return [(id, value)]
        else:
            raise NuskellEnvError("""The left hand side of an assignment should be
      either an identifier or a list-tree consisting only of identifiers.""")

    @staticmethod
    def remove_id_tags(l):
        """Helper function to remove all tags from the given id_list."""
        kwd = l[0]
        if kwd == "idlist":
            return list(map(NuskellFunctions.remove_id_tags, l[1:]))
        elif kwd == "id":
            return l[1]
        else:
            raise NuskellEnvError("""The left hand side of an assignment should be
      either an identifier or a list-tree consisting only of identifiers.""")


class NuskellEnvironment(NuskellExpressions):
    """The Nuskell language environment.

    :obj:`NuskellEnvironment()` interprets the Nuskell language using low-level
    instructions, and executes the formal and main functions of translation
    schemes. Inherits from :obj:`NuskellExpressions()` and initializes built-in
    functions provided by :obj:`NuskellFunctions()`.

    Raises:
      NuskellEnvError

    """

    def __init__(self):
        # Setup the builtin functions.
        self._env = [{}]
        self._bfunc_Obj = self._init_builtin_functions()

    # Public functions #
    def interpret(self, code):
        """Build the Environment (the final namespace).

        Creates bindings for variables and functions defined in the body of the
        translation scheme.

        Note:
          All keywords (class, function, macro, module) are treated the same, only
          the **global** keyword is special, as global expressions are
          **interpreted first** and then bound.
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

                for key, value in self._bfunc_Obj.asgn_pattern_match(
                        id_list, value):
                    self._create_binding(key, value)
            else:
                # create binding to the function, without interpretation
                # e.g. kwd = module; id = 'rxn'; args = ['r']; body_ =
                # [where...]
                assert body[0][0] == "id"
                id = body[0][1]
                # remove 'id' tags from args
                args = list(map(lambda x: x[1], body[1]))
                body_ = body[2]  # the ['where' [...]] part
                self._create_binding(id, NusFunction(args, body_))

    def translate_formal_species(self, fs_list):
        """ Apply the formal() function to the formal species in the input CRN.

        First, the bindings for formal species are created, then the formal()
        function is applied to every formal species in the input CRN.

        Returns:
          [dict()] A dictionary of key=name, value=:obj:`NusComplex()`
        """
        formal_species_objects = list(map(Species, fs_list))

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

    def translate_reactions(self, crn_parsed, modular=False):
        """Execute the main() function of the translation scheme.

        The input CRN is replaced with previously initialized formal species objects.

        Args:
          crn_paresed (List[Lists]): A crn in crn_parser format.

        Raises:
          NuskellEnvError: If the compiled formal species cannot be found

        Returns:
            dict: fuel species 
        """
        if not self.formal_species_dict:
            raise NuskellEnvError(
                'Could not find the compiled formal species!')

        # replace every fs (str) with fs(NusComplex())
        crn_remap = list(map(lambda rxn: [rxn.k_rev] + 
                list(map(lambda y: list(map(lambda z: self.formal_species_dict[z], y)), rxn[:2])), crn_parsed))

        crn_object = list(map(
            lambda x: Reaction(x[1], x[2], x[0] != 0), crn_remap))

        # main(__crn__)
        modules = []
        if modular:
            for module in crn_object:
                self._create_binding("__crn__", [module])
                self._env, self.constant_species_solution = self.interpret_expr(
                    ["trailer", ["id", "main"], ["apply", ["id", "__crn__"]]])
                modules.append(self.constant_species_solution)

            # NOTE: The modular reaction-by-reaction translation does not *sum* over
            # the whole CRN, but sums over every reaction. The common memory management
            # of NusComplexes ensures that species wich are shared between reactions 
            # have the same name.
            unified_solution = sum(modules)
            modules.insert(0, unified_solution)

        else:
            self._create_binding("__crn__", crn_object)
            self._env, self.constant_species_solution = self.interpret_expr(
                ["trailer", ["id", "main"], ["apply", ["id", "__crn__"]]])
            modules = [self.constant_species_solution]

        return modules

    # Private Functions #
    def _init_builtin_functions(self):
        """ Initialized built-in Nuskell functions.

        Add a binding for every builtin-fuction to the environment.

        Returns:
          [:obj:`NuskellFunctions()`] Execution object for built-in Nuskell
          functions.
        """
        self._create_binding("print", NusFunction(["s"], "print"))
        self._create_binding("abort", NusFunction(["s"], "abort"))
        self._create_binding("tail", NusFunction(["l"], "tail"))
        self._create_binding("flip", NusFunction(["l", "n"], "flip"))
        self._create_binding("long", NusFunction([], "long"))
        self._create_binding("short", NusFunction([], "short"))
        self._create_binding("infty", NusFunction(["species"], "infty"))
        self._create_binding("unique", NusFunction(["l"], "unique"))
        #self._create_binding("complement", NusFunction(["l"], "complement"))
        self._create_binding(
            "rev_reactions", NusFunction(
                ["crn"], "rev_reactions"))
        self._create_binding(
            "irrev_reactions", NusFunction(
                ["crn"], "irrev_reactions"))
        self._create_binding("empty", [])
        return NuskellFunctions()

    def _create_binding(self, name, value):
        """ Create binding of a Nuskell function.

        Adds a binding to the last level of function bindings.

        Args:
          name (str): Name of the function.
          value (...): A built-in data type (Function, dict, Species,
            Domain, Reaction, Complex, etc ...), either as single value or list

        Returns:
          An updated environment inclding the function binding in the top-level.
        """

        bindings = (NusFunction, dict, Species, NuskellDomain, Reaction, NusComplex,
                    void, int, list)
        if isinstance(value, list):
            assert all(isinstance(s, bindings) for s in value)
        else:
            assert isinstance(value, bindings)

        self._env[-1][name] = value

    # Private environment modification functions #
    def _create_level(self):
        """
        Create a new level for function bindings.

        Note:
          This commonly initializes a level for a 'where' statement or the
          evaluation of a not built-in function call.
        """
        # interpret_expr -> _where
        # interpret_expr -> _trailer -> _apply -> _eval_func
        self._env.append({})

    def _destroy_level(self):
        """
        Revert to previous level of function bindings.

        Note:
          This commonly closes a 'where' environment or the evaluation of a
          not built-in function call.
        """
        # interpret_expr -> _where
        # interpret_expr -> _trailer -> _apply -> _eval_func
        self._env.pop()

    def _ref_binding(self, name):
        """
        Search levels of function bindings from last to first to find a reference.

        Args:
          name (str) : Name of a function

        Returns:
          Function binding (self._env[?][name])
        """
        for level in reversed(self._env):
            if name in level.keys():
                return level[name]
        raise NuskellEnvError("Cannot find a binding for `" + name + "'.")

    def _eval_func(self, f, args):
        """
        Evaluate a function.

        The function can either be one of the built-in functions or a module,
        class, macro, etc. as specified in the translation scheme.

        Args:
          f (:obj:`NusFunction()`): The function to evaluate
          args (list): A list of arguments for the function.
        """
        if not isinstance(f, NusFunction):
            raise NuskellEnvError(str(f) + "is not a function.")

        if isinstance(f.body, str):  # the function is a built-in function
            return self._env, self._bfunc_Obj.eval_builtin_functions(
                f.body, args)

        def hardcopy_list(l):
            if not isinstance(l, list):
                return copy(l)
            return list(map(hardcopy_list, l))

        if len(f.args) > len(args):
            raise NuskellEnvError("The function `" + f.name + "' requires at least "
                                  + str(len(f.args)) + " arguments but only found " + str(args) +
                                  " arguments.")

        self._create_level()
        for i in range(len(f.args)):
            arg_name = f.args[i]
            arg_value = args[i]
            self._create_binding(arg_name, arg_value)

        self._env, value = self.interpret_expr(hardcopy_list(f.body))
        self._destroy_level()
        return self._env, value
