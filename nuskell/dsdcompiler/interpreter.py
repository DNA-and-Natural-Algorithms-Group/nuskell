#
#  nuskell/dsdcompiler/interpreter.py
#  NuskellCompilerProject
#
# Copyright (c) 2010-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
""" The nuskell programming language interpreter.

Nuskell code is interpreted using public functions of the **Environment**
class.  All other classes and functions are needed internally to set up the
interpreter environment.  Built-in functions and expressions are encapuslated
in the respective namespace.

.. When modifying this file, follow the Google Python Style Guide:
.. http://google.github.io/styleguide/pyguide.html
"""
import logging
log = logging.getLogger(__name__)

from copy import deepcopy
from .objects import NuskellDomain, NuskellComplex, SingletonError

class NuskellExit(SystemExit):
    pass

class NuskellEnvError(Exception):
    pass

class void:
    """ An empty object returned by the print function. """
    pass

class Species:
    """ A species of a CRN.

    Args:
      name (str): The name of a species.
    """
    def __init__(self, name):
        self.name = name

class Reaction:
    """ A reversible or irreversible reaction.

    Args:
      reactands (List[str]) : a list of educts, eg. ['A', 'B']
      products (List[str]) : a list of products, eg. ['C']
      reversible (bool) : True or False
    """
    def __init__(self, reactants, products, reversible):
        self.reactants = reactants
        self.products = products
        self.reversible = reversible

class NusFunction:
    """ A function of the nuskell programming language.

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

class ComplexFragment:
    """ A nucleic acid complex, i.e. a (sequence, structure) pair.

    This is an extension of dsdobjects.core.DSD_Complex(), which stores an
    additional dictionary of attributes. This dictionary maps domain names as
    specified in the translation-scheme to the unique Domain() object.

    Args:
      attr (Dict[name]=Domain): A mapping between names in the translation scheme
                                and Domain objects
      sequence (List[str]): Domain Objects
      structure (List[str]): A list of dot-parens characters ".(.+)."
    """
    def __init__(self, sequence, structure, attr = None):
        self.sequence = sequence
        self.structure = structure
        self.attributes = dict() if attr is None else attr

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
            elif isinstance(s, ComplexFragment):
                assert (c == '~')
                return s.sequence, s.structure
            elif isinstance(s, list):
                assert (c == '~')
                # map returns a list of [[seq, str], ..]
                # zip transposes this into tuples
                # ([seq1,seq2,...]),([str1,str2,...])
                return [list(a) for a in zip(*map(lambda i: resolve(s[i], '~'), range(len(s))))]
            elif s == '?':
                assert (c in ['(', '.', ')'])
                return NuskellDomain(prefix = 'h', dtype = 'long'), c
            else:
                raise NotImplementedError

        newseq = []
        newstr = []
        for i in range(len(self.sequence)):
            if self.sequence[i] == [] and self.structure[i] == '~':
                continue
            se, ss = resolve(self.sequence[i], self.structure[i])
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

        self.sequence = newseq
        self.structure = newstr
        return self.sequence, self.structure

class NuskellExpressions:
    """ Builtin expressions for Nuskell translation schemes. """

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
            'id': lambda content: self.ref_binding(content[0]),
            'num': lambda content: int(content[0]),
            'quote': lambda content: content[0],
        }

        tag = expr[0]
        content = expr[1:]

        if tag in operators.keys():
            operand1 = self.interpret_expr(content[0])
            operand2 = self.interpret_expr(content[1])
            return operators[tag](operand1, operand2)

        elif tag == 'trailer':
            # Only the result of the last evaluation is returned.  That allows
            # for calculating things in the earlier runs and return only the
            # intersting part at last. => see translate_formal_species
            head = self.interpret_expr(content[0])
            body = content[1:]
            for x in content[1:]:
                key, args = x[0], x[1:]
                assert key in ('apply', 'index', 'attribute')
                head = getattr(self, key)(head, args)
            return head

        elif tag in keywords.keys():
            return keywords[tag](content)
        else:
            assert tag in ('if', 'or', 'and', 'dict', 'list', 'dna', 'uminus', 'where')
            return getattr(self, '_'+tag)(content)

    # Context-dependent
    def _if(self, content):
        """ Evaluate a (nested) if clause. """
        while len(content) > 1:
            test = self.interpret_expr(content[0])
            if test:
                value = self.interpret_expr(content[1])
                return value
            else:
                content = content[2:]
        value = self.interpret_expr(content[0])
        return value

    def _or(self, content):
        operand1 = self.interpret_expr(content[0])
        if operand1:
            return operand1
        operand2 = self.interpret_expr(content[1])
        return operand1 or operand2

    def _and(self, content):
        operand1 = self.interpret_expr(content[0])
        if not operand1:
            return False
        operand2 = self.interpret_expr(content[1])
        return operand1 and operand2

    def _dict(self, content):
        kwargs = {}
        for asign in content:
            if asign[1][0] == 'num':
                asign[1][1] = int(asign[1][1])
            if asign[1][0] == 'quote':
                asign[1][1] = asign[1][1][1:-1]
            kwargs[asign[0][1]] = asign[1][1]
        return kwargs

    def _list(self, content):
        for i in range(len(content)):
            content[i] = self.interpret_expr(content[i])
        return content

    def _dna(self, content):
        domains = content[0]
        dotparen = content[1]
        attributes = {}
        for i in range(len(domains)):
            if domains[i] != "?" and domains[i] != "+":
                starred = (len(domains[i]) == 2)
                dom = domains[i][0]
                # Get the binding for e.g.: ['id', 'd13']
                dom_value = self.interpret_expr(dom)
                attributes[dom[1]] = dom_value
                if starred:
                    dom_value = ~dom_value
                domains[i] = dom_value
        return ComplexFragment(domains, dotparen, attributes)

    def _uminus(self, content):
        value = self.interpret_expr(content[0])
        if not isinstance(value, int):
            raise NuskellEnvError("The unary minus operator can only be used with integers.")
        return -value

    def _where(self, content):
        self._create_level()
        if len(content) > 1:
            for asgn in content[1]:
                assert len(asgn) == 2 # to make the next line prettier
                id_list, value = asgn[0], asgn[1]
                value = self.interpret_expr(value)
                id_list = self._fun.remove_id_tags(id_list)

                for key, value in self._fun.asgn_pattern_match(id_list, value):
                    self.create_binding(key, value)

        value = self.interpret_expr(content[0])
        self._destroy_level()
        return value

    # trailer functions
    def apply(self, head, args):
        for i in range(len(args)):
            args[i] = self.interpret_expr(args[i])
        return self._eval_func(head, args)

    def index(self, head, args):
        subscript = self.interpret_expr(args[0])
        if not isinstance(head, list):
            raise NuskellEnvError("Only lists can be indexed.")
        if not isinstance(subscript, int):
            raise NuskellEnvError("Subscript should be an integer.")
        try:
            head = head[subscript]
        except IndexError:
            raise NuskellEnvError(
                "Error in translation scheme, expected element but got empty list.")
        return head

    def attribute(self, head, args):
        identifier = args[0][1]  # strip the tag
        if isinstance(head, ComplexFragment):
            head = head.attributes[identifier]
        elif identifier in head.__dict__:
            head = head.__dict__[identifier]
        else:
            raise NuskellEnvError(f"The attribute '{identifier}' could not be found.")
        return head

class NuskellFunctions:
    """ Builtin functions of Nuskell translation schemes.

    This object is initialized within the :obj:`NuskellEnvironment()` after the
    respective functions have been bound.  Most methods do not require any class
    variables and are therefore declared as staticmethods.  As a consequence, the
    functions can be accessed without prior initialization of the Object.
    """
    ################################
    ### Builtin-function section ###
    ################################
    def __init__(self, env):
        env.create_binding("print", NusFunction(["s"], "print"))
        env.create_binding("abort", NusFunction(["s"], "abort"))
        env.create_binding("tail", NusFunction(["l"], "tail"))
        env.create_binding("flip", NusFunction(["l", "n"], "flip"))
        env.create_binding("long", NusFunction([], "long"))
        env.create_binding("short", NusFunction([], "short"))
        env.create_binding("infty", NusFunction(["species"], "infty"))
        env.create_binding("unique", NusFunction(["l"], "unique"))
        env.create_binding("complement", NusFunction(["l"], "complement"))
        env.create_binding("rev_reactions", NusFunction(["crn"], "rev_reactions"))
        env.create_binding("irrev_reactions", NusFunction(["crn"], "irrev_reactions"))
        env.create_binding("empty", [])
 
    def eval_builtin_functions(self, f, args):
        """Evaluate built-in functions.

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
            raise NuskellEnvError(f"Function '{f}' could not be found.")

    def unique(self, args):
        """ A function that returns branch-migration domains. """
        if type(args[0]) != int:
            raise RuntimeError("The first argument of `unique' should be an integer.")
        dom = NuskellDomain(length = args[0], prefix = 'u')
        cdm = ~dom
        return dom

    def long(self, args):
        """ Returns a long domain. """
        if args:
            raise NuskellEnvError('Did not expect arguments in long function: {args}')
        dom = NuskellDomain(dtype = 'long', prefix = 'd')
        cdm = ~dom
        return dom

    def short(self, args):
        """ Returns a short domain. """
        if args:
            raise NuskellEnvError('Did not expect arguments in short function: {args}')
        dom = NuskellDomain(dtype = 'short', prefix = 't')
        cdm = ~dom
        return dom

    def _print(self, args):
        """ Print statment, primarily to debug nuskell scripts. """
        for a in args:
            if isinstance(a, ComplexFragment) or isinstance(a, NuskellComplex):
                print(list(map(str, a.sequence)))
            else:
                print(str(a), end = '')
        print()
        return void()

    def _abort(self, args):
        """ Raise NuskellExit and print the error message. """
        raise NuskellExit(args[0])

    @staticmethod
    def tail(args):
        """ Returns the list minus the first element. """
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
        """ Return a fuel complex.

        This function assigns infinite concentration to the complex, which
        is a placeholder until we have come to a more reasonable solution.
        """
        if not isinstance(args[0], ComplexFragment):
            raise NuskellEnvError("The argument of 'infty' should be a complex")

        cplxL = []
        frag = args[0]
        frag.flatten_cplx
        if frag.sequence == [] or frag.structure == []:
            assert frag.sequence == [] and frag.structure == []
            pass
        else:
            try: 
                fuel = NuskellComplex(frag.sequence, 
                                      frag.structure, 
                                      prefix = 'f')
                fuel.concentration = ('constant', float('inf'), 'nM')
                cplxL.append(fuel)
            except SingletonError as e:
                pass
        return cplxL

    @staticmethod
    def complement(args):
        """ Returns the complement of a given input. 
        
        It can interpret Complex, Domains and dot-bracket lists, as well as
        lists of lists.

        Args:
            args: An expression in form of an array, where the first element
                in the array contains the data of interest. 
        """
        x = args[0]
        if isinstance(x, list):
            # args[0] forces us to introduce additional lists ...
            for i in range(len(x)):
                x[i] = [x[i]]
            return list(reversed(map(self.complement, x)))
        elif isinstance(x, Domain):
            log.warning(f'Untested function: {~x}')
            return ~x
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
        """ Split reversible rxns into a pair of irreversible rxns. """
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
        """ Pattern matching for list assignments.

        For example: [a, b, c] = [1, 2, 3]
        """
        if isinstance(id, list):
            if not isinstance(value, list) or len(id) != len(value):
                raise NuskellEnvError(
                        f"Pattern matching failed while assigning values to {str(id)}.")
            result = []
            for i in range(len(id)):
                result += NuskellFunctions.asgn_pattern_match(id[i], value[i])
            return result
        elif isinstance(id, str):
            return [(id, value)]
        else:
            raise NuskellEnvError("Pattern matching failed.")

    @staticmethod
    def remove_id_tags(l):
        """ Helper function to remove all tags from the given id_list. """
        assert l[0] in ('idlist', 'id')
        return l[1] if l[0] == 'id' else list(map(NuskellFunctions.remove_id_tags, l[1:]))

class NuskellEnvironment(NuskellExpressions):
    """ The Nuskell language environment.

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
        self._fun = NuskellFunctions(self)

    # Public functions #
    def interpret(self, code):
        """ Setup the environment (the final namespace).

        Creates bindings for variables and functions defined in the body of the
        translation scheme. All keywords (class, function, macro, module) are
        treated the same, only the **global** keyword is special, as global
        expressions are **interpreted first** and then bound.
        """
        for stmt in code:
            kwd = stmt[0]
            body = stmt[1:]

            if kwd == "global":
                # create binding to interpretation of the expression
                id_list = body[0]
                id_list = self._fun.remove_id_tags(id_list)

                value = body[1]
                value = self.interpret_expr(value)

                for key, value in self._fun.asgn_pattern_match(id_list, 
                                                                     value):
                    self.create_binding(key, value)
            else:
                # create binding to the function, without interpretation
                # e.g. kwd = module; id = 'rxn'; args = ['r']; body_ =
                # [where...]
                assert body[0][0] == "id"
                id = body[0][1]
                # remove 'id' tags from args
                args = list(map(lambda x: x[1], body[1]))
                body_ = body[2]  # the ['where' [...]] part
                self.create_binding(id, NusFunction(args, body_))

    def translate_formal_species(self, fs_list):
        """ Apply the formal() function to the formal species in the input CRN.

        First, the bindings for formal species are created, then the formal()
        function is applied to every formal species in the input CRN.

        Returns:
          [dict()] A dictionary of key=name, value=:obj:`ComplexFragemnt()`
        """
        formal_species_objects = list(map(Species, fs_list))

        # compile the formal species
        self.create_binding("__formalspecies__", formal_species_objects)

        # map(formal, __formalspecies__)
        fs_result = self.interpret_expr(
                # tag         head           content
                #                            key       args
                ["trailer", ["id", "map"], ["apply", ["id", "formal"], 
                                                     ["id", "__formalspecies__"]]])

        self.formal_species_dict = {}
        for i in range(len(fs_list)):
            self.formal_species_dict[fs_list[i]] = fs_result[i]

        return self.formal_species_dict

    def translate_reactions(self, crn_parsed, modular=False):
        """ Execute the main() function of the translation scheme.

        The input CRN is replaced with previously initialized formal species objects.

        Args:
          crn_paresed (List[Lists]): A crn in crn_parser format.

        Raises:
          NuskellEnvError: If the compiled formal species cannot be found

        Returns:
            list: fuel species 
        """
        if not self.formal_species_dict:
            raise NuskellEnvError('Could not find the compiled formal species!')

        # replace every fs (str) with fs(NusComplex())
        crn_objects = []
        for (r, p, fw, rv) in crn_parsed:
            r = [self.formal_species_dict[x] for x in r]
            p = [self.formal_species_dict[x] for x in p]
            crn_objects.append(Reaction(r, p, rv != 0))

        # main(__crn__)
        modules = []
        if modular:
            for module in crn_objects:
                self.create_binding("__crn__", [module])
                self.constant_species_solution = self.interpret_expr(
                    ["trailer", ["id", "main"], ["apply", ["id", "__crn__"]]])
                modules.append(self.constant_species_solution)

            # NOTE: The modular reaction-by-reaction translation returns one
            # module per Reaction object.  The memory management of
            # NusComplexes ensures that species wich are shared between
            # reactions have the same name.
            unified_solution = list(set(flatten(modules)))
            modules.insert(0, unified_solution)

        else:
            self.create_binding("__crn__", crn_objects)
            self.constant_species_solution = self.interpret_expr(
                ["trailer", ["id", "main"], ["apply", ["id", "__crn__"]]])
            modules = [self.constant_species_solution]

        return modules

    def create_binding(self, name, value):
        """ Create binding of a Nuskell function in the last level.
        Args:
          name (str): Name of the function.
          value (...): anything, really.
        """
        self._env[-1][name] = value

    # Private environment modification functions #
    def _create_level(self):
        """ Create a new level for function bindings.

        This typically starts the encapulation of contents for a 'where'
        statement or a 'trailer' function.
        """
        # interpret_expr -> _where
        # interpret_expr -> _trailer -> apply -> _eval_func
        self._env.append({})

    def _destroy_level(self):
        """ Revert to previous level of function bindings.

        This commonly closes a 'where' environment or a 'trailer' function.
        """
        # interpret_expr -> _where
        # interpret_expr -> _trailer -> apply -> _eval_func
        self._env.pop()

    def ref_binding(self, fname):
        """ Search levels (reversed) for function bindings given the reference.

        Args:
            fname (str): Name of a function

        Returns:
            Function binding (self._env[?][fname])
        """
        for level in reversed(self._env):
            if fname in level.keys():
                return level[fname]
        raise NuskellEnvError(f"Cannot find a function binding for `{fname}'.")

    def _eval_func(self, f, args):
        """ Evaluate a function.

        A function (or module, class, macro, etc.) as specified in the 
        translation scheme.

        Args:
          f (:obj:`NusFunction()`): The function to evaluate.
          args (list): A list of arguments for the function.
        """
        if not isinstance(f, NusFunction):
            raise NuskellEnvError(f"`{f}' cannot be evaluated.")

        if isinstance(f.body, str): # the function is a built-in function
            return self._fun.eval_builtin_functions(f.body, args)
        else:
            if len(f.args) != len(args): # Used to be > ... why?
                raise NuskellEnvError(f"`{f.name}' requires {len(f.args)} arguments but got: {args}.")

            self._create_level()
            for i in range(len(f.args)):
                name, value = f.args[i], args[i]
                self.create_binding(name, value)
            value = self.interpret_expr(deepcopy(f.body))
            self._destroy_level()
            return value

