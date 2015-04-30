#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#
#
# The interpreter for translation schemes.
#

import ts_parser 
from copy import copy
import DNAObjects

# default lengths for short() and long() built-in functions
short_domain_length = 6
long_domain_length = 15

# built-in types for the ts language
class Species:
    def __init__(self, name):
        self.name = name

class Reaction:
    def __init__(self, reactants, products, reversible):
        self.reactants = reactants
        self.products = products
        self.reversible = reversible

class Function:
    def __init__(self, args, body):
        self.args = args
        self.body = body

class Domain:
    domain_id = 0
    domain_length = {}
    def __init__(self, length, id = None):
        if id == None:
            Domain.domain_id += 1
            id = Domain.domain_id
            Domain.domain_length[id] = length
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
        if not isinstance(other, Structure): return False
        return self.id == other.id
    def __ne__(self, other):
        return not (self == other)

class Structure:
    def __init__(self, domains, dotparens, attr):
        self.domains = domains
        self.dotparens = dotparens
        self.attributes = attr

class Solution:
    def __init__(self, initial):
        self.molecules = initial
    def __add__(self, other):
        return Solution(self.molecules.union(other.molecules))

# exceptions
class RuntimeError(Exception):
    def __init__(self, value):
        self.value = "Runtime error : " + value
    def __str__(self):
        return self.value

# environment management functions
def env_init():
    return [{}]

def env_create_binding(env, name, value):
    env[-1][name] = value
    return env

def env_ref_binding(env, name):
    for level in reversed(env):
        if name in level.keys():
            return level[name]
    raise RuntimeError("Cannot find a binding for `" + name + "'.")

def env_create_level(env):
    env.append({})
    return env

def env_destroy_level(env):
    env.pop()
    return env

def interpret_main(env, code):
    """Creates bindings for variables and functions defined in the body of
       the code. Returns the environment (the final namespace).
    """
    for stmt in code:
        kwd = stmt[0]
        body = stmt[1:]
        if kwd == "global":
            id_list = body[0]
            value = body[1]
            env, value = interpret_expr(env, value)
            asgn_list = asgn_pattern_match(trim_id_list(id_list), value)
            for key, value in asgn_list:
                env = env_create_binding(env, key, value)
        else:
            # function declaration
            assert body[0][0] == "id"
            id = body[0][1]
            args = map(lambda x: x[1], body[1])
            body_ = body[2]
            env = env_create_binding(env, id, Function(args, body_))
    return env

def asgn_pattern_match(id, value):
    """Does the pattern matching for list assignments.
       ex) [a, b, c] = [1, 2, 3]
    """
    if type(id) == list:
        if type(value) != list or len(id) != len(value):
            raise RuntimeError("Pattern matching failed while assigning \
values to " + str(id) + ".")
        result = []
        for i in range(len(id)):
            result += asgn_pattern_match(id[i], value[i])
        return result
    elif type(id) == str:
        return [(id, value)]
    else:
        raise RuntimeError("The left hand side of an assignment should be \
either an identifier or a list-tree consisting only of identifiers.")

def trim_id_list(l):
    """Removes all the tags from the given id_list."""
    kwd = l[0]
    if kwd == "idlist":
        return map(trim_id_list, l[1:])
    elif kwd == "id":
        return l[1]
    else:
        raise RuntimeError("The left hand side of an assignment should be \
either an identifier or a list-tree consisting only of identifiers.")

def interpret_expr(env, v):
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

    tag = v[0]
    content = v[1:]
    if tag == "id":
        return env, env_ref_binding(env, content[0])
    elif tag == "num":
        return env, int(content[0])
    elif tag == "list":
        for i in range(len(content)):
            env, content[i] = interpret_expr(env, content[i])
        return env, content
    elif tag == "dna":
        domains = content[0]
        dotparen = content[1]
        attributes = {}
        for i in range(len(domains)):
            if domains[i] != "?" and domains[i] != "+":
                starred = (len(domains[i]) == 2)
                dom = domains[i][0]
                env, dom_value = interpret_expr(env, dom)        
                attributes[dom[1]] = dom_value
                if starred:
                    dom_value = ~dom_value
                domains[i] = dom_value
        return env, Structure(domains, dotparen, attributes)
    elif tag == "trailer":
        # function call, list indexing, or accessing attribute of an object
        env, head = interpret_expr(env, content[0])
        for x in content[1:]:
            if x[0] == "apply":
                arguments = x[1:]
                for i in range(len(arguments)):
                    env, arguments[i] = interpret_expr(env, arguments[i])
                env, head = eval_function(env, head, arguments)
            elif x[0] == "index":
                env, subscript = interpret_expr(env, x[1])
                if type(head) != list:
                    raise RuntimeError("Only lists can be indexed.")
                if type(subscript) != int:
                    raise RuntimeError("Subscript should be an integer.")
                head = head[subscript]
            elif x[0] == "attribute":
                identifier = x[1][1] # strip the tag
                if isinstance(head, Structure):
                    head = head.attributes[identifier]
                elif identifier in head.__dict__:
                    head = head.__dict__[identifier]
                else:
                    raise RuntimeError("The attribute `"+identifier+"' \
could not be found.")
        return env, head
    elif tag == "uminus":
        env, value = interpret_expr(env, content[0])
        if type(value) != int:
            raise RuntimeError("The unary minus operator can only be used \
with integers.")
        return env, -value
    elif tag in operators.keys():
        env, operand1 = interpret_expr(env, content[0])
        env, operand2 = interpret_expr(env, content[1])
        return env, operators[tag](operand1, operand2)
    elif tag == "and":
        env, operand1 = interpret_expr(env, content[0])
        if not operand1:
            return env, False
        env, operand2 = interpret_expr(env, content[1])
        return env, operand1 and operand2
    elif tag == "or":
        env, operand1 = interpret_expr(env, content[0])
        if operand1:
            return env, operand1
        env, operand2 = interpret_expr(env, content[1])
        return env, operand1 or operand2
    elif tag == "where":
        env = env_create_level(env)
        if len(content) > 1:
            # if this is not a trivial where clause 
            for asgn in content[1]:
                id_list = asgn[0]
                value = asgn[1]
                env, value = interpret_expr(env, value)
                asgn_list = asgn_pattern_match(trim_id_list(id_list), value)
                for key, value in asgn_list:
                    env = env_create_binding(env, key, value)
            # evaluate the final value
        env, value = interpret_expr(env, content[0])
        env = env_destroy_level(env)
        return env, value
    elif tag == "if":
        while len(content) > 1:
            env, test = interpret_expr(env, content[0])
            if test:
                env, value = interpret_expr(env, content[1])
                return env, value
            else:
                content = content[2:]
        env, value = interpret_expr(env, content[0])
        return env, value
    else:
        raise RuntimeError("Unknown expression `" + tag + "' was found.")
        return env, None

def eval_function(env, f, args):
    if type(f.body) == str: # the function is a built-in function
        return eval_builtin_function(env, f.body, args)

    def hardcopy_list(l):
        if type(l) != list:
            return copy(l)
        return map(hardcopy_list, l)

    env = env_create_level(env)

    if not isinstance(f, Function):
        raise RuntimeError(str(f) + "is not a function.")

    if len(f.args) > len(args):
        raise RuntimeError("The function `" + f.name + "' requires at \
least " + str(len(f.args)) + " arguments but only found " + str(args) +\
" arguments.")

    for i in range(len(f.args)):
        arg_name = f.args[i]
        arg_value = args[i]
        env = env_create_binding(env, arg_name, arg_value)

    # hardcopy the function body so that it does not change during
    # evaluation.
    env, value = interpret_expr(env, hardcopy_list(f.body))
    env = env_destroy_level(env)
    return env, value

def eval_builtin_function(env, f, args):
    if f == "tail":
        if type(args[0]) != list:
            raise RuntimeError("`tail' should have a list as its argument.")
        return env, args[0][1:]
    elif f == "complement":
        def complement(x):
            if type(x) == list:
                return reversed(map(complement, x))
            elif isinstance(x, Domain):
                return ~x
            elif isinstance(x, Structure):
                return Structure(complement(x.domains),
                                 complement(x.dotparens),
                                 dict(x.attributes))
            elif x == "(":
                return ")"
            elif x == ")":
                return "("
            else:
                return x
        return env, complement(args[0])
    elif f == "infty":
        if not isinstance(args[0], Structure):
            raise RuntimeError("The argument of `infty' should be a \
structure")
        return env, Solution(set([args[0]]))
    elif f == "short":
        return env, Domain(short_domain_length)
    elif f == "long":
        return env, Domain(long_domain_length)
    elif f == "unique":
        if type(args[0]) != int:
            raise RuntimeError("The first argument of `unique' should be an\
 integer.")
        return env, Domain(args[0])
    elif f == "flip":
        if type(args[0]) != list:
            raise RuntimeError("The first argument of `flip' should be a\
 list.")
        if type(args[1]) != int:
            raise RuntimeError("The second argument of `flip' should be an\
 integer.")
        return env, flip(args[0], args[1])
    elif f == "rev_reactions":
        if type(args[0]) != list:
            raise RuntimeError("The argument of `rev_reactions' should be a\
 list.")
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
        return env, new_crn
    elif f == "irrev_reactions":
        if type(args[0]) != list:
            raise RuntimeError("The argument of `irrev_reactions' should be a\
 list.")
        crn = args[0]
        new_crn = []
        for r in crn:
            if r.reversible:
                new_crn.append(Reaction(r.reactants, r.products, False))
                new_crn.append(Reaction(r.products, r.reactants, False))
            else:
                new_crn.append(r) 
        return env, new_crn
        
    raise RuntimeError("`" + f + "' could not be resolved.")

# the `flip' built-in function
def flip(l, n):
    res = [0] * n
    for i in range(n):
        res[i] = [0] * len(l)
    for i in range(len(l)):
        if len(l[i]) != n:
            raise RuntimeError("The argument of `flip' does not have a \
required form.");
        for j in range(n):
            res[j][i] = l[i][j]
    return res

def init_builtin_functions(env):
    env = env_create_binding(env, "tail", Function(["l"], "tail"))
    env = env_create_binding(env, "complement", \
                             Function(["l"], "complement"))
    env = env_create_binding(env, "infty", Function(["species"], "infty"))
    env = env_create_binding(env, "short", Function([], "short"))
    env = env_create_binding(env, "long", Function([], "long"))
    env = env_create_binding(env, "unique", Function(["l"], "unique"))
    env = env_create_binding(env, "flip", Function(["l", "n"], "flip"))
    env = env_create_binding(env, "empty", Solution(set()))
    env = env_create_binding(env, "rev_reactions", \
                             Function(["crn"], "rev_reactions"))
    env = env_create_binding(env, "irrev_reactions", \
                             Function(["crn"], "irrev_reactions"))
    return env

def flatten(x):
    """Takes a Structure instance and convert it to the following format:
       [(a, "("), (b, "("), ..., (-a, ")")]
    """
    if isinstance(x, Structure):
        l = flip([x.domains, x.dotparens], len(x.domains))
        l = map(lambda x: (x[0], x[1]), l)
        return flatten(l)
    if type(x) == list:
        def concat(l):
            res = []
            for x in l:
                res += x
            return res
        return concat(map(flatten, x))
    else:
        x, y = x
        if y != "~":
            return [(x, y)]
        else:
            return flatten(x)

def strip_consecutive_strandbreaks(l):
    flag = 1
    res = []
    for (x, y) in l:
        if x == "+" and flag == 1: continue
        if x == "+":
            flag = 1
        else:
            flag = 0
        res.append((x, y))
    (x, y) = res[-1]
    if x == "+": res = res[:-1]
    return res

def rotate(complex):
    def find(l, key):
        for i in range(len(l)):
            if l[i] == key:
                return i
        return None

    def hardcopyList(l):
        if type(l) != list:
            return l
        return map(hardcopyList, l)

    complex = hardcopyList(complex)

    if "+" not in complex[1]:
        return complex        
    else:
        dom = complex[0][1:] + complex[0][:1]

        # change parentheses appropriately
        p = find(complex[1], "+")
        dpr = complex[1]
        stack = []
        for i in range(p):
            if dpr[i] == "(": stack.append(i)
            elif dpr[i] == ")": stack.pop()
        for i in stack:
            dpr[i] = ")"
        stack = []
        for i in reversed(range(p + 1, len(dpr))):
            if dpr[i] == ")": stack.append(i)
            elif dpr[i] == "(": stack.pop()
        for i in stack:
            dpr[i] = "("
        
        dpr = dpr[p + 1:] + ["+"] + dpr[:p]
        return [dom, dpr]

def interpret(translation_scheme, crn, formal_species):
    # initialize the namespace
    env = env_init()

    # initialize the built-in functions
    env = init_builtin_functions(env)
    header = ts_parser.parse("""
function range(x) = if x == 0 then [] else range(x - 1) + [x - 1] ;
function sum(x) = if len(x) == 0 then empty elseif len(x) == 1 then x[0] else x[0] + sum(tail(x)) ;
function len(x) = if x == [] then 0 else 1 + len(tail(x)) ;
function reverse(x) = if x == [] then [] else reverse(tail(x)) + [x[0]] ;
function rxn_degree(x, r) = if len(x) == 0 then [] elseif len(x[0].reactants) == r then [x[0]] + rxn_degree(tail(x), r) else rxn_degree(tail(x), r) ;
function unirxn(x) = if len(x) == 0 then [] elseif len(x[0].reactants) == 1 then [x[0]] + unirxn(tail(x)) else unirxn(tail(x)) ;
function birxn(x) = if len(x) == 0 then [] elseif len(x[0].reactants) == 2 then [x[0]] + birxn(tail(x)) else birxn(tail(x)) ;
function map(f, x) = if len(x) == 0 then [] else [f(x[0])] + map(f, tail(x)) ;
function map2(f, y, x) = if len(x) == 0 then [] else [f(y, x[0])] + map2(f, y, tail(x))
    """)
    env = interpret_main(env, header)

    # create formal species
    formal_species_objects = map(Species, formal_species)

    # run the translation scheme
    env = interpret_main(env, translation_scheme)

    # compile the formal species
    env_create_binding(env, "__formalspecies__", formal_species_objects)
    env, formal_species_result = interpret_expr(env, ["trailer", \
        ["id", "map"], ["apply", ["id", "formal"], ["id", \
        "__formalspecies__"]]])
    # final tranlsation of formal species
    formal_species_dict = {}
    for i in range(len(formal_species)):
        formal_species_dict[formal_species[i]] = formal_species_result[i]

    # compile reactions
    crn_remap = map(lambda x: [x[0]] + map(\
         lambda y: map(lambda z: formal_species_dict[z], y), x[1:]), crn)
    crn_object = map(lambda x: Reaction(x[1], x[2], x[0] == "reversible"), crn_remap)
    env_create_binding(env, "__crn__", crn_object)
    env, solution = interpret_expr(env, \
          ["trailer", ["id", "main"], ["apply", ["id", "__crn__"]]])

    # flatten outputs
    for k in formal_species_dict.keys():
        formal_species_dict[k] = \
          strip_consecutive_strandbreaks(flatten(formal_species_dict[k]))

    solution_as_list = list()
    for m in solution.molecules:
        solution_as_list.append(strip_consecutive_strandbreaks(flatten(m)))

    # convert output to DNA Objects
    domains = []
    strands = []
    formal_species = []
    constant_species = []

    def get_domain(x):
        domain_name = str(x)
        starred = False
        if domain_name[-1] == "*":
            domain_name = domain_name[:-1]
            starred = True
        for d in domains:
            if d.name == domain_name:
                if starred: return d.C
                return d
        new_dom = DNAObjects.Domain(name = domain_name, length = x.length)
        domains.append(new_dom)
        if starred: new_dom = new_dom.C
        return new_dom

    def get_strand(strand):
        for x in strands:
            if x.domain_list == strand:
                return x
        new_strand = DNAObjects.Strand(domains = strand)
        strands.append(new_strand)
        return new_strand

    wildcard_domain = DNAObjects.Domain(name = "?")

    # convert formal species
    for fs_name in formal_species_dict:
        # add a strand break at the end for convenience
        complex_old_format = formal_species_dict[fs_name] + [("+", "+")]
        complex = []
        strand = []
        structure = ""
        for (x, y) in complex_old_format:
            if y == "?": y = "."
            structure += y
            if x == "+":
                strand = get_strand(strand)
                complex.append(strand)
                strand = []
                continue
            if x == "?":
                x = wildcard_domain
            else:
                x = get_domain(x)
            strand.append(x)
        # remove the strand break that was added at the beginning
        structure = structure[:-1]
        formal_species.append( \
            DNAObjects.Complex(name = fs_name, \
                               strands = complex, \
                               structure = structure))

    previous = []
    # convert constant species
    for cs in solution_as_list:
        # add a strand break at the end for convenience
        complex_old_format = cs + [("+", "+")]
        complex = []
        strand = []
        structure = ""
        for (x, y) in complex_old_format:
            if y == "?": y = "."
            structure += y
            if x == "+":
                strand = get_strand(strand)
                complex.append(strand)
                strand = []
                continue
            if x == "?":
                x = wildcard_domain
            else:
                x = get_domain(x)
            strand.append(x)
        # remove the strand break that was added at the beginning
        structure = structure[:-1]
        # has it been already added?
        flag = False
        for y in previous:
            x = [complex, list(structure)]
            p = rotate(x)
            while True:
                if p == y:
                    flag = True
                if p == x:
                    break
                p = rotate(p)
        if not flag:
            previous.append([complex, list(structure)])
            constant_species.append( \
                DNAObjects.Complex(strands = complex, \
                                   structure = structure))

    return domains, strands, formal_species, constant_species
