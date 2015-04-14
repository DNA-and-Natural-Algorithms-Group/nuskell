import unittest
# from nose.tools import *

import input
from utils import Complex, Domain, Strand, index_parts
from reactions import ReactionPathway

from condense import *

# ----------------------------------------------------------------------------
# Utils

def rsort(lst):
    return sorted(map(sorted,lst))

def pluck(dct,lst):
    return [dct[key] for key in lst]

def is_fast(reaction):
    """
    Current heuristic to determine if reaction is fast: unimolecular in reactants
    AND rate constant > k_fast
    """
    k_fast = 0.0
    return (reaction.arity == (1,1) or reaction.arity == (1,2)) and reaction.rate() > k_fast

def print_dict(d,key_format=str,value_format=str):
    print '{'
    for key in d:
        print key_format(key) + ' : ' + value_format(d[key])
    print '}'

def print_set(s,fmt=str):
    print '{ ' + ', '.join(map(fmt,s)) + ' }'

def print_condensed_reaction(reaction):
    return reaction.name + ' ' + str(map(str,reaction.reactants)) + ' -> ' + str(map(str,reaction.products))

def assert_dict_eq(d1,d2):
    print '-------------'
    print 'Comparing dicts:'
    print 'd1 = '
    print_dict(d1)
    print
    print 'd2 = '
    print_dict(d2) 
    d1_keys = d1.keys()
    d2_keys = d2.keys()
    assert sorted(d1_keys) == sorted(d2_keys)
    print '{'
    for key in d1:
        print str(key) + ': '
        print '\td1: ' + repr(d1[key]) + ', d2: ' + repr(d2[key]) 
        if not (d1[key] == d2[key]):
            d1k = d1[key]
            d2k = d2[key]
            pass
        assert d1[key] == d2[key]
    print '}'
    print '-------------'


# Small Enumerator class shim for testing
class Enum(object):
    def __init__(self,complexes,reactions):
        self.complexes = complexes
        self.reactions = reactions

class CondenseTests(unittest.TestCase):
    def setUp(self):
          
        self.domains = domains = {
                   'a':Domain('a','long'),
                   'b':Domain('b','long'),
                   'c':Domain('c','long'),
                   'd':Domain('d','long'),
                   'e':Domain('e','long'),
                   'f':Domain('f','long'),
                   'g':Domain('g','long'),
                   'h':Domain('h','long'),
                   'i':Domain('i','long'),
                   'j':Domain('j','long'),
                   }
        self.strands = strands = {
                   's1':Strand('s1',[domains['a']]),
                   's2':Strand('s2',[domains['b']]),
                   's3':Strand('s3',[domains['c']]),
                   's4':Strand('s4',[domains['d']]),
                   's5':Strand('s5',[domains['e']]),
                   's6':Strand('s6',[domains['f']]),
                   's7':Strand('s7',[domains['g']]),
                   's8':Strand('s8',[domains['h']]),
                   's9':Strand('s9',[domains['i']]),
                   's10':Strand('s10',[domains['j']]),
                   }
        self.complexes = complexes = {
                     'A':Complex('A',[strands['s1']],[[None]]),
                     'B':Complex('B',[strands['s2']],[[None]]),
                     'C':Complex('C',[strands['s3']],[[None]]),
                     'D':Complex('D',[strands['s4']],[[None]]),
                     'E':Complex('E',[strands['s5']],[[None]]),
                     'F':Complex('F',[strands['s6']],[[None]]),
                     'G':Complex('G',[strands['s7']],[[None]]),
                     'H':Complex('H',[strands['s8']],[[None]]),
                     'I':Complex('I',[strands['s9']],[[None]]),
                     'J':Complex('J',[strands['s10']],[[None]])
                     }
        self.reactions = reactions = {

                     # cycle 1: A -> B -> C -> D -> A
                     'A->B':ReactionPathway('bind11',[complexes['A']],[complexes['B']]),
                     'B->C':ReactionPathway('bind11',[complexes['B']],[complexes['C']]),
                     'C->D':ReactionPathway('bind11',[complexes['C']],[complexes['D']]),
                     'D->A':ReactionPathway('bind11',[complexes['D']],[complexes['A']]),
                     
                     # cycle 2: E -> F -> G ->E
                     'E->F':ReactionPathway('bind11',[complexes['E']],[complexes['F']]),
                     'F->G':ReactionPathway('bind11',[complexes['F']],[complexes['G']]),
                     'G->E':ReactionPathway('bind11',[complexes['G']],[complexes['E']]),
                     
                     # fast transition: cycle 1 -> cycle 2
                     'D->E':ReactionPathway('bind11',[complexes['D']],[complexes['E']]),

                     # non-deterministic choice:
                     'E->H':ReactionPathway('bind11',[complexes['E']],[complexes['H']]),
                     'E->I':ReactionPathway('bind11',[complexes['E']],[complexes['I']]),

                     # bimolecular reaction
                     'C+D->E':ReactionPathway('bind21',[complexes['C'],complexes['D']],[complexes['E']]),

                     # reversible reactions
                     'B->A':ReactionPathway('bind11',[complexes['B']],[complexes['A']]),
                     'C->B':ReactionPathway('bind11',[complexes['C']],[complexes['B']]),
                     'D->C':ReactionPathway('bind11',[complexes['D']],[complexes['C']]),
                     'A->D':ReactionPathway('bind11',[complexes['A']],[complexes['D']]),
                     'A->E':ReactionPathway('bind11',[complexes['A']],[complexes['E']]),
                     'E->A':ReactionPathway('bind11',[complexes['E']],[complexes['A']]),
                     'E->C':ReactionPathway('bind11',[complexes['E']],[complexes['C']]),
                     'E->F':ReactionPathway('bind11',[complexes['E']],[complexes['F']]),
                     'F->D':ReactionPathway('bind11',[complexes['F']],[complexes['D']])

                     }
        
        enum = Enum(complexes.values(),reactions.values())
        self.enumerator = self.enum = enum
        self.neighborhood_abcd = pluck(complexes,['A','B','C','D'])
        self.neighborhood_e = [complexes['E']]


    def testGetReactionsConsuming(self):
        complexes = pluck(self.complexes,['A','B','C','D','E','F','G'])
        print complexes
        print
        
        reactions = pluck(self.reactions,['A->B','B->C','C->D','D->A','D->E','E->F','F->G','G->E'])
        print reactions
        print
        
        enum = Enum(complexes,reactions)
        reactions_consuming = get_reactions_consuming(enum.complexes,enum.reactions)    
        expected = {
                    self.complexes['B']: [self.reactions['B->C'],], 
                    self.complexes['F']: [self.reactions['F->G'],], 
                    self.complexes['G']: [self.reactions['G->E'],], 
                    self.complexes['A']: [self.reactions['A->B'],], 
                    self.complexes['C']: [self.reactions['C->D'],], 
                    self.complexes['E']: [self.reactions['E->F'],], 
                    self.complexes['D']: [self.reactions['D->A'],self.reactions['D->E'],]
        }
        print reactions_consuming
        assert reactions_consuming == expected

    
    def testTarjans(self):
        tarjans_reactions = pluck(self.reactions,['A->B','B->C','C->D','D->A','D->E'])
        print "Reactions:"
        print tarjans_reactions
        print
        
        tarjans_complexes = self.neighborhood_abcd+self.neighborhood_e
        print "Complexes:"
        print tarjans_complexes
        print 
        
        reactions_consuming_abcd_e = get_reactions_consuming(tarjans_complexes,tarjans_reactions)
        print "Reactions consuming:"
        print reactions_consuming_abcd_e
        print 
        
        SCCs = [self.neighborhood_abcd,self.neighborhood_e]
        
        assert rsort(tarjans(tarjans_complexes,tarjans_reactions,reactions_consuming_abcd_e, is_fast)) == rsort(SCCs)

    def testTarjans2(self):
        reactions = pluck(self.reactions,['A->B','B->A','B->C','C->B','C->D','D->C','D->A'])
        complexes = self.neighborhood_abcd
        reactions_consuming = get_reactions_consuming(complexes,reactions)
        SCCs = [complexes]

        assert rsort(tarjans(complexes, reactions, reactions_consuming, is_fast)) == rsort(SCCs)

    def testTarjans3(self):
        reactions = pluck(self.reactions,['A->B','B->A','B->C','C->B','C->D','D->C','D->A', 'A->E', 'E->A', 'E->F', 'F->D'])
        complexes = pluck(self.complexes,['A','B','C','D','E','F'])
        reactions_consuming = get_reactions_consuming(complexes,reactions)
        SCCs = [complexes]

    def testTarjans4(self):
        reactions = pluck(self.reactions,['A->B','B->A','B->C','C->B','D->C','A->D','D->A', 'A->E', 'E->A', 'E->C'])
        complexes = pluck(self.complexes,['A','B','C','D','E'])
        reactions_consuming = get_reactions_consuming(complexes,reactions)
        SCCs = [complexes]

        assert rsort(tarjans(complexes, reactions, reactions_consuming, is_fast)) == rsort(SCCs)

    def testIsOutgoing(self):
        assert is_outgoing(self.reactions['D->E'],set(self.neighborhood_abcd))
        assert not is_outgoing(self.reactions['A->B'],set(self.neighborhood_abcd))
    
    def testTupleSum(self):
        # make sure a simple example works
        assert tuple_sum([(1,2,3),(4,),(5,6,7)]) == (1, 2, 3, 4, 5, 6, 7)

        # make sure there's not an extra level of summing going on
        assert tuple_sum([((1,2),3),(4,),((5,6),)]) == ((1, 2), 3, 4, (5, 6))

    def testCartesianSum(self):
        assert cartesian_sum( [ [(1,), (2,3)], [(4,5,6), (7,8)] ] ) == \
            [(1, 4, 5, 6), (1, 7, 8), (2, 3, 4, 5, 6), (2, 3, 7, 8)]
        
    # def testCartesianMultiset(self):
    #     outgoing_1_n = [[{ (1,2), (3,) }, { (4,), (5,6) }], 
    #                     [{ (7,), (8,9) }, 
    #                      { (10,), (11,) }, 
    #                      { (12,), (13,), (14,) }]]
        
    #     ms = sorted(cartesian_multisets(outgoing_1_n))
    #     expected_ms = [sorted([(1, 2, 5, 6), (1, 2, 4), (3, 5, 6), (3, 4)]), 
    #                    sorted([(7,10,12),  (7,10,13),  (7,10,14),   (7,11,12),  (7,11,13),  (7,11,14),
    #                            (8,9,10,12),(8,9,10,13),(8,9,10,14), (8,9,11,12),(8,9,11,13),(8,9,11,14) ])]
    #     print ms
    #     print expected_ms
    #     assert ms == expected_ms
    
    def testCondenseGraph(self):
        reactions = pluck(self.reactions,['A->B','B->C','C->D','D->A','D->E','E->F','F->G','G->E'])
        complexes = pluck(self.complexes,['A','B','C','D','E','F','G'])
        enum = Enum(complexes,reactions)
        out = condense_graph(enum) # dict((k,v) for (k,v) in condense_graph(enum))
        resting_states = out['resting_state_map']
        resting_state_targets = out['resting_state_targets']
        reactions = out['condensed_reactions'] 
        
        # Resting states
        print "Resting states: "
        print str(resting_states)
        
        print "Expected resting states: "
        expected_resting_states = {
            # frozenset([self.complexes['B'], self.complexes['C'], self.complexes['D'], self.complexes['A']]): RestingState('2', [self.complexes['A'], self.complexes['B'], self.complexes['C'], self.complexes['D']]), 
            frozenset([self.complexes['G'], self.complexes['E'], self.complexes['F']]): RestingState('1', [self.complexes['E'], self.complexes['F'], self.complexes['G']])
        }
        print str(expected_resting_states)
        assert_dict_eq(resting_states,expected_resting_states)
        
        # Resting state targets
        print "Resting state targets: "
        print_dict(resting_state_targets)    
        print
        
        print "Expected resting state targets: "
        resting_state = (RestingState('1', [self.complexes['E'], self.complexes['F'], self.complexes['G']]),)
        expected_resting_state_targets = {
            self.complexes['G']: SetOfFates([ resting_state ]), 
            self.complexes['E']: SetOfFates([ resting_state ]), 
            self.complexes['D']: SetOfFates([ resting_state ]), 
            self.complexes['A']: SetOfFates([ resting_state ]), 
            self.complexes['B']: SetOfFates([ resting_state ]), 
            self.complexes['F']: SetOfFates([ resting_state ]), 
            self.complexes['C']: SetOfFates([ resting_state ]),
        }
        print_dict(expected_resting_state_targets)
        assert_dict_eq(resting_state_targets, expected_resting_state_targets)
        
        # Reactions
        print "Reactions: "
        print reactions
        print
        
        print "Expected reactions:"
        expected_reactions = []
        print expected_reactions
        
        assert reactions == expected_reactions

    def testCondenseGraph2(self):
        reactions = pluck(self.reactions,['A->B','B->C','C->D','D->A','E->F','F->G','G->E','C+D->E'])
        complexes = pluck(self.complexes,['A','B','C','D','E','F','G'])
        enum = Enum(complexes,reactions)        
        out = condense_graph(enum)
        resting_states = out['resting_state_map']
        resting_state_targets = out['resting_state_targets']
        reactions = out['condensed_reactions']
        
        # Resting states: 
        print "Resting states: "
        print str(resting_states)
        print
        
        resting_state_A = (RestingState('3', [self.complexes['A'], self.complexes['B'], self.complexes['C'], self.complexes['D']]),)
        resting_state_E = (RestingState('4', [self.complexes['E'], self.complexes['F'], self.complexes['G']]),)
        
        expected_resting_states = {
            frozenset([self.complexes['B'], self.complexes['C'], self.complexes['D'], self.complexes['A']]): RestingState('3', [self.complexes['A'], self.complexes['B'], self.complexes['C'], self.complexes['D']]), 
            frozenset([self.complexes['G'], self.complexes['E'], self.complexes['F']]): RestingState('4', [self.complexes['E'], self.complexes['F'], self.complexes['G']])
        }
        assert resting_states == expected_resting_states

        # Resting state targets: 
        print "Resting state targets: "
        print_dict(resting_state_targets)    
        print

        print "Expected resting state targets: "
        expected_resting_state_targets =  {
            self.complexes['A']: SetOfFates([ resting_state_A ]), 
            self.complexes['B']: SetOfFates([ resting_state_A ]), 
            self.complexes['C']: SetOfFates([ resting_state_A ]),
            self.complexes['D']: SetOfFates([ resting_state_A ]), 
            self.complexes['E']: SetOfFates([ resting_state_E ]), 
            self.complexes['F']: SetOfFates([ resting_state_E ]), 
            self.complexes['G']: SetOfFates([ resting_state_E ]), 
        }
        print_dict(expected_resting_state_targets)
        print
        
        assert_dict_eq(resting_state_targets,expected_resting_state_targets)      
        assert resting_state_targets == expected_resting_state_targets
        
        # Reactions: 
        print "Reactions: "
        print_set(reactions,fmt=repr)
        print "Expected Reactions: "
        expected_reactions = set([ ReactionPathway('condensed',[resting_state_A[0],resting_state_A[0]],[resting_state_E[0]]) ])
        print_set(expected_reactions,fmt=repr)
        print
        assert sorted(reactions) == sorted(expected_reactions)


    # def testCondenseGraph3(self):
    #     from input import input_enum
    #     enum = input_enum('test_files/examples/sarma2010/fig4.in')
    #     enum.enumerate()
        
    #     complexes = dict()
    #     for complex in enum.complexes:
    #         complexes[complex.name] = complex
        
    #     out = condense_graph(enum)
    #     resting_states = out['resting_state_map']
    #     resting_state_targets = out['resting_state_targets']
    #     reactions = out['condensed_reactions']
        
    #     print "Detailed reactions"
    #     for r in enum.reactions:
    #         print repr(r)
    #     print
        
    #     print "Complexes"
    #     print_dict(complexes)
    #     for c in enum.complexes:
    #         print repr(c)
    #     print
        
    #     # Resting states
    #     print "Resting states"
    #     print_dict(resting_states)
    #     print
        
    #     rs = {
    #         'com2' : RestingState('com2',[complexes['com2']]),
    #         '38' : RestingState('38',[complexes['38']]),
    #         'bot' : RestingState('bot',[complexes['bot']]),
    #         'com1' : RestingState('com1',[complexes['com1']]),
    #         '39' : RestingState('39',[complexes['39']]),
    #         'top' : RestingState('top',[complexes['top']]),
    #         '15' : RestingState('15',[complexes['15']]),
    #         '40' : RestingState('40',[complexes['40']]),
    #         '18' : RestingState('18',[complexes['18']]),
    #         '25' : RestingState('25',[complexes['25']]),
    #     }
        
    #     expected_resting_states = {
    #         frozenset([complexes['com2']]) : rs['com2'],
    #         frozenset([complexes['38']]) : rs['38'],
    #         frozenset([complexes['bot']]) : rs['bot'],
    #         frozenset([complexes['com1']]) : rs['com1'],
    #         frozenset([complexes['39']]) : rs['39'],
    #         frozenset([complexes['top']]) : rs['top'],
    #         frozenset([complexes['15']]) : rs['15'],
    #         frozenset([complexes['40']]) : rs['40'],
    #         frozenset([complexes['18']]) : rs['18'],
    #         frozenset([complexes['25']]) : rs['25'],
    #     }

    #     print "Expected resting states"
    #     print_dict(expected_resting_states)
        
    #     assert expected_resting_states == resting_states
        
        
    #     # Resting state targets
    #     print "Resting state targets"
    #     print_dict(resting_state_targets)
    #     print
        
    #     expected_resting_state_targets = {
    #         complexes['55'] : SetOfFates([ ( rs['25'], rs['40'],) ]),
    #         complexes['47'] : SetOfFates([ ( rs['com1'], rs['39'],), ( rs['25'], rs['40'], rs['18'],) ]),
    #         complexes['com2'] : SetOfFates([ ( rs['com2'],) ]),
    #         complexes['38'] : SetOfFates([ ( rs['38'], ) ]),
    #         complexes['top'] : SetOfFates([ ( rs['top'], ) ]),
    #         complexes['39'] : SetOfFates([ ( rs['39'], ) ]),
    #         complexes['67'] : SetOfFates([ ( rs['25'], rs['40'] ) ]),
    #         complexes['17'] : SetOfFates([ ( rs['25'], rs['15'] ) ]),
    #         complexes['11'] : SetOfFates([ ( rs['com1'], rs['com2'] ), ( rs['25'], rs['15'], rs['18'] ) ]),
    #         complexes['44'] : SetOfFates([ ( rs['com1'], rs['39'] ), ( rs['25'], rs['40'], rs['18'] ) ]),
    #         complexes['49'] : SetOfFates([ ( rs['com1'], rs['39'] ), ( rs['25'], rs['40'], rs['18'] ) ]),
    #         complexes['bot'] : SetOfFates([ ( rs['bot'], ) ]),
    #         complexes['40'] : SetOfFates([ ( rs['40'], ) ]),
    #         complexes['23'] : SetOfFates([ ( rs['25'], rs['18'] ) ]),
    #         complexes['com1'] : SetOfFates([ ( rs['com1'], ) ]),
    #         complexes['18'] : SetOfFates([ ( rs['18'], ) ]),
    #         complexes['15'] : SetOfFates([ ( rs['15'], ) ]),
    #         complexes['8'] : SetOfFates([ ( rs['com1'], rs['com2'] ), ( rs['25'], rs['15'], rs['18'] ) ]),
    #         complexes['48'] : SetOfFates([ ( rs['com1'], rs['39'] ), ( rs['25'], rs['40'], rs['18'] ) ]),
    #         complexes['16'] : SetOfFates([ ( rs['25'], rs['18'] ) ]),
    #         complexes['25'] : SetOfFates([ ( rs['25'], ) ]),
    #     }
    #     assert resting_state_targets == expected_resting_state_targets
        
        
    #     # Condensed reactions
    #     print "Condensed Reactions"
    #     print "\n".join(map(print_condensed_reaction,reactions))
        
    #     expected_reactions = [ 
    #         ReactionPathway('condensed', [ rs['com2'], rs['bot'] ], [ rs['39'] ]),
    #         ReactionPathway('condensed', [ rs['com1'], rs['39'] ], [ rs['25'], rs['40'], rs['18'] ]),
    #         ReactionPathway('condensed', [ rs['top'], rs['bot'] ], [ rs['38'] ]),
    #         ReactionPathway('condensed', [ rs['15'], rs['bot'] ], [ rs['40'] ]),
    #         ReactionPathway('condensed', [ rs['com1'], rs['com2'] ], [ rs['25'], rs['15'], rs['18'] ]), 
    #     ]
        
    #     assert sorted(reactions) == sorted(expected_reactions)


    # def testCondenseGraph3(self):
        
    #     reactions = pluck(self.reactions,['A->B','B->C','C->D','D->A','E->F','F->G','G->E','C+D->E'])
    #     enum = Enum(self.complexes.values(),reactions)
    #     out = condense_graph(enum)
    #     resting_states = out['resting_state_map']
    #     resting_state_targets = out['resting_state_targets']
    #     reactions = out['condensed_reactions']
        
    #     # Resting states: 
    #     print "Resting states: "
    #     print str(resting_states)
    #     print
        
    #     resting_state_A = (RestingState('3', [self.complexes['A'], self.complexes['B'], self.complexes['C'], self.complexes['D']]),)
    #     resting_state_E = (RestingState('4', [self.complexes['E'], self.complexes['F'], self.complexes['G']]),)
        
    #     expected_resting_states = {
    #         frozenset([self.complexes['B'], self.complexes['C'], self.complexes['D'], self.complexes['A']]): RestingState('3', [self.complexes['A'], self.complexes['B'], self.complexes['C'], self.complexes['D']]), 
    #         frozenset([self.complexes['G'], self.complexes['E'], self.complexes['F']]): RestingState('4', [self.complexes['E'], self.complexes['F'], self.complexes['G']])
    #     }
    #     assert resting_states == expected_resting_states

    #     # Resting state targets: 
    #     print "Resting state targets: "
    #     print_dict(resting_state_targets)    
    #     print

    #     print "Expected resting state targets: "
    #     expected_resting_state_targets =  {
    #         self.complexes['A']: SetOfFates([ resting_state_A ]), 
    #         self.complexes['B']: SetOfFates([ resting_state_A ]), 
    #         self.complexes['C']: SetOfFates([ resting_state_A ]),
    #         self.complexes['D']: SetOfFates([ resting_state_A ]), 
    #         self.complexes['E']: SetOfFates([ resting_state_E ]), 
    #         self.complexes['F']: SetOfFates([ resting_state_E ]), 
    #         self.complexes['G']: SetOfFates([ resting_state_E ]), 
    #     }
    #     print_dict(expected_resting_state_targets)
    #     print
        
    #     assert_dict_eq(resting_state_targets,expected_resting_state_targets)      
    #     assert resting_state_targets == expected_resting_state_targets
        
    #     # Reactions: 
    #     print "Reactions: "
    #     print_set(reactions,fmt=repr)
    #     print "Expected Reactions: "
    #     expected_reactions = set([ ReactionPathway('condensed',[resting_state_A[0],resting_state_A[0]],[resting_state_E[0]]) ])
    #     print_set(expected_reactions,fmt=repr)
    #     print
    #     assert sorted(reactions) == sorted(expected_reactions)


    def testCondenseGraph4(self):

        # import pdb; pdb.set_trace() 

        self.bounded_dendrimer = input.input_enum('test_files/examples/bounded-dendrimer.enum')
        
        self.bounded_dendrimer.MAX_COMPLEX_SIZE = 10
        self.bounded_dendrimer.MAX_REACTION_COUNT = 1000
        self.bounded_dendrimer.MAX_COMPLEX_COUNT = 200
        self.bounded_dendrimer.RELEASE_CUTOFF = 8

        self.bounded_dendrimer.enumerate()

        out = condense_graph(self.bounded_dendrimer)
        
    def testCondenseGraph5(self):

        # import pdb; pdb.set_trace() 

        self.fate_example = input.input_pil('test_files/examples/fate-example.pil')
        
        self.fate_example.MAX_COMPLEX_SIZE = 10
        self.fate_example.MAX_REACTION_COUNT = 1000
        self.fate_example.MAX_COMPLEX_COUNT = 200
        self.fate_example.RELEASE_CUTOFF = 7

        self.fate_example.enumerate()

        out = condense_graph(self.fate_example)


    def testCondenseGraph6(self):

        self.fate_example = input.input_pil('test_files/examples/fate-example.pil')

        self.fate_example.MAX_COMPLEX_SIZE = 10
        self.fate_example.MAX_REACTION_COUNT = 1000
        self.fate_example.MAX_COMPLEX_COUNT = 200
        self.fate_example.RELEASE_CUTOFF = 7

        self.fate_example.enumerate()

        enumerator = self.fate_example
        condensed = condense_resting_states(self.fate_example)

        # Domains 
        domains = { 
            '2' : Domain('2', 8, is_complement=False, sequence='NNNNNNNN'),
            '2*' : Domain('2', 8, is_complement=True, sequence='NNNNNNNN'),
            '3' : Domain('3', 8, is_complement=False, sequence='NNNNNNNN'),
            '3*' : Domain('3', 8, is_complement=True, sequence='NNNNNNNN'),
            'a' : Domain('a', 8, is_complement=False, sequence='NNNNNNNN'),
            'a*' : Domain('a', 8, is_complement=True, sequence='NNNNNNNN'),
            'b' : Domain('b', 8, is_complement=False, sequence='NNNNNNNN'),
            'b*' : Domain('b', 8, is_complement=True, sequence='NNNNNNNN'),
            'c' : Domain('c', 8, is_complement=False, sequence='NNNNNNNN'),
            'c*' : Domain('c', 8, is_complement=True, sequence='NNNNNNNN'),
            't' : Domain('t', 4, is_complement=False, sequence='NNNN'),
            't*' : Domain('t', 4, is_complement=True, sequence='NNNN')
        }
        assert set(domains.values()) == set(enumerator.domains)

        # Strands 
        strands = { 
            '3a' : Strand('3a', [domains['3*'], domains['a*']]),
            '23' : Strand('23', [domains['2'], domains['3']]),
            'gate' : Strand('gate', [domains['a*'], domains['b*'], domains['c'], domains['b'], domains['a'], domains['2*'], domains['t*']]),
            't23' : Strand('t23', [domains['t'], domains['2'], domains['3']])
        }
        assert set(strands.values()) == set(enumerator.strands)

        # Complexes 
        complexes = { 
            '2' : Complex('2', [strands['23'], strands['3a'], strands['gate']], [[(2, 5), (1, 0)], [(0, 1), None], [(2, 4), (2, 3), None, (2, 1), (2, 0), (0, 0), None]]),
            '5' : Complex('5', [strands['23'], strands['3a'], strands['gate'], strands['t23']], [[(2, 5), (1, 0)], [(0, 1), (2, 4)], [None, (2, 3), None, (2, 1), (1, 1), (0, 0), (3, 0)], [(2, 6), None, None]]),
            '6' : Complex('6', [strands['23'], strands['3a'], strands['gate'], strands['t23']], [[(2, 5), (1, 0)], [(0, 1), None], [(2, 4), (2, 3), None, (2, 1), (2, 0), (0, 0), (3, 0)], [(2, 6), None, None]]),
            '12' : Complex('12', [strands['gate'], strands['t23']], [[(0, 4), (0, 3), None, (0, 1), (0, 0), (1, 1), (1, 0)], [(0, 6), (0, 5), None]]),
            '13' : Complex('13', [strands['23'], strands['3a']], [[None, (1, 0)], [(0, 1), None]]),
            '17' : Complex('17', [strands['23'], strands['3a'], strands['gate'], strands['t23']], [[None, (1, 0)], [(0, 1), (2, 4)], [None, (2, 3), None, (2, 1), (1, 1), (3, 1), (3, 0)], [(2, 6), (2, 5), None]]),
            '20' : Complex('20', [strands['3a'], strands['gate'], strands['t23']], [[(2, 2), (1, 4)], [None, (1, 3), None, (1, 1), (0, 1), (2, 1), (2, 0)], [(1, 6), (1, 5), (0, 0)]]),
            '21' : Complex('21', [strands['23']], [[None, None]]),
            '26' : Complex('26', [strands['3a'], strands['gate'], strands['t23']], [[(2, 2), None], [(1, 4), (1, 3), None, (1, 1), (1, 0), (2, 1), (2, 0)], [(1, 6), (1, 5), (0, 0)]]),
            'gate' : Complex('gate', [strands['23'], strands['3a'], strands['gate']], [[(2, 5), (1, 0)], [(0, 1), (2, 4)], [None, (2, 3), None, (2, 1), (1, 1), (0, 0), None]]),
            't23' : Complex('t23', [strands['t23']], [[None, None, None]])
        }
        assert set(complexes.values()) == set(enumerator.complexes)

        # Reactions 
        reactions = { 
            ReactionPathway('bind21', [complexes['gate'], complexes['t23']], [complexes['5']]),
            ReactionPathway('bind21', [complexes['2'], complexes['t23']], [complexes['6']]),
            ReactionPathway('branch_3way', [complexes['gate']], [complexes['2']]),
            ReactionPathway('branch_3way', [complexes['2']], [complexes['gate']]),
            ReactionPathway('branch_3way', [complexes['6']], [complexes['13'], complexes['12']]),
            ReactionPathway('branch_3way', [complexes['6']], [complexes['5']]),
            ReactionPathway('branch_3way', [complexes['5']], [complexes['17']]),
            ReactionPathway('branch_3way', [complexes['5']], [complexes['6']]),
            ReactionPathway('branch_3way', [complexes['17']], [complexes['21'], complexes['20']]),
            ReactionPathway('branch_3way', [complexes['17']], [complexes['13'], complexes['12']]),
            ReactionPathway('branch_3way', [complexes['17']], [complexes['5']]),
            ReactionPathway('branch_3way', [complexes['20']], [complexes['26']]),
            ReactionPathway('branch_3way', [complexes['26']], [complexes['20']]),
            ReactionPathway('open', [complexes['6']], [complexes['2'], complexes['t23']]),
            ReactionPathway('open', [complexes['5']], [complexes['gate'], complexes['t23']])
        }
        assert set(reactions) == set(enumerator.reactions)

        # Resting states 
        resting_states = { 
            '40' : RestingState('40', [complexes['12']]),
            '42' : RestingState('42', [complexes['13']]),
            '43' : RestingState('43', [complexes['21']]),
            '41' : RestingState('41', [complexes['26'], complexes['20']]),
            '38' : RestingState('38', [complexes['2'], complexes['gate']]),
            '39' : RestingState('39', [complexes['t23']])
        }
        assert set(resting_states.values()) == set(condensed['resting_states'])

        # Condensed Reactions 
        condensed_reactions = { 
            ReactionPathway('condensed', [resting_states['38'], resting_states['39']], [resting_states['43'], resting_states['41']]),
            ReactionPathway('condensed', [resting_states['38'], resting_states['39']], [resting_states['42'], resting_states['40']])
        }
        assert set(condensed_reactions) == set(condensed['reactions'])

    def testCondenseGraphCRN(self):

        complexes, reactions = dict(), set()

        def cplx(name):
            """
            Dummy function for generating formal species
            """
            name = name.strip()
            if name in complexes: 
                return complexes[name]
            else:
                complexes[name] = Complex(name, [Strand(name, [])], [])
                return complexes[name]

        def rxn(string):
            """
            Dummy function for generating reactions between formal species
            """
            reactants, products = string.split('->')
            reactants = reactants.split('+')
            products = products.split('+')

            reactants = [cplx(x) for x in reactants]
            products  = [cplx(x) for x in products]
            reactions.add(ReactionPathway('dummy', reactants, products))

        def rs(name): 
            return resting_states[cplx(name)]

        # CRN #1
        complexes, reactions = dict(), set()
        rxn('A -> B + C')
        rxn('B -> D + E')
        rxn('C -> F + G')
        enum = Enum(complexes.values(), list(reactions))

        out = condense_graph(enum) 
        resting_states = out['resting_state_map']
        resting_state_targets = out['resting_state_targets']
        reactions = out['condensed_reactions'] 
        print out
        print '----------'
        # print resting_states
        print resting_state_targets[cplx('A')]
        # assert resting_state_targets[cplx('A')] == SetOfFates([ rs('D') ])
        print "========================================================"

        complexes, reactions = dict(), set()
        rxn('A -> B')
        rxn('A -> C')
        rxn('B -> D')
        rxn('B -> E')
        rxn('C -> F')
        rxn('C -> G')
        enum = Enum(complexes.values(), list(reactions))

        out = condense_graph(enum) 
        resting_states = out['resting_state_map']
        resting_state_targets = out['resting_state_targets']
        reactions = out['condensed_reactions'] 
        print resting_state_targets[cplx('A')]

        # assert False





