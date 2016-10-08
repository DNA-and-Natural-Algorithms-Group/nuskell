# import unittest
# import nuskell.interpreter as ni
# from nuskell.interpreter.environment import Domain, Structure
# 
# class FunctionTesting(unittest.TestCase):
#   def setUp(self):
#     """ unit tests are run upon successful setUp """
#     pass
#   
#   def tearDown(self):
#     """ is run after the unittest (even if it was not successfull, but it
#     cleans up the *setUp* 
#     """
#     pass
# 
#   def test_flatten_func(self):
#     domains = map(Domain, [5, 5, 5])
#     dotparens = ['.', '.', '.']
#     attributes = {}
#     #attributes['d0'] = domains[0]
#     #attributes['d1'] = domains[1]
#     #attributes['d2'] = domains[2]
#     [[[['id', 'k']], '+', [['id', 'x']], '+', [['id', 'y']], [['id', 'l']], [['id', 't'], '*']], ['~', '+', '~', '+', '~', '~', '.']]
# 
# 
#     s = Structure(domains, dotparens, attributes)
#     outF = ni.flatten(s)
#     outT = [(domains[0], dotparens[0]), (domains[1], dotparens[1]), (domains[2], dotparens[2])]
#     self.assertEqual(outF, outT, 'flatten Structure object')
# 
#     outT = [('a', '('), ('b', '('), ('-a', ')')]
#     outF = ni.flatten(outT)
#     self.assertEqual(outF, outT, 'flatten flat object')
# 
#   def test_strip_consecutive_strandbreaks(self):
#     flatlist = [('+', '+'), ('a', '('), ('+', '+'), ('+', '+'), ('-a', ')'),('+', '+')]
# 
#     outF = ni.strip_consecutive_strandbreaks(flatlist)
#     outT = [('a', '('), ('+', '+'), ('-a', ')')]
#     self.assertEqual(outF, outT, 'strip consecutive strandbreaks')
# 
#   def test_rotate(self):
#     """ requires DNAObjects.Strand """
#     pass
# 
#   def test_interpret(self):
#     """ not ready yet, test environment first... """
#     pass


#class funky_Testing(unittest.TestCase):
#  def setUp(self):
#    pass
#
#  def test_shortest_cycle(self):
#    edges = [(1,2),(2,3),(2,5),(1,4),(4,5),(1,3)]
#    G = nx.Graph()
#    G.add_edges_from(edges)
#
#    sc = fug.shortest_cycle(G, 3)
#    self.assertEqual(sorted(sc), [1,2,3])
#    sc = fug.shortest_cycle(G, 1)
#    self.assertEqual(sorted(sc), [1,2,3])
#    sc = fug.shortest_cycle(G, 2)
#    self.assertEqual(sorted(sc), [1,2,3])
#
#    sc = fug.shortest_cycle(G, 4)
#    self.assertEqual(sorted(sc), [1,2,4,5])
#    sc = fug.shortest_cycle(G, 5)
#    self.assertEqual(sorted(sc), [1,2,4,5])
#
#class WidgetTestCase(unittest.TestCase):
#  def setUp(self):
#    """ unit tests are run upon successful setUp """
#    self.widget = Widget('The widget')
#
#  def tearDown(self):
#    """ is run after the unittest (even if it was not successfull, but it
#    cleans up the *setUp* 
#    """
#    self.widget.dispose()
#    self.widget = None
#
#  def test_default_size(self):
#    self.assertEqual(self.widget.size(), (50,50),
#        'incorrect default size')
#
#  def test_resize(self):
#    self.widget.resize(100,150)
#    self.assertEqual(self.widget.size(), (100,150), 'wrong size after resize')
#
## As an alternative to the functions specified above
# class DefaultWidgetSizeTestCase(SimpleWidgetTestCase):
#   def runTest(self):
#     self.assertEqual(self.widget.size(), (50,50),
#         'incorrect default size')
# 
# class WidgetResizeTestCase(SimpleWidgetTestCase):
#   def runTest(self):
#     self.widget.resize(100,150)
#     self.assertEqual(self.widget.size(), (100,150), 'wrong size after resize')


#def suite():
#  suite = unittest.TestSuite()
#  suite.addTest(WidgetTestCase('test_default_size'))
#  suite.addTest(WidgetTestCase('test_resize'))
#  return suite
 

if __name__ == '__main__':
  unittest.main()

