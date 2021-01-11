#
#  nuskell/objects.py
#  NuskellCompilerProject
#
import logging
log = logging.getLogger(__name__)

import gc
from itertools import chain

from dsdobjects import (SingletonError, clear_singletons, show_singletons)
from .dsdcompiler.objects import NuskellDomain, NuskellComplex
from dsdobjects import MacrostateS, ReactionS 

def clear_memory():
    # This is for unittests only!  The memory of a Singleton clears
    # automatically if there is no hard reference to the object.
    gc.collect()
    if list(show_memory()):
        log.warning('Could not collect singleton objects, trying other means.')
    clear_singletons(NuskellDomain)
    #clear_singletons(NuskellStrand)
    clear_singletons(NuskellComplex)
    clear_singletons(NuskellMacrostate)
    clear_singletons(NuskellReaction)

def show_memory():
    # This is for unittests only!  The memory of a Singleton clears
    # automatically if there is no hard reference to the object.
    for x in chain(show_singletons(NuskellDomain),
                   #show_singletons(NuskellStrand),
                   show_singletons(NuskellComplex),
                   show_singletons(NuskellMacrostate),
                   show_singletons(NuskellReaction)):
        yield x

class NuskellReaction(ReactionS):
    @property
    def name(self):
        """ str: name of the complex object. """
        return self._name

    @name.setter
    def name(self, value):
        assert self.__class__._instanceNames[self._name] == self
        assert value not in self.__class__._instanceNames
        del self.__class__._instanceNames[self._name] 
        assert self._name not in self.__class__._instanceNames
        self.__class__._instanceNames[value] = self
        self._name = value

    @property
    def auto_name(self):
        return "[{}] {} -> {}".format(self.rtype,
                                      " + ".join([x.name for x in sorted(self.reactants, 
                                                   key = lambda y: y.canonical_form)]), 
                                      " + ".join([x.name for x in sorted(self.products, 
                                                   key = lambda y: y.canonical_form)]))
 
class NuskellMacrostate(MacrostateS):
    pass


