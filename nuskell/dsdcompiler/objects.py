#
#  nuskell/dsdcompiler/objects.py
#  NuskellCompilerProject
#
import logging
log = logging.getLogger(__name__)

import gc
from dsdobjects import (SingletonError, clear_singletons)
from dsdobjects import (DomainS, ComplexS)

class NuskellDomain(DomainS):
    pass

class NuskellComplex(ComplexS):
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

def clear_memory():
    gc.collect()
    clear_singletons(NuskellDomain)
    clear_singletons(NuskellComplex)

