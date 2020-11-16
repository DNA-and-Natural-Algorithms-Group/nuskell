#
#  nuskell/dsdcompiler/objects.py
#  NuskellCompilerProject
#
import logging
log = logging.getLogger(__name__)

import gc
from dsdobjects import (SingletonError, clear_singletons)
from dsdobjects import DomainS as NuskellDomain
from dsdobjects import ComplexS as NuskellComplex

def clear_memory():
    gc.collect()
    clear_singletons(NuskellDomain)
    clear_singletons(NuskellComplex)

