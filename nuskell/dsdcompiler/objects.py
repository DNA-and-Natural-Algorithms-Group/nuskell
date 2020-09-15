#
#  nuskell/dsdcompiler/objects.py
#  NuskellCompilerProject
#
import logging
log = logging.getLogger(__name__)

from dsdobjects import (clear_memory,
                        DSDObjectsError,
                        DSDDuplicationError)
from dsdobjects.core import DL_Domain, DSD_Complex
from dsdobjects.prototypes import Complex
from dsdobjects.utils import natural_sort
#import dsdobjects.objectio as oio

class NuskellObjectError(Exception):
    pass

class NuskellDomain(DL_Domain):
    """ Nucleic acid domain sequence.

    Inherits from :obj:`dsdobjects.base_classes.DL_Domain()`.

    A domain is a sequence of consecutive nucleotides. Sequences of domains can
    form the same secondary structures as sequences of nucleotides.  Several
    options for specifying domain properties are allowed. Domains might have an
    explicit integer (bp) length, or may be designated as short or long. If the
    latter method is used, the code will use the relevant constant as the
    integer domain length.

    Globals:
      ID (int): An automatically assigned ID that is used to find automated
      species names.

    Args:
      name (str, optional): Name of this domain. If not specified, an automatic
        name is generated, defaults to ''
      prefix (str, optional): A prefix for automated naming of Domains. Has to be
        set if no name is provided. Defaults to 'd' for 'domain'. Usually, the
        prefix 't' is used for 'toehold domains', and 'h' for 'histoy domains'.
      dtype (bool, optional): One of *short* or *long*. Defaults to None, i.e. 
        guess from length parameter.
      length (int, optional): Set length of a domain. Defaults to None, i.e. set
        default length for *dtype* parameter.

    Raises:
      NuskellObjectError: Domain prefix must not be empty!
      NuskellObjectError: NuskellDomain prefix must not be empty!
      NuskellObjectError: NuskellDomain prefix must not end with a digit!

    """
    ID = 0          # ID is used to assign names automatically
    def __new__(cls, name = '', prefix = 'd', dtype = None, length = None):
        # The new method returns the present instance of an object, if it exists
        self = DL_Domain.__new__(cls)
        if name == '':
            if prefix == '':
                raise NuskellObjectError(
                        'NuskellDomain prefix must not be empty!')
            elif prefix[-1].isdigit():
                raise NuskellObjectError(
                        'NuskellDomain prefix must not end with a digit!')
            name = prefix + str(NuskellDomain.ID)
            NuskellDomain.ID += 1
        try:
            super(NuskellDomain, self).__init__(name, dtype, length)
        except DSDDuplicationError as e:
            other = e.existing
            if dtype and (other.dtype != dtype) :
                raise DSDObjectsError('Conflicting dtype assignments',
                        f'for {name}: "{dtype}" vs. "{other.dtype}"')
            elif length and (other.length != length) :
                raise DSDObjectsError('Conflicting length assignments',
                        f'for {name}: "{length}" vs. "{other.length}"')
            return e.existing
        self.nucleotides = None
        return self

    def __init__(self, name='', prefix='d', dtype=None, length=None):
        # Remove default initialziation to get __new__ to work
        pass

    @property
    def complement(self):
        """ NuskellDomain: retuns the complementary domain object. """
        if self._complement is None:
            cname = self._name[:-1] if self.is_complement else self._name + '*'
            if cname in DL_Domain.MEMORY:
                self._complement = DL_Domain.MEMORY[cname]
            else :
                self._complement = NuskellDomain(name=cname, dtype=self.dtype, length=self.length)
        return self._complement

    @property
    def is_complement(self):
        """ bool: The domain is a complement: (e.g. A* rather than A). """
        return self._name[-1:] == '*'

class NuskellComplex(Complex):
    pass

