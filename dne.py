"""Does Not Exist
   Object to signify answer doesn't exist (e.g. division by zero)
"""

from jamenson.runtime.multimethod import defmethod, defboth_wrapper
from jamenson.runtime.atypes import anytype

from . import algebra as A


class DNEType(A.AlgebraBase):
    def __str__(self):
        return 'd.n.e.'

dne = DNEType()

@defmethod(A.mm_unop_base, [DNEType])
def meth(op):
    return dne

@defboth_wrapper(A.mm_binop_base, [DNEType, anytype])
def meth(a, b):
    return dne
