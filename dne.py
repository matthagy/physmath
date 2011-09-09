
from jamenson.runtime.multimethod import defmethod
from jamenson.runtime.atypes import anytype

from . import algebra as A


class DNEType(A.AlgebraBase):
    def __str__(self):
        return 'd.n.e.'

dne = DNEType()

@defmethod(A.mm_unop_base, [DNEType])
def meth(op):
    return dne

@defmethod(A.mm_binop_base, [DNEType, anytype])
@defmethod(A.mm_binop_base, [anytype, DNEType])
def meth(a, b):
    return dne
