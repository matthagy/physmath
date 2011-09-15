'''Physical numbers combine a quantity, a unit, and a possible name
'''

from __future__ import division
from __future__ import absolute_import

import re
from decimal import Decimal
from functools import partial

from hlab.bases import AutoRepr
from hlab.memorize import memorize

from jamenson.runtime import atypes
from jamenson.runtime.atypes import anytype, as_optimized_type, typep, Seq, as_type
from jamenson.runtime.multimethod import defmethod, MultiMethod
from jamenson.runtime.as_string import as_string

from . import algebra as A
from .sigfig import SigFig
from .ratio import Ratio, as_ratio
from .dne import DNEType, dne
from .units import as_unit, BaseUnit, ex_parse_unit, dimensionless

name_type = as_optimized_type((str,unicode,type(None)))
lossless_number_type = as_optimized_type((int,long,Ratio,Decimal,SigFig,DNEType))


class PhysNum(A.DivAlgebraBase, AutoRepr):

    def __init__(self, quantity, unit=None, name=None):
        assert typep(quantity, lossless_number_type), \
               'bad quantity %r(%s)' % (quantity, type(quantity).__name__)
        self.quantity = quantity
        self.unit = as_unit(unit)
        assert typep(name, name_type)
        self.name = name

    def repr_args(self):
        yield self.quantity
        if self.unit != dimensionless:
            yield self.unit
        if self.name is not None:
            if self.unit == dimensionless:
                yield None
            yield self.name

    def __str__(self):
        parts = [self.quantity]
        if self.unit != dimensionless:
            parts.append(self.unit)
        if self.name is not None:
            parts.append(self.name)
        return ' '.join(map(str, parts))

    def __hash__(self):
        if self.unit==dimensionless and not self.name:
            return hash(self.quantity)
        return hash(self.quantity) ^ hash(self.unit) ^ hash(self.name) ^ 893409342

    def split_posneg(self):
        unit = self.unit
        if isinstance(unit, PrimitiveUnit):
            unit = primunit_to_compound(unit)
        u_pos,u_neg = unit.split_posneg()
        if isinstance(self.quantity, Ratio):
            q_pos,q_neg = self.quantity.num, self.quantity.den
        else:
            q_pos,q_neg = self.quantity, 1
        return [self.__class__(q_pos, u_pos, self.name),
                self.__class__(q_neg, u_neg)]


as_physnum = MultiMethod('as_physnum')

@defmethod(as_physnum, [PhysNum])
def meth(pn):
    return pn

@defmethod(as_physnum, [lossless_number_type])
def meth(n):
    return PhysNum(n)

@defmethod(as_physnum, [str])
def meth(s):
    return parse_physical_number(s)

@defmethod(as_physnum, [unicode])
def meth(u):
    return as_physnum(str(u))

@memorize
def parse_physical_number(bytes, quantity_class=None, create_unit=False):
    parts = bytes.strip().split(None, 1)
    number, rest = (parts if len(parts)==2 else (parts[0], ''))
    unit,extra = ex_parse_unit(rest, create=create_unit) if rest else (dimensionless, None)
    name = extra and extra.strip()
    if number.endswith('d'):
        number = number[:-1]
        quantity_class = quantity_class or Decimal
    if number.endswith('s'):
        number = number[:-1]
        quantity_class = quantity_class or SigFig
    if quantity_class is None:
        quantity_class = int if not re.search('[.eE]', number) else Decimal
    return PhysNum(quantity_class(number), unit, name)

@defmethod(A.mm_eq, [PhysNum, PhysNum])
def meth(a, b):
    return (self.quantity == other.quantity and
            self.unit == other.unit and
            self.name == other.name)

@A.defboth_mm_eq([PhysNum, lossless_number_type])
def meth(p, o):
    if p.unit != dimensionless or p.name is not None:
        return NotImplemented
    return p.quantity == o

@defmethod(A.mm_neg, [PhysNum])
def meth(p):
    return PhysNum(-p.quantity, p.unit)

@defmethod(A.mm_pow, [PhysNum, (int,long)])
def meth(p,pow):
    return PhysNum((as_ratio(p.quantity)**pow
                             if pow < 0 and isinstance(p, (int,long)) else
                          p.quantity) ** pow,
                       p.unit**pow)

@defmethod(A.mm_mul, [PhysNum, PhysNum])
def meth(a, b):
    return PhysNum(a.quantity * b.quantity, a.unit * b.unit)

@A.defboth_mm_mul([PhysNum, lossless_number_type])
def meth(p, a):
    return PhysNum(p.quantity * a, p.unit)

def xdiv(a, b):
    try:
        return a/b
    except ZeroDivisionError:
        return dne
@defmethod(A.mm_div, [PhysNum, PhysNum])
def meth(a, b):
    return PhysNum(xdiv(a.quantity, b.quantity), a.unit / b.unit)

@defmethod(A.mm_div, [PhysNum, lossless_number_type])
def meth(p, a):
    return PhysNum(xdiv(p.quantity, a), p.unit)

@defmethod(A.mm_div, [anytype, PhysNum])
def meth(a, p):
    return PhysNum(xdiv(a, p.quantity), p.unit ** -1)

def add_unit_check(verb, a, b):
    if a.unit != b.unit:
        raise ValueError("cannot %s %s and %s; units are incompatible" % (verb,a,b))
    return a.unit

def add_dimensionless_check(p, other):
    if p.unit != dimensionless:
        raise ValueError("cannot add/subtract dimensionless %s and dimensional %s" % (other, p))
    return dimensionless

@defmethod(A.mm_add, [PhysNum, PhysNum])
def meth(a, b):
    return PhysNum(a.quantity + b.quantity, add_unit_check('add', a, b))

@A.defboth_mm_add([PhysNum, lossless_number_type])
def meth(p, a):
    return PhysNum(p.quantity + a, add_dimensionless_check(p, a))

@defmethod(A.mm_sub, [PhysNum, PhysNum])
def meth(a, b):
    return PhysNum(a.quantity - b.quantity, add_unit_check('sub', a, b))

@defmethod(A.mm_sub, [PhysNum, lossless_number_type])
def meth(p, a):
    return PhysNum(p.quantity - a, add_dimensionless_check(p, a))

@defmethod(A.mm_sub, [anytype, PhysNum])
def meth(a, p):
    return PhysNum(a - p.quantity, add_dimensionless_check(p, a))

@defmethod(A.mm_pow, [PhysNum, (Decimal, int, long)])
def meth(p, po):
    return PhysNum(p.quantity ** po, p.unit ** po)


# # # # # # # # # #
# Algebric Types  #
# # # # # # # # # #

class PhysNumInnerType(atypes.TypeBase):
    '''Matches inner matchers on the quantity and unit
       of a PhysNum
    '''
    __slots__ = ['quantity_inner','unit_inner']
    def __init__(self, quantity_inner=anytype, unit_inner=anytype):
        self.quantity_inner = atypes.as_type(quantity_inner)
        self.unit_inner = atypes.as_type(unit_inner)

# Type Methods
@defmethod(as_string, [PhysNumInnerType])
def meth(op):
    return '(physical_number_inner_type %s %s)' % (op.quantity_inner, op.unit_inner)

@defmethod(atypes.eq_types, [PhysNumInnerType, PhysNumInnerType])
def meth(a, b):
    return a.quantity_inner == b.quantity_inner and a.unit_inner == b.unit_inner

@defmethod(atypes.hash_type, [PhysNumInnerType])
def meth(op):
    return hash(op.quantity_inner) ^ hash(op.unit_inner) ^ 83412734

# Reductions and optimizations

# INCLUDE PAIRWISE REDUCTIONS

@defmethod(atypes.optimize_type, [PhysNumInnerType])
def meth(pn):
    quantity_inner = atypes.optimize_type(pn.quantity_inner)
    unit_inner = atypes.optimize_type(pn.unit_inner)
    if quantity_inner==anytype and unit_inner==anytype:
        return as_optimized_type((PhysNum,))
    return PhysNumInnerType(quantity_inner, unit_inner)

# typep
@defmethod(atypes.typep, [object, PhysNumInnerType])
def meth(op, pn):
    return (isinstance(op, PhysNum) and
            typep(op.quantity, pn.quantity_inner) and
            typep(op.unit, pn.unit_inner))

# Type Keyers
class PhysNumInnerKeyer(atypes.KeyerBase):
    def __init__(self, quantity_keyer, unit_keyer):
        self.quantity_keyer = quantity_keyer
        self.unit_keyer = unit_keyer
    def __eq__(self, other):
        if not isinstance(other, PhysNumInnerKeyer):
            return NotImplemented
        return (self.quantity_keyer == other.quantity_keyer and
                self.unit_keyer == other.unit_keyer)
    def __hash__(self):
        return hash(self.quantity_keyer) ^ hash(self.unit_keyer) ^ 348923432

@defmethod(atypes.get_type_keyer, [PhysNumInnerType])
def meth(op):
    return PhysNumInnerKeyer(atypes.get_type_keyer(op.quantity_inner),
                             atypes.get_type_keyer(op.unit_inner))

def physical_number_inner_keyer(quantity_keyer_func, unit_keyer_func, op):
    if not isinstance(op, PhysNum):
        return None
    return quantity_keyer_func(op.quantity), unit_keyer_func(op.unit)

@defmethod(atypes.keyer_getfunc, [PhysNumInnerKeyer])
def meth(pnk):
    return partial(physical_number_inner_keyer,
                   atypes.keyer_getfunc(pnk.quantity_keyer),
                   atypes.keyer_getfunc(pnk.unit_keyer))

def physical_number_scorer(quantity_scorer, unit_scorer, pn_key):
    if pn_key is None:
        return atypes.worst_score
    quantity_key, unit_key = pn_key
    quantity_score = quantity_scorer(quantity_key)
    unit_score = unit_scorer(unit_key)
    if atypes.no_score in [quantity_score, unit_score]:
        return atypes.no_score
    return int(round((float(quantity_score) + float(unit_score)) / 2.0))

@defmethod(atypes.get_key_scorer, [PhysNumInnerType])
def meth(pn):
    return partial(physical_number_scorer,
                   atypes.get_key_scorer(pn.quantity_inner),
                   atypes.get_key_scorer(pn.unit_inner))

