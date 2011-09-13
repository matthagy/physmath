'''algebric relationship framework using multimethods
'''

from __future__ import absolute_import

from functools import partial

try:
    from numpy import number as numpy_number_type
except ImportError:
    numpy_number_type = None

from jamenson.runtime.multimethod import MultiMethod, defmethod, defboth_wrapper
from jamenson.runtime.atypes import anytype, as_optimized_type
from jamenson.runtime.func import identity

__all__ = '''
   unop_names binop_names AlgebraBase scalar_number_type
'''.split()

scalar_number_type = as_optimized_type(tuple(filter(None, (int,long,float,numpy_number_type))))

unop_names = '''
neg pos float abs
'''.split()

__all__ += unop_names

binop_names = '''
add mul sub div truediv mod pow
lshift rshift or and
lt le eq ge gt ne
'''.split()

__all__ += binop_names

def construct_multimethods():
    '''Create all of the algebric mutlimethods (e.g. mm_add)
    '''
    gbls = globals()

    global mm_unop_base
    mm_unop_base = MultiMethod(name='mm_unop_base')

    for name in unop_names:
        mm_name = 'mm_' + name
        gbls[mm_name] = MultiMethod(name=mm_name,
                                    doc='''multimethod for unary operation %s
                                    ''' % (name,),
                                    inherit_from=[mm_unop_base])

    global mm_binop_base
    mm_binop_base = MultiMethod(name='mm_binop_base')

    for name in binop_names:
        mm_name = 'mm_' + name
        gbls[mm_name] = MultiMethod(name=mm_name,
                                    doc='''multimethod for binary operation %s
                                    ''' % (name,),
                                    inherit_from=[mm_binop_base])
        gbls['defboth_mm_' + name] = partial(defboth_wrapper, gbls[mm_name])

construct_multimethods()
del construct_multimethods


class AlgebraBase(object):
    '''Base class that translates Python algebric operations (e.g. __add__)
       into multimethod calls (e.g. mm_add)
    '''

    def construct_methods(locs=locals()):

        gbls = globals()

        def make_wrapper(name, func):
            meth_name = '__%s__' % name
            func.func_name = meth_name
            locs[meth_name] = func

        def make_unop_wrapper(name):
            mm = gbls['mm_%s' % name]
            make_wrapper(name, lambda op: mm(op))

        def make_binop_wrapper(name):
            mm = gbls['mm_%s' % name]
            make_wrapper(name, lambda lop, rop: mm(lop, rop))

        def make_binrop_wrapper(name):
            mm = gbls['mm_%s' % name]
            make_wrapper('r'+name, lambda rop, lop: mm(lop, rop))

        for name in unop_names:
            make_unop_wrapper(name)

        for name in binop_names:
            make_binop_wrapper(name)
            make_binrop_wrapper(name)

    construct_methods()
    del construct_methods


# Default multimethod behavior

@defmethod(mm_unop_base, [anytype])
def meth(op):
    return NotImplemented

@defmethod(mm_binop_base, [anytype, anytype])
def meth(lop, rop):
    return NotImplemented

@defmethod(mm_eq, [AlgebraBase, AlgebraBase])
def meth(a,b):
    return a is b

@defmethod(mm_ne, [AlgebraBase, AlgebraBase])
def meth(a,b):
    return not mm_eq(a,b)

@defmethod(mm_sub, [AlgebraBase, scalar_number_type])
def meth(a,s):
    return mm_add(a, -s)

@defmethod(mm_div, [AlgebraBase, scalar_number_type])
def meth(a,s):
    return mm_mul(a, 1.0 / float(s))


class DivAlgebraBase(AlgebraBase):
    '''Base class where truediv is delegated out to div
       (i.e. truediv is semantically equivalent to div)
    '''

@defboth_wrapper(mm_truediv, [DivAlgebraBase, anytype])
def meth(a, b):
    return mm_div(a, b)


