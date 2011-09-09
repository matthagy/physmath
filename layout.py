from __future__ import absolute_import

import operator
from decimal import Decimal

from jamenson.runtime.multimethod import MultiMethod, defmethod
from jamenson.runtime.atypes import anytype, Seq, as_optimized_type, IsType, typep, eq_types, union
from jamenson.runtime.struct import make_struct, BaseStruct, no_default

from .types import lossless_number_type
from .sigfig import SigFig
from . import units
from .physnum import PhysNum, as_physnum
from .dne import DNEType

class BaseML(object):

    pass

class BaseMLStruct(BaseML, BaseStruct):
    pass

ml_type = as_optimized_type((BaseML, units.BaseUnit))
ml_seq_type = as_optimized_type(Seq(ml_type))

as_ml = MultiMethod('as_ml')

get_ml_json =  MultiMethod('get_ml_json')

json_ml = MultiMethod('json_ml')

@defmethod(as_ml, [ml_type])
def meth(ml):
    return ml

def defml(names, *decs):
    if isinstance(names, str):
        names = [names]
    for name in names:
        globals()[name] = struct = make_struct(name, base=BaseMLStruct, *decs)
        # set proper module to be pickleable
        struct.__module__ = __name__

@defmethod(get_ml_json, [object])
def meth(o):
    return get_ml_json(as_ml(o))

@defmethod(get_ml_json, [BaseMLStruct])
def meth(ml):
    data = json_ml(ml)
    assert isinstance(data, dict), 'bad data for %s; %r' % (ml.__class__.__name__, data)
    assert 'cls' not in data
    data['cls'] = ml.__class__.__name__
    return data

@defmethod(json_ml, [ml_type])
def meth(a):
    data = {}
    for n,d,tp in a._properties:
        v = getattr(a, n)
        if eq_types(tp, ml_type):
            assert typep(v, ml_type), 'invalid value %r for %s of %s w/ tp %s' % (v, n, a.__class__.__name__, tp)
            v = get_ml_json(v)
        elif tp == ml_seq_type:
            assert typep(v, ml_seq_type)
            v = map(get_ml_json, v)
        elif typep(v, ml_type):
            v = get_ml_json(v)
        data[n] = v
    return data


# # # # #
# null  #
# # # # #

defml("null")

@defmethod(as_ml, [IsType(None)])
def meth(n):
    return null()


# # # # #
# units #
# # # # #

unit_json_cache = {}
@defmethod(get_ml_json, [units.BaseUnit])
def meth(unit):
    try:
        return unit_json_cache[unit]
    except KeyError:
        json = unit_json_cache[unit] = make_unit_json(unit)
        return json

make_unit_json = MultiMethod('make_unit_json')

primitive_unit_unicode = {
    units.temperatures.C : [u"\u00B0C", "centigrade"],
    units.temperatures.F : [U"\u00B0F", "fahrenheit"],
    }

@defmethod(make_unit_json, [units.PrimitiveUnit])
def meth(pu):
    try:
        abbrev, name =  primitive_unit_unicode[pu]
    except KeyError:
        abbrev = unicode(pu.get_abbrev())
        name = unicode(pu.get_name())
    return dict(cls='primitive_unit', abbrev=abbrev, name=name)

prefix_unicode = {
    units.prefixes.micro : u'\u03BC',
    }
def get_prefix_unicode(p):
    try:
        return prefix_unicode[p]
    except KeyError:
        if p == units.prefixes.no_prefix:
            return u''
        return unicode(p.get_abbrev())

@defmethod(make_unit_json, [units.CompoundUnit])
def meth(cu):
    compound_unit = dict(cls='compound_unit',
                         prefix=get_prefix_unicode(cu.prefix),
                         mls_and_powers=[[get_ml_json(unit), power]
                                         for unit,power in cu.atoms_and_powers])
    prefix,name,abbrev = cu.get_name_abbrev_prefix()
    if name is not None:
        return dict(cls='named_unit',
                    prefix=get_prefix_unicode(prefix), name=unicode(name),
                    abbrev=abbrev and unicode(abbrev),
                    compound=compound_unit)
    return compound_unit

defml('error',
      ['text', no_default, unicode])


# # # # # #
# numbers #
# # # # # #

defml('dne')

@defmethod(as_ml, [DNEType])
def meth(d):
    return dne()

zero_division_error = error(u'division by zero leads to a number that does not exist (D.N.E.)')

defml('decimal',
      ['bytes', no_default, unicode],
      ['sigfigs', None, (int,long,type(None))],
      ['lsp', None, (int,long,type(None))])

@defmethod(as_ml, [(int,long)])
def meth(op):
    return decimal(unicode(op))

defml('exp10_decimal',
      ['decimal', no_default, decimal],
      ['power', no_default, (int,long)])

defml('exact',
      ['inner', no_default, ml_type])

@defmethod(as_ml, [Decimal])
def meth(d):
    if d._isnan():
        return decimal(bytes=u'NaN')
    elif d._isinfinity():
        return decimal(bytes=u'\u221E')
    elif -1000 <= d <= 1000:
        op = decimal(unicode(d))
    else:
        base,exp = SigFig(d).get_format_args()
        dec = decimal(unicode(base))
        op = dec if exp is None else exp10_decimal(dec, exp)
    return exact(op)

@defmethod(as_ml, [SigFig])
def meth(op):
    base,exp = op.get_format_args()
    dec = decimal(unicode(base), sigfigs=op.sigfigs, lsp=op.least_significant_place)
    return dec if exp is None else exp10_decimal(dec, exp)

# # # # # # # # # # # # # # #
# simple (linear) formating #
# # # # # # # # # # # # # # #

defml('constant',
      ['name', no_default, (str, unicode)])

defml('compound_name',
      ['bytes', no_default, (str, unicode)])

defml('span',
      ['mls', no_default, Seq(ml_type)])

def make_span(*args):
    return span(map(as_ml, args))

defml("crossed",
      ["ml", no_default, ml_type])

defml('boxed',
      ['ml',  no_default,  ml_type],
      ['color', u'blue',     (str, unicode)])

defml('boxed_below',
      ['ml',  no_default,  ml_type],
      ['ml_below',  no_default,  ml_type],
      ['color', u'blue',     (str, unicode)])

defml('equals',
        ['l_ml', no_default, ml_type],
        ['r_ml', no_default, ml_type],
        ['eqsign', u'=',       (str, unicode)])

defml('labeled',
      ['ml',   no_default, ml_type],
      ['label',  no_default, (str, unicode)])


# # # # # # # # # #
# physical number #
# # # # # # # # # #

defml('crossed_unit',
      ['unit', no_default, units.BaseUnit])

@defmethod(get_ml_json, [crossed_unit])
def meth(cu):
    return get_ml_json(crossed(cu.unit))

@defmethod(units.as_unit, [crossed_unit])
def meth(cu):
    return cu.unit

defml('crossed_name',
      ['name',     no_default,  (unicode,compound_name)])

@defmethod(get_ml_json, [crossed_name])
def meth(cn):
    return get_ml_json(crossed(compound_name(cn.name) if isinstance(cn.name,unicode) else cn.name))

defml('x_physnum',
      ['quantity', no_default,  union(lossless_number_type, dne)],
      ['units',    (),          Seq((units.BaseUnit, crossed_unit))],
      ['name',     None,        (type(None), unicode, compound_name, crossed_name)])

def V(quantity, unit=None, name=None, crossed_unit=False, crossed_name=False):
    if isinstance(quantity, PhysNum):
        if unit is None:
            unit = quantity.unit
        #assert unit == quantity.unit
        if name is None:
            name = quantity.name
        #assert name == quantity.name
        quantity = quantity.quantity
    if unit is not None and crossed_unit:
        unit = globals()['crossed_unit'](unit)
    if name is not None:
        if isinstance(name, str):
            name = unicode(name)
        name = compound_name(name)
        if crossed_name:
            name = globals()['crossed_name'](crossed_name)
    return x_physnum(quantity, (unit,) if unit else (), name)

@defmethod(as_ml, [PhysNum])
def meth(op):
    return x_physnum(op.quantity, (op.unit,), op.name)

@defmethod(as_physnum, [x_physnum])
def meth(op):
    return PhysNum(op.quantity,
                   reduce(operator.mul,
                          map(units.as_unit, op.units))
                   if op.units else units.dimensionless,
                   op.name)

@defmethod(get_ml_json, [x_physnum])
def meth(op):
    return get_ml_json(make_span(*
               ([op.quantity] +
                list(op.units) +
                ([compound_name(op.name) if isinstance(op.name, unicode) else op.name]
                   if op.name is not None else []))))


# # # # # # # # # # # # # # #
# Mathmatical Relationships #
# # # # # # # # # # # # # # #

defml('power',
      ["base",    no_default,  ml_type],
      ["power",   no_default,  ml_type])

defml(['pos','neg','invert'],
      ['ml',    no_default,  ml_type])

defml(['add','mul'],
      ['mls',   (),    Seq(ml_type)])

defml('group',
      ['ml',         no_default,  ml_type],
      ['bracket_type', '(',         ['(','{','[']])

defml('bar_divide',
      ['top',    no_default,  ml_type],
      ['bottom',    no_default,  ml_type])

defml('call',
      ['func', no_default, unicode],
      ['args', no_default, Seq(ml_type)])

# # # # # # # # # # #
# unit convertions  #
# # # # # # # # # # #

defml('convertion',
      ['terms',  no_default, Seq(Seq(ml_type))])

def make_convertion(*terms):
    c = convertion(terms=[])
    for term in terms:
        add_convertion_term(c, term)
    return c

def add_convertion_term(c, num, den=None):
    if den is None:
        try:
            num,den = num
        except (AttributeError, ValueError, TypeError):
            pass
    c.terms.append([as_ml(num), as_ml(den)])

@defmethod(json_ml, [convertion])
def meth(c):
    return dict(terms=[map(get_ml_json, term) for term in c.terms])

defml('equation_set',
        ['text',   no_default,  unicode],
        ['mls',    no_default,  Seq(ml_type)])

defml('calculations',
        ['mls',  no_default, Seq(ml_type)])

# # # # # # # #
# yield table #
# # # # # # # #

defml('yield_',
      ['name', no_default,     unicode],
      ['quantity', no_default, (decimal, exp10_decimal)],
      ['unit', no_default,     units.BaseUnit],
      ['limiting', False,      bool])

def make_yield(name, quantity, unit=None, limiting=False):
    if isinstance(quantity, units.PhysicalNumber):
        if unit is None:
            unit = quantity.unit
        assert unit == quantity.unit
        quantity = quantity.quantity
    return yield_(unicode(name), as_ml(quantity), as_ml(unit), limiting)

defml('yield_table',
      ['name', no_default, unicode],
      ['yields', no_default, Seq(yield_)])

@defmethod(json_ml, [yield_table])
def meth(yt):
    return dict(name=yt.name,
                yields=list(dict(name=y.name,
                                 quantity=get_ml_json(y.quantity),
                                 unit=get_ml_json(y.unit),
                                 limiting=y.limiting)
                            for y in yt.yields))

