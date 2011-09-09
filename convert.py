'''Utilities for performing unit convertions
'''

import re
import operator
from decimal import Decimal

from jamenson.runtime.multimethod import defmethod, MultiMethod
from jamenson.runtime.atypes import as_optimized_type, typep, anytype

from .sigfig import SigFig
from .dne import dne
from . import units as U
from . import layout
from . import physnum
from .physnum import PhysNum, as_physnum
from .layout import V
from .annotator import annotator
from .types import lossless_number_type

convert = MultiMethod('convert',
                      '''
                      High level unit converter
                      ''')

@defmethod(convert, [PhysNum, U.BaseUnit])
def meth(num, to_unit):
    '''High level unit converter
    '''
    if num.unit == to_unit:
        return num
    if num.unit.without_prefix() == to_unit.without_prefix():
        return convert_unit_prefix(num, to_unit)
    return x_convert(num, to_unit)

def convert_factor(number, num, den=None, power=1):
    '''Convert using a series of conversion factors
    '''
    return Converter(number).factor_convert(num, den, power).finish()

def convert_unit_prefix(num, to_unit):
    '''perform power convertions (convertion of unit prefix)
         ie. km -> mm or mol -> nmol
    '''
    if num.unit == to_unit:
        return num
    return Converter(num).prefix_convert(to_unit).finish()

def convert_by_path(num, to_unit):
    '''Perform convertions using a path of conversion factors
         ie. m -> cm -> in -> ft
       Also handles prefix convertions as needed
    '''
    if num.unit == to_unit:
        return run
    return Converter(num).path_convert(to_unit).finish()

x_convert = MultiMethod('x_convert')

@defmethod(x_convert, [anytype, anytype])
def meth(num, unit):
    '''Fallback on convert_by_path when nothing more specific
    '''
    return convert_by_path(num, unit)


# # # # # # # # # # # # # #
# Dimensional Convertions #
# # # # # # # # # # # # # #

def defdimconvert(dim):
    dimt = U.UnitDimensionalityType(dim)
    def wrap(func):
        return defmethod(x_convert, [physnum.PhysNumInnerType(unit_inner=dimt), dimt])(func)
    return wrap

def temperature_factors():
    F = U.temperatures.F
    C = U.temperatures.C
    K = U.temperatures.K
    return {(K,C): [(1,1), '-273.15'],
            (K,F): [(9,5), '-459.67'],
            (C,F): [(9,5), '32']}
temperature_factors = temperature_factors()

@defdimconvert('temperature')
def meth(num, to_unit):
    '''Convert between different temperatures
    '''
    try:
        (mn,md),b = temperature_factors[num.unit, to_unit]
        b = Decimal(b)
    except KeyError:
        (mn,md),b = temperature_factors[to_unit, num.unit]
        b = Decimal(b)
        b = -md * b / mn
        mn,md = md,mn
    try:
        result = PhysNum(num.quantity * Decimal(mn) / Decimal(md) + b, to_unit)
    except ZeroDivisionError:
        result = dne
    if annotator.annotating:
        annotator.annotate(
          layout.equals(
            layout.add([layout.make_convertion(
                                [V(num, crossed_unit=True),
                                 None],
                                [V(mn, to_unit),
                                 V(md, num.unit, crossed_unit=True)]),
                        V(b, to_unit) if b>0 else layout.neg(V(abs(b), unit))]),
            V(result)))
        if result is dne:
            annotator.annotate(layout.zero_division_error)
    return result

def volume_systems():
    data = '''
    liter: L
    metric: m3
    imperial_gas: in3 ft3 yd3
    imperial_liquid: oz pt qt gal
    '''
    for line in data.split('\n'):
        line = line.strip()
        if not line:
            continue
        system, units = line.split(':')
        for unit in map(U.parse_unit, units.strip().split()):
            yield unit, system
volume_systems = dict(volume_systems())

@defdimconvert('volume')
def meth(num, to_unit):
    '''Convert between different volumes.
       This algorithm is a bit a heristics hack, but likely the most logical way to handle this insanity.
    '''

    from_system = volume_systems[num.unit.cannonicalized().without_prefix()]
    to_system = volume_systems[to_unit.cannonicalized().without_prefix()]

    # definitions of howto go from on system to another using the
    # following function-based shorthand
    cnv = Converter(num)
    factor = cnv.factor_convert
    prefix = cnv.prefix_convert
    path = cnv.path_convert

    def rfactor(a,b,power=1): return factor(b,a,power)
    def liter_cc(): prefix('mL'); rfactor('1 mL', '1 cc')
    def cc_liter(): prefix('cm3', 3); rfactor('1 cc', '1 mL')
    def cc_in(): prefix('cm3', 3); rfactor('2.54 cm', '1 in', 3)
    def in_cc(): prefix('in3', 3); rfactor('1 in', '2.54 cm', 3)
    def liter_in(): liter_cc(); cc_in()
    def mL_il(to_unit=to_unit): prefix('mL'), rfactor('28.41 mL', '1 floz'), path(to_unit)
    def ig_cc(): path('in3', 3), in_cc()
    def il_in3(): path('floz'), rfactor('1 floz', '1.739 in3')
    def il_cc(): il_in3(), in_cc()

    {'ll': lambda : prefix(to_unit),
     'lm': lambda : (liter_cc(), prefix(to_unit, 3)),
     'lig': lambda : (liter_in(), path(to_unit, 3)),
     'lil': lambda : mL_il(),

     'mm': lambda : prefix(to_unit, 3),
     'ml': lambda : (cc_liter(), prefix(to_unit)),
     'mig': lambda : (cc_in(), path(to_unit, 3)),
     'mil': lambda : (cc_liter(), mL_il()),

     'igl': lambda : (ig_cc(), cc_liter(), path(to_unit)),
     'igm': lambda : (ig_cc(), prefix(to_unit, 3)),
     'igig': lambda : prefix(to_unit, 3),
     'igil': lambda : (path('in3', 3), rfactor('1.739 in3', '1 floz'), path(to_unit)),

     'ill': lambda : (il_cc(), cc_liter()),
     'ilm': lambda : (il_cc(), path(to_unit, 3)),
     'ilig' : lambda : (il_in3(), path(to_unit, 3)),
     'ilil' : lambda : prefix(to_unit)
     }[''.join(x[0] for name in [from_system, to_system] for x in name.split('_'))
       ]() #exectue the proper convertion
    return cnv.finish()


def sign(op):
    if op>0: return 1
    if op<0: return -1
    return 0

class Converter(object):
    '''Utility class for construction a sequence of factor convertions and
       rendering as `layout.convertion`
    '''

    error = None
    def __init__(self, num, seed=True):
        self.terms = []
        num = as_physnum(num)
        if seed:
            self.current_value = 1
            self.add_term(num)
        else:
            self.current_value = num
        self.name = num.name

    def factor_convert(self, num, den, power=1):
        num,den = map(self.x_as_physum, [num,den])
        self.add_term(num, den, power)
        if self.current_value.unit == num.unit:
            self.current_value.unit = num.unit #force to have same form when equivalent
        return self

    def prefix_convert(self, to_unit, power=1):
        from_unit = self.powered_unit(self.current_value.unit, power).cannonicalized()
        to_unit = self.powered_unit(to_unit, power).cannonicalized()
        assert from_unit.without_prefix() == to_unit.without_prefix(), \
               '%s to %s' % (from_unit, to_unit)
        dp = (to_unit / from_unit).cannonicalized().prefix.power
        quantity = self.current_value.quantity * Decimal(10) ** (-dp * power)
        unit = from_unit
        s = sign(dp)
        cnvs = [[1000, abs(dp)//3], [10**(abs(dp)%3), 1 if dp%3 else 0]]
        if from_unit.prefix.power % 3:
            cnvs.reverse()
        for factor,iterations in cnvs:
            for _ in xrange(iterations):
                new_unit = (unit * factor**s)
                num = (factor if s==-1 else 1, new_unit)
                den = (1 if s==-1 else factor, unit)
                self.add_term(num, den, power)
                unit = new_unit
        assert unit == to_unit
        self.current_value = PhysNum(quantity, to_unit ** power)
        return self

    def path_convert(self, to_unit, power=1):
        from_unit = self.powered_unit(self.current_value.unit, power)
        to_unit = self.powered_unit(to_unit, power)
        n_from_unit = convertion_graph.normalize_unit(from_unit)
        n_to_unit = convertion_graph.normalize_unit(to_unit)
        dp_from, node_from = convertion_graph.get_node(n_from_unit)
        dp_to, node_to = convertion_graph.get_node(n_to_unit)
        if dp_from!=0:
            self.prefix_convert(node_from.unit**power, power)
        dp_cnv, path = convertion_graph.find_best_convertion_path(n_from_unit, n_to_unit)
        assert dp_from - dp_to == dp_cnv
        assert path, 'not path from %s -> %s' % (from_unit, to_unit)
        unit = node_from.unit
        for arc in path.arcs:
            if not self.add_term((1 if arc.invert_factor else arc.factor, arc.node.unit),
                                 (arc.factor if arc.invert_factor else 1, unit),
                                 power):
                return self
            unit = arc.node.unit
        if dp_to!=0:
            self.prefix_convert(to_unit**power, power)
        return self

    @staticmethod
    def powered_unit(op, power):
        assert power != 0
        op = U.as_unit(op)
        if power == 1: #everything is a power of 1
            return op
        if isinstance(op, U.PrimitiveUnit):
            op = U.primunit_to_compound(op)
        op = op.ordered()
        if len(op.atoms_and_powers) != 1 or op.atoms_and_powers[0][1] != power:
            raise ValueError("%s is not a unit of power %d" % (op,power))
        base = op.atoms_and_powers[0][0]
        assert op.prefix.power % power == 0
        return U.CompoundUnit([[base, 1]], U.Prefix.from_power(op.prefix.power // power))

    @staticmethod
    def x_as_physum(op):
        if op is None:
            return None
        if isinstance(op, tuple):
            num,unit = op
            op = PhysNum(num, unit)
        op = as_physnum(op)
        if isinstance(op.quantity, (int,long)):
            op = 1*op
            op.quantity = Decimal(op.quantity)
        return op

    def add_term(self, num, den=None, power=1):
        if den is None:
            den = as_physnum(Decimal(1))
        num,den = map(self.x_as_physum, [num,den])
        self.terms.append([num, den, power])
        factor = num / den
        if power != 1:
            assert not isinstance(factor.quantity, (int,long))
            factor = factor ** power
        self.current_value = self.current_value * factor
        #print 'XXX', self.current_value
        self.current_value.unit = self.current_value.unit.cannonicalized()
        if self.current_value.quantity is dne:
            if not annotator.annotating:
                raise
            self.error = layout.zero_division_error
            return False
        return True

    def finish(self):
        result = 1 * self.current_value
        result.name = self.name
        if annotator.annotating:
            annotator.annotate(
                layout.equals(
                  self.make_convertion(),
                  V(result)))
            if self.error is not None:
                annotator.annotate(self.error)
        return result

    def make_convertion(self):
        cnv = []
        n_terms = len(self.terms)
        cnv = []
        for i_term,term in enumerate(self.terms):
            parts = []
            num,den,power = term
            for i_part,part in enumerate([num,den]):
                if part is None:
                    parts.append(None)
                    continue
                op = V(part, crossed_unit=not (i_term==len(self.terms)-1 and i_part==0))
                if power != 1:
                    op = layout.power(layout.group(op), power)
                #if i_term==0 and i_part==0 and self.name is not None:
                #    op = layout.make_span(op, layout.compound_name(unicode(self.name)))
                parts.append(op)
            cnv.append(parts)
        return layout.make_convertion(*cnv)


# # # # # # # # # # # # #
# Path-Based Convertion #
# # # # # # # # # # # # #

# Defines a graph where the nodes are units and the arcs describe pair-wise convertion
# between two different units. To convert from a source to target destination, find the
# optimal path through this graph connecting the source and target nodes. Then
# perform the combined convertion given by the path.

class NoSuchConvertionError(LookupError):

    def __init__(self, unit_from, unit_to):
        self.unit_from = unit_from
        self.unit_to = unit_to
    def __str__(self):
        return 'no convertion from %s to %s' % (self.unit_from, self.unit_to)

class ConvertionNode(object):

    def __init__(self, unit):
        self.unit = unit
        self.convertion_arcs = []

class ConvertionArc(object):

    def __init__(self, node, factor, invert_factor=False, weight=1):
        self.node = node
        self.factor = factor
        self.invert_factor = invert_factor
        self.weight = weight

class ConvertionPath(object):

    def __init__(self, arcs):
        self.arcs = arcs

    def calculate_weight(self):
        return reduce(operator.mul, (arc.weight for arc in self.arcs), 1)

    def calculate_factor(self):
        #minimize divisions
        num,den = (reduce(operator.mul,
                       (arc.factor for arc in self.arcs if arc.invert_factor==invert), 1)
                   for invert in [False,True])
        num,den = (Decimal(op) if isinstance(op, (int,long)) else op
                   for op in [num,den])
        factor = num / den
        assert typep(factor, lossless_number_type)
        return factor

class ConvertionGraph(object):

    def __init__(self):
        self.unit_nodes = {}

    @staticmethod
    def normalize_unit(unit):
        if isinstance(unit, U.PrimitiveUnit):
            unit = U.primunit_to_compound(unit)
        return unit.ordered()

    def find_best_convertion_path(self, unit_from, unit_to):
        power_delta, paths = self.find_convertion_paths(unit_from, unit_to)
        if not paths:
            return power_delta, None
        paths.sort(key=lambda path: path.calculate_weight())
        return power_delta, paths[-1]

    def register(self, unit_from, unit_to, factor, weight=None):
        unit_from = U.as_unit(unit_from)
        unit_to = U.as_unit(unit_to)
        if isinstance(factor, str):
            factor = Decimal(factor)
        if not typep(factor, lossless_number_type):
            raise TypeError("invalid factor type %r" % (factor,))
        if weight is None:
            weight = self.calculate_factor_weight(factor)
        dpower_from, node_from = self.get_node(self.normalize_unit(unit_from))
        dpower_to, node_to = self.get_node(self.normalize_unit(unit_to))
        if node_from is node_to:
            raise ValueError("register prefix convertion; done automatically")
        factor = factor * Decimal(10) ** (dpower_to - dpower_from)
        node_from.convertion_arcs.append(ConvertionArc(node_to, factor,
                                                       invert_factor=False,
                                                       weight=weight))
        node_to.convertion_arcs.append(ConvertionArc(node_from, factor,
                                                     invert_factor=True,
                                                     weight=weight))
    @staticmethod
    def calculate_factor_weight(op):
        if isinstance(op, (int,long)):
            return 100
        if isinstance(op, Decimal):
            return 10
        return op.sigfigs

    def get_node(self, unit):
        prefix,name,abbrev = unit.get_name_abbrev_prefix()
        if name is not None:
            key = abbrev
            node_unit = unit * 10 ** - prefix.power
        else:
            key = unit.without_prefix()
            prefix = unit.prefix / key.prefix
            node_unit = key
        power_delta = prefix.power
        try:
            node = self.unit_nodes[key]
        except KeyError:
            node = self.unit_nodes[key] = ConvertionNode(node_unit)
        return power_delta, node

    def find_convertion_paths(self, unit_from, unit_to):
        dpower_from, node_from = self.get_node(unit_from)
        dpower_to, node_to = self.get_node(unit_to)
        return [dpower_from - dpower_to,
                list(ConvertionPath(arcs) for arcs in
                     self.iter_paths_between(node_from, node_to, (), set()))]

    def iter_paths_between(self, node, end_node, arc_path, seen_nodes):
        if node in seen_nodes:
            return
        seen_nodes.add(node)
        if node is end_node:
            yield arc_path
        for arc in node.convertion_arcs:
            for arc_path_x in self.iter_paths_between(arc.node, end_node, arc_path + (arc,), seen_nodes):
                yield arc_path_x

    def calculate_convertion_factor_dimensionally(self, unit_from, unit_to):
        dimensionally_from = unit_from.get_dimensionality().cannonicalized()
        dimensionally_to = unit_from.get_dimensionality().cannonicalized()
        if dimensionally_from != dimensionally_to:
            raise ValueError("cannot convert %s(%s) to %s(%s); different dimensionallities"
                             % (unit_from, dimensionally_from,
                                unit_to, dimensionally_to))
        #
        #   could implement this later
        #
        raise NoSuchConvertionError(unit_from, unit_to)

def convertion_graph():
    def parse_op(op):
        op = op.strip()
        m = re.match(r''' \( ( [\deE.]+ \s+ [a-zA-Z_] +) \) \*\* (\d+) ''',
                     op, re.VERBOSE)
        if m:
            factor_unit,power = m.groups()
            return parse_op(factor_unit) ** int(power)
        factor,unit = op.split()
        factor = (Decimal if re.search('[.eE]', factor) else int)(factor)
        unit = U.parse_unit(unit)
        return PhysNum(factor, unit)
    converter = ConvertionGraph()
    factors = '''
    # # # # # # # # # # # #
    # Convertion Factors  #
    # # # # # # # # # # # #

    #lengths
    1 in = 2.54 cm
    1 ft = 12 in
    1 yd = 3 ft
    1 furlong = 220 yd
    1 mile = 8 furlong
    #mass
    1 oz = 28.349523 g
    1 lb = 16 oz
    1 st = 14 lb
    1 ton = 2240 lb
    #volume
    1 mL = 1 cc
    1 yd3 = (3 ft)**3
    1 ft3 = (12 in)**3
    1 in3 = (2.54 cm)**3
    1 floz = 28.4130625 mL
    1 floz = 1.7339 in3
    1 pt = 20 floz
    1 qt = 2 pt
    1 gal = 4 qt
    #times
    1 min = 60 s
    1 hour = 60 min
    1 day = 24 hours
    1 week = 7 days
    1 year = 365.25 days
    #pressure
    1 bar = 100000 Pa
    1 atm = 101325 Pa
    1 torr = 133.322 Pa
    1 torr = 1 mmHg
    1 psi = 6.394e3 Pa
    #quantities
    1 mol = 6.02e23 quantity


    '''
    for line in factors.split('\n'):
        line = line.split('#')[0].strip()
        if not line:
            continue
        lpn,rpn = map(parse_op, line.split('='))
        converter.register(lpn.unit, rpn.unit, rpn.quantity / lpn.quantity
                           if lpn.quantity != 1 else rpn.quantity)
    return converter
convertion_graph = convertion_graph()
