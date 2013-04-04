'''Collections of various calculations used in chemistry/physics
'''

from hlab.sigfigs import SigFig
from hlab import xunits as U

from jamenson.runtime.atypes import anytype, Seq, as_optimized_type, typep
from jamenson.runtime.multimethod import defmethod
from jamenson.runtime.struct import defstruct

import .mathlayout as ML
from .mathlayout import V
from .convert import convert_temperature, convert_gas_volume, convert_unit_prefix, convert_by_path


# utilities to use generators to construct calculations

def gen_expand(gen):
    itr = iter(gen)
    acc = []
    no_send = object()
    send = no_send
    while True:
        try:
            value = itr.next() if send is no_send else itr.send(send)
        except StopIteration:
            break
        if isinstance(value, generator_type):
            sub_acc = gen_expand(value)
            send = sub_acc.pop(-1)
            acc.extend(sub_acc)
            if send is None:
                acc.append(None)
                break
        else:
            send = None
            acc.append(value)
    return acc

def equation_set_iterator(func):
    def wrap(*args, **kwds):
        mls = gen_expand(func(*args, **kwds))
        text = mls.pop(0)
        result = mls.pop(-1)
        mls = list(x for x in mls if x is not None)
        if mls:
            for arg in args:
                if hasattr(arg, 'add_calculation'):
                    break
            else:
                raise 'hell'
            arg.add_calculation(ML.equation_set(text, mls))
            for ml in mls:
                if isinstance(ml, ML.error):
                    arg.flag_error()
        return result
    return wrap

def x_convert_unit_prefix(num, to_unit, name=None):
    cnv,answer = convert_unit_prefix(num, to_unit, name)
    yield cnv
    yield answer

def x_convert_by_path(num, to_unit, name=None):
    cnv,answer = convert_by_path(num, to_unit, name)
    yield cnv
    yield answer

def x_convert_temperature(temperature, unit):
    cnv,answer = convert_temperature(temperature, unit)
    yield cnv
    yield answer

def x_convert_gas_volume(num, to_unit, name=None):
    cnv,answer = convert_gas_volume(num, to_unit, name)
    yield cnv
    yield answer


# converters

@equation_set_iterator
def convert_mols_to_mols(calculator, mols, to_unit=None, name=None):
    if to_unit is None:
        to_unit = U.quantities.mol
    yield simple_convertion_title(name, mols.unit, to_unit)
    mols = yield x_convert_unit_prefix(mols, to_unit, name)
    yield mols

@equation_set_iterator
def convert_mass_to_mols(calculator, mw, mass, to_unit=None, name=None):
    assert mw.unit == U.mws.g_mol
    if to_unit is None:
        to_unit = U.quantities.mol
    yield simple_convertion_title(name, mass.unit, to_unit)
    mass = yield x_convert_by_path(mass, U.masses.g, name)
    mol = yield simple_factor_convert(mass, mw, name=name, invert=True)
    mol = yield x_convert_unit_prefix(mol, to_unit, name)
    yield mol

@equation_set_iterator
def convert_mols_to_mass(calculator, mw, mols, to_unit=None, name=None):
    assert mw.unit == U.mws.g_mol
    assert mols.unit.without_prefix() == U.quantities.mol
    if to_unit is None:
        to_unit = U.masses.g
    yield simple_convertion_title(name, mw.unit, to_unit)
    mols = yield x_convert_unit_prefix(mols, U.quantities.mol, name)
    mass = yield simple_factor_convert(mols, mw, name=name)
    mass = yield x_convert_unit_prefix(mass, to_unit, name)
    yield mass

@equation_set_iterator
def convert_liquid_volume_to_mols(calculator, concentration, volume, to_unit=None, name=None):
    assert concentration.unit.without_prefix() == U.concentrations.mol_L
    assert volume.unit.get_dimensionality() == U.dimensionalities.volume
    if to_unit is None:
        to_unit = U.quantities.mol
    yield simple_convertion_title(name, volume.unit, to_unit)
    concentration = yield x_convert_unit_prefix(concentration, U.concentrations.mol_L)
    concentration = concentration * 1
    concentration.unit = concentration.unit.cannonicalized()
    volume = yield x_convert_by_path(volume, U.liquid_volumes.L)
    mol = yield simple_factor_convert(volume, concentration, name)
    mol = yield x_convert_unit_prefix(mol, to_unit, name)
    yield mol

@equation_set_iterator
def convert_mols_to_liquid_concentration(calculator, volume, mols, to_unit=None, name=None):
    assert volume.unit.get_dimensionality() == U.dimensionalities.volume
    assert mols.unit.without_prefix() == U.quantities.mol
    if to_unit is None:
        to_unit = U.concentrations.M
    yield simple_convertion_title(name, mols.unit, to_unit)
    mols = yield x_convert_unit_prefix(mols, U.quantities.mol, name)
    volume = yield x_convert_by_path(volume, U.liquid_volumes.L)
    concentration = yield simple_factor_convert(mols, volume, name, invert=True)
    concentration = yield x_convert_unit_prefix(concentration, to_unit, name)
    yield concentration

R = U.PhysicalNumber(Decimal('8.314'),
                     U.pressures.Pa * U.gas_volumes.m3 / U.temperatures.K / U.quantities.mol)

@equation_set_iterator
def convert_gas_PVT_to_mols(calculator, pressure, volume, temperature, to_unit=None, name=None):
    assert pressure.unit.get_dimensionality() == U.dimensionalities.pressure
    assert volume.unit.get_dimensionality() == U.dimensionalities.volume
    assert temperature.unit.get_dimensionality() == U.dimensionalities.temperature
    if to_unit is None:
        to_unit = U.quantities.mol
    yield simple_convertion_title(name, volume.unit, to_unit, u'using the ideal gas law')
    volume = yield x_convert_gas_volume(volume, U.gas_volumes.m3, name)
    pressure = yield x_convert_by_path(pressure, U.pressures.Pa)
    temperature = yield x_convert_temperature(temperature, U.temperatures.K)
    try:
        mols = pressure * volume / (R * temperature)
    except ZeroDivisionError:
        error = ML.error(u'division by zero leads to a number that does not exist (D.N.E.)')
        mols = None
        mols_ml = V(ML.dne(), mols_unit, name=name)
    else:
        error = None
        mols.unit = mols.unit.cannonicalized()
        assert mols.unit == U.quantities.mol
        mols_ml = V(mols, name=name)
    cnv = ML.make_convertion(
        [V(volume, name=name, crossed_unit=True), V(1)],
        [V(pressure, crossed_unit=True), V(1)],
        [V(1), V(temperature, crossed_unit=True)],
        [ML.x_physical_number(1, [ML.crossed_unit(U.temperatures.K), U.quantities.mol]),
         ML.x_physical_number(R.quantity, [ML.crossed_unit(U.pressures.Pa),
                                           ML.crossed_unit(U.gas_volumes.m3)])])
    yield ML.equals(cnv, mols_ml)
    if error:
        yield error
    else:
        mols = yield x_convert_unit_prefix(mol, to_unit, name)
    yield mols

@equation_set_iterator
def convert_gas_PnT_to_volume(calculator, pressure, mols, temperature, to_unit=None, name=None):
    assert pressure.unit.get_dimensionality() == U.dimensionalities.pressure
    assert mols.unit.without_prefix() == U.quantities.mol
    assert temperature.unit.get_dimensionality() == U.dimensionalities.temperature
    if to_unit is None:
        to_unit = U.gas_volumes.L
    yield simple_convertion_title(name, mols.unit, to_unit, u'using the ideal gas law')
    mols = yield x_convert_p(mols, U.quantities.mol, name)
    pressure = yield x_convert_by_path(pressure, U.pressure.Pa)
    temperature = yield x_convert_temperature(temperature, U.temperatures.K)
    try:
        volume = mols * temperature * R / pressure
    except ZeroDivisionError:
        error = ML.error(u'division by zero leads to a number that does not exist (D.N.E.)')
        volume = None
        volume_ml = V(ML.dne(), mols_unit, name=name)
    else:
        error = None
        volume.unit = volume.unit.cannonicalized()
        volume_ml = V(volume, name=name)
    raise 'not yet finished'



# convertion utilities

def simple_factor_convert(num, factor, name=None, cannonicalize=True, invert=False):
    answer_unit = (num.unit / factor.unit) if invert else (num.unit * factor.unit)
    if cannonicalize:
        answer_unit = answer_unit.cannonicalized()
    n,d = factor.split_posneg()
    if invert:
        n,d = d,n
    try:
        answer = num*n/d
    except ZeroDivisionError:
        error = ML.error(u'division by zero leads to a number that does not exist (D.N.E.)')
        answer = None
        answer_ml = V(ML.dne(), answer_unit, name=name)
    else:
        error = None
        answer.unit = answer_unit
        answer_ml = V(answer, name=name)
    #need to fix unit crossing
    cnv = ML.make_convertion([V(num, crossed_unit=True, name=name),
                              None],
                             ([V(n),
                               V(d, crossed_unit=True)]))
    yield ML.equals(cnv, V(answer, name=name))
    if error:
        yield error
    yield answer

def simple_convertion_title(name, from_unit, to_unit, extra=None):
    return u'Convert%s from %s to %s%s.' % (name_fix(name),
                                          plural_unit_name(from_unit),
                                          plural_unit_name(to_unit),
                                            u' '+extra if extra else u'')
def name_fix(name):
    if not name:
        return u''
    return u' ' + unicode(name)

special_unit_names = {
    U.gas_volumes.m3: u'cubic meters',
    U.metric.cm**3: u'cubic centimeters',
    U.imperial.ft**3: u'cubic feet',
    U.imperial.in_**3: u'cubic inches',
    U.imperial.yd**3: u'cubic yards'
}

def plural_unit_name(unit):
    try:
        return special_unit_names[unit]
    except KeyError:
        pass
    name = unicode(unit.get_name())
    if name.endswith(u'y'):
        return name[:-1:] + u'ies'
    return name + u's'


class Calculator(ML.BaseML):

    def __init__(self):
        self.calculations = []
        self.error = False

    def add_calculation(self, calc):
        if calc is not None:
            self.calculations.append(ML.as_ml(calc))

    def flag_error(self):
        self.error = True

@defmethod(ML.as_ml, [Calculator])
def meth(c):
    return c

@defmethod(ML.get_ml_json, [Calculator])
def meth(c):
    return ML.get_ml_json(ML.calculations(mls=c.calculations))

def solve_expression(*terms):
    pass


mols = solve_expression([pressure,    'Pa',     1],
                        [volume,      'm3',     1],
                        [None,        to_unit, -1]
                        [temperature, 'K',     -1],
                        [R,           None,    -1])



if 0:
    calc = Calculator()
    conc = U.PhysicalNumber(SigFig('1.07'), U.concentrations.mmol_L)
    vol = U.PhysicalNumber(SigFig('250.3'), U.liquid_volumes.gal)
    mols = convert_liquid_volume_to_mols(calc, conc, vol, 1e-6 * U.quantities.mol, name=u'agent x')
    print conc
    print vol
    print mols
    show_math(calc)

ppn = U.parse_physical_number
pu = U.parse_unit

if 1:
    calc = Calculator()
    volume = ppn('34.05 gallons')
    mols = ppn('1.30e8 nmols')
    concentration = convert_mols_to_liquid_concentration(calc, volume, mols, pu('mM'), name=u'stuff')
    show_math(calc)

if 0:
    calc = Calculator()
    pressure = ppn('1.007 atm')
    volume = ppn('3.54e5 yd3')
    temperature = ppn('23.4 C')
    mols = convert_gas_PVT_to_mols(calc, pressure, volume, temperature, U.quantities.mol, name=u'stuff')
    show_math(calc)


#from pprint import pprint
#pprint(ML.get_ml_json(calc))
