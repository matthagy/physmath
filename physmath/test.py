
from reload import fixup
fixup(0)

from physmath.units import parse_unit as pu
from physmath.physnum import parse_physical_number as ppn
from physmath.convert import convert, convert_factor
from physmath.annotator import annotator
from physmath import layout

import browsemath.client
reload(browsemath.client)
from browsemath.client import show_math

if 0:
    label = annotator.push('stuff')
    T1 = ppn('34.54s C')
    T2 = convert(T1, pu('K'))
    print T1, T2
    acc = annotator.pop(label)
    print acc
    es = layout.equation_set(u'Temperature convertion', acc)

    show_math(es)

if 0:
    label = annotator.push('stuff')
    n1 = ppn('34.54s mmol C')
    n2 = convert(n1, pu('mol'))
    print n1, n2
    acc = annotator.pop(label)
    es = layout.equation_set(u'Molar convertion', acc)
    es = layout.as_ml(es)
    show_math(es)

if 1:
    label = annotator.push('stuff')
    n1 = ppn('1.03e5s yd3 deadly gas')
    n2 = convert(n1, pu('gal'))
    print n1, n2
    acc = annotator.pop(label)
    es = layout.equation_set(u'Volume convertion', acc)
    es = layout.as_ml(es)
    show_math(es)


if 0:
    label = annotator.push('stuff')
    n1 = ppn('1.03e5s mi run')
    n2 = convert(n1, pu('nm'))
    print n1, n2
    acc = annotator.pop(label)
    es = layout.equation_set(u'Length convertion', acc)
    es = layout.as_ml(es)
    show_math(es)

if 0:
    label = annotator.push('stuff')
    n1 = ppn('1.03s mM run')
    n2 = convert(n1, pu('nmol/L'))
    print n1, n2
    acc = annotator.pop(label)
    es = layout.equation_set(u'Concentration convertion', acc)
    es = layout.as_ml(es)
    show_math(es)

if 0:
    r = convert_factor(ppn('10.3s g carbon'), ppn('1 mol'), ppn('12.01d g'))
    acc = annotator.pop(label)
    es = layout.equation_set(u'Concentration convertion', acc)
    es = layout.as_ml(es)
    show_math(es)
    print r
