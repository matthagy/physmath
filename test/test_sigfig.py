
import operator
from decimal import Decimal

from physmath.sigfig import SigFig

def parse_lines(data, n=None, split=lambda line: line.split()):
    for line in data.split('\n'):
        line = line.split('#',1)[0]
        line = line.strip()
        if line:
            parts = split(line)
            if n is not None and len(parts) != n:
                raise ValueError("bad line %r; split into %d parts when expecting %d" %
                                 (line,len(parts),n))
            yield parts

def create_seq_checks(name, func, seq):
    def wrap():
        for args in seq:
            yield (func,) + tuple(args)
    wrap.func_name = name
    globals()[name] = wrap


def check_from_string(expr, canonical_str, eq):
    sf = SigFig(expr)
    assert str(sf) == canonical_str

    assert sf == eq
    assert sf <= eq
    assert sf >= eq

    assert eq == sf
    assert eq <= sf
    assert eq >= sf

string_checks = [
    [      '1',       '1',         1],
    [      '1',       '1',         1],
    [     '10',      '10',         10],
    [     '10',      '10',         10],
    [     '10',      '10',         10],
    [     '-3',      '-3',        -3],
    [    '1.0',     '1.0',         1],
    [    '1.4',     '1.4',         Decimal('1.4')],
    [   '-2.7',    '-2.7',         Decimal('-2.7')],
    [    '-18',     '-18',         -18],
    [    '120',     '120',         120],
    [    '73.30',  '73.30',        Decimal('73.30')],
]

create_seq_checks('test_from_string', check_from_string, string_checks)


def check_bounded(expr, lt, gt):
    sf = SigFig(expr)
    assert lt < sf
    assert lt <= sf
    assert not (lt >= sf)

    assert sf > lt
    assert sf >= lt
    assert not (sf <= lt)

    assert gt > sf
    assert gt >= sf
    assert not (gt <= sf)

    assert sf < gt
    assert sf <= gt
    assert not (sf >= gt)

bound_checks = [
    ['1',     0,       10],
    ['5',     0,       20],
    ['20',    5,       40],
]

create_seq_checks('test_bounded', check_bounded, bound_checks)


def test_comparisions():
    assert SigFig('10') > 0
    assert SigFig('10') > 3

    assert not (SigFig('10') > 7)
    assert SigFig('10.') > 7


def check_n_sigfigs(expr, n):
    assert expr
    assert not expr.startswith('-')
    assert SigFig(expr).sigfigs == n
    assert SigFig('-' + expr).sigfigs == n

n_sigfigs_checks = [
    [       '1', 1],
    [      '10', 1],
    [     '10.', 2],
    [      '23', 2],
    [     '230', 2],
    [    '7400', 2],
    [    '6000', 1],
    [     '532', 3],
    [     '200', 1],
    [     '850', 2],
    [     '45.', 2],
    [    '45.0', 3],
    [   '45.00', 4],
    [  '45.005', 5],
    [   '45.05', 4],
    [  '45.050', 5],
    [    '2030', 3],
    [   '20301', 5],
    [   '20300', 3],
    [  '20300.', 5],
    [ '20300.0', 6],
    [     '1e0', 1],
    [   '1.0e0', 2],
    [     '2e7', 1],
    [   '2.5e3', 2],
    [ '8.700e8', 4],
    [  '3.23e4', 3],
    [ '2.00e-6', 3],
]

create_seq_checks('test_n_sigfigs', check_n_sigfigs, n_sigfigs_checks)


old_check_data = '''
#input text      sigfigs     least significant place      text representation
1                1           0                            1
6                1           0                            6
6.               1           0                            6
6.0              2           -1                           6.0
6.07             3           -2                           6.07
0.07             1           -2                           0.07
.07              1           -2                           0.07
7e-2             1           -2                           0.07
0.034            2           -3                           0.034
0.0340           3           -4                           0.0340
3.4e-2           2           -3                           0.034
3.40e-2          3           -4                           0.0340
3.4e-4           2           -5                           3.4e-4
3.40e-4          3           -6                           3.40e-4
10               1           1                            10
1e1              1           1                            10
1.0e1            2           0                            1.0e1
10.              2           0                            1.0e1
200              1           2                            200
200.             3           0                            2.00e2
2.00e2           3           0                            2.00e2
2.0e2            2           1                            2.0e2
0                1           0                            0
0.0              1           -1                           0.0
185              3           0                            185
'''

def test_old():
    for inp,sigfigs,lsp,rep in parse_lines(old_check_data, n=4):
        sigfigs,lsp = map(int, [sigfigs, lsp])
        yield check_old, inp, sigfigs, lsp, rep
        if any(d in inp for d in '123456789'):
            yield check_old, '-'+inp, sigfigs, lsp, '-'+rep

def check_old(inp, sigfigs, lsp, rep):
    s = SigFig(inp)
    assert s.sigfigs == sigfigs
    assert s.least_significant_place == lsp
    assert str(s) == rep, '%s != %s' % (str(s), rep)

def parse_number(s):
    if s.endswith('x'):
        return Decimal(s[:-1:])
    return SigFig(s)

def check_unop(operand, op, expected_answer):
    actual_answer = {'+' : operator.pos,
                     '-' : operator.neg}[op](parse_number(operand))
    expected_answer = parse_number(expected_answer)
    assert actual_answer == expected_answer, '%s != %s' % (actual_answer, expected_answer)

unop_checks = '''
    - 2.34 = -2.34
    + 2.34 = 2.34
    - -1.087 = 1.087
    + -1.087 = -1.087
    + 1 = 1
    - 1 = -1
    - -1 = 1
    + 300 = 300
    - 300 = -300
    - -300 = 300
    + -300 = -300
'''

def test_unops():
    for op, operand, eq, expected_answer in parse_lines(unop_checks, n=4):
        assert eq == '='
        assert op in '+-'
        yield check_unop, operand, op, expected_answer

def check_binop(lop, op, rop, expected_answer):
    actual_answer = {'+'  : operator.add,
                     '-'  : operator.sub,
                     '/'  : operator.div,
                     '*'  : operator.mul,
                     '**' : operator.pow}[op](parse_number(lop),
                                             parse_number(rop))
    expected_answer = parse_number(expected_answer)
    assert actual_answer == expected_answer, '%s != %s' % (actual_answer, expected_answer)

binop_checks = '''
    1 + 1 = 2
    1 + 0 = 1
    1 + 10 = 10
    1 + 10. = 11
    1.0 + 10.0 = 11.0
    1.0 + 10. = 11
    5 + 6 = 11
    2 + 3 = 5
    5 + 5 = 1.0e1
    12 + 8 = 2.0e1

    20 + 30 = 50
    20 + 32 = 50
    20 + 34 = 50
    20 + 35 = 60
    20 + 36 = 60
    20 + 45 = 60
    21 + 45 = 66

    18 + 6 = 24
    18 + 1.0 = 19
    18 + 1.3 = 19
    18 + 1.5 = 2.0e1
    18 + 1.7 = 2.0e1

    28 + 1.4 = 29
    28 + 1.5 = 3.0e1
    28 + 1.6 = 3.0e1

    10 - 2 = 10
    10. - 2 = 8
    1.0e1 - 2 = 8
    10 - 5 = 5

    10 + -2 = 10
    10. + -2 = 8
    1.0e1 + -2 = 8
    10 + -5 = 5

    103 + 7 = 1.10e2
    252 + 8 = 2.60e2
    130 + 10 = 140
    130. + 10. = 1.40e2

    20 + -10 = 10
    20 - 10 = 10

    1.52 + 2.48 = 4.00
    2.32 + 1.43 = 3.75
    2.80 + 0.10 = 2.90

    5.3x + 0.05x = 5.35x
    5.3x + 0.05x = 5.35

    5.3 + 0.05 = 5.4
    5.3 + 0.04 = 5.3
    5.3 + 0.06 = 5.4

    5.3 + 0.05x = 5.4
    5.3 + 0.04x = 5.3
    5.3 + 0.06x = 5.4

    5.3x + 0.05 = 5.4
    5.3x + 0.04 = 5.3
    5.3x + 0.06 = 5.4

    5.2x + 0.05x = 5.25x
    5.2x + 0.05x = 5.25

    5.2 + 0.05 = 5.2
    5.2 + 0.04 = 5.2
    5.2 + 0.06 = 5.3

    5.2 + 0.05x = 5.2
    5.2 + 0.04x = 5.2
    5.2 + 0.06x = 5.3

    5.2x + 0.05 = 5.2
    5.2x + 0.04 = 5.2
    5.2x + 0.06 = 5.3

    3.4x + -0.2x = 3.2x
    3.4 + -0.2 = 3.2
    3.4 - 0.2 = 3.2

    3.4x + -0.25x = 3.15x

    3.4 + -0.25 = 3.2
    3.4 + -0.24 = 3.2
    3.4 + -0.26 = 3.1

    3.4 + -0.25x = 3.2
    3.4 + -0.24x = 3.2
    3.4 + -0.26x = 3.1

    3.4x + -0.25 = 3.2
    3.4x + -0.24 = 3.2
    3.4x + -0.26 = 3.1

    3.4 - 0.25 = 3.2
    3.4 - 0.24 = 3.2
    3.4 - 0.26 = 3.1

    3.4 - 0.25x = 3.2
    3.4 - 0.24x = 3.2
    3.4 - 0.26x = 3.1

    3.4x - 0.25 = 3.2
    3.4x - 0.24 = 3.2
    3.4x - 0.26 = 3.1


    1 * 1 = 1
    1 * 2 = 2
    1 * 10 = 10

    1 * 1.0e1 = 1.0e1
    1 * 10. = 1.0e1

    1 * -20 = -20
    1 * -2.0e1 = -2.0e1

    2x * 2x = 4x
    2 * 2 = 4
    2.0 * 2.0 = 4.0
    2 * 2.0 = 4
    2x * 2 = 4
    2x * 2.0 = 4.0
    2x * 2.00 = 4.00

    3x * -3x = -9x
    3 * -3 = -9
    3.0 * -3.0 = -9.0
    3 * -3.0 = -9
    3x * -3 = -9
    3x * -3.0 = -9.0
    3x * -3.00 = -9.00
'''

def test_binops():
    for l,o,r,e,a in parse_lines(binop_checks):
        assert e == '='
        assert o in ['+', '-', '/', '*', '**']
        yield check_binop, l, o, r, a


