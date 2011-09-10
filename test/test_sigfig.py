
from decimal import Decimal

from physmath.sigfig import SigFig

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
    for line in old_check_data.strip().split('\n'):
        line = line.split('#',1)[0].strip()
        if not line:
            continue
        inp,sigfigs,lsp,rep = line.split()
        sigfigs,lsp = map(int, [sigfigs, lsp])
        yield check_old, inp, sigfigs, lsp, rep

def check_old(inp, sigfigs, lsp, rep):
    s = SigFig(inp)
    assert s.sigfigs == sigfigs
    assert s.least_significant_place == lsp
    assert str(s) == rep
