
from __future__ import absolute_import

from collections import defaultdict
from hlab.lexing import Lexer, LexicalError
from hlab.bases import AutoRepr

from .sigfig import test_sigfigs

class SigFig(AutoRepr):

    def __init__(self, arg):
        if isinstance(arg, unicode):
            arg = str(arg)
        if isinstance(arg, str):
            arg = parse_string(arg)
        elif isinstance(arg, (int,long)):
            arg = convert_integer(arg)
        elif isinstance(arg, SigFig):
            arg = (arg.sign, arg.digits, arg.power)
        (self.sign,
         self.digits,
         self.power
         ) = arg

    def repr_args(self):
        return [str(self)]

    #legacy method
    def as_scientific(self):
        return self

    def as_tuple(self):
        return self.sign, self.digits, self.power

    @property
    def sigfigs(self):
        if self.digits[0] != 0:
            return len(self.digits)
        return max(1, len(self.digits) - 1)

    @property
    def most_significant_place(self):
        return self.power

    @property
    def least_significant_place(self):
        return 1 + self.power - len(self.digits)

    def __str__(self):
        base, exp = self.get_format_args()
        if exp is None:
            return base
        return '%se%s' % (base, exp)

    def get_format_args(self, min_exp_power=3):
        if self.digits[0] == 0:
            return self._get_zero_format_args()

        sigfigs = self.sigfigs
        power = self.power

        if min_exp_power!=None and abs(power) > max(sigfigs, min_exp_power):
            return self._get_exp_format_args()

        digs = map(str, self.digits)

        #deal with trailing zeros
        if sigfigs > power:
            insignificant_zeros = 0
            #prevet trailing decimal, i.e. 10.
            if power > 0 and len(digs) == power+1 and digs[-1] == '0':
                return self._get_exp_format_args()
        else:
            #check for significant trailing zeros
            if self.digits[-1::] == (0,):
                return self._get_exp_format_args()
            #add insignificant trailing zeros
            insignificant_zeros = (1 + power - sigfigs)
        digs += ['0'] * insignificant_zeros

        #special case to prevent decimal points in zero
        if digs == ['0']:
            return '0', None

        #place deicmal point
        if power >= 0:
            if not insignificant_zeros and power+1 < len(digs):
                digs.insert(1+power, '.')
        else:
            digs = ['0.'] + ['0'] * (-1 - power) + digs

        return ['%s%s' % ('-' if self.sign else '', ''.join(digs)),
                None]

    def _get_zero_format_args(self):
        assert set(self.digits) == set([0])
        if self.power > 0:
            return self._get_exp_format_args()
        assert self.digits == (0,)
        if self.power == 0:
            return '0', None
        return '0.' + '0' * -self.power, None

    def _get_exp_format_args(self):
        return ['%s%d%s' % ('-' if self.sign else '',
                            self.digits[0],
                            '.' + ''.join(map(str, self.digits[1::])) if
                            len(self.digits) > 1 else ''),
                self.power]

    def __neg__(self):
        return self.__class__((1 if self.sign==0 else 1,
                               self.digits, self.power))

    def __sub__(self, other):
        return self + (-other)

    def __add__(self, other):
        assert isinstance(other, SigFig)

        li,lf,ld = self.create_arithmetic_array()
        ri,rf,rd = other.create_arithmetic_array()

        si = 1 + max(li, ri)
        sf = max(lf, rf) - 2
        sd = [0] * (1 + si - sf)

        for oi,of,od in [[li,lf,ld], [ri,rf,rd]]:
            for i,d in zip(xrange(oi, of-1, -1), od):
                if i>=sf:
                    sd[-(i-si)] += d

        for i in xrange(len(sd)-1, 0, -1):
            if sd[i] > 9:
                sd[i] -= 10
                sd[i-1] += 1
            if sd[i] < -9:
                sd[i] += 10
                sd[i-1] -= 1

        while len(sd)>3 and sd[0] == 0:
            sd.pop(0)

        sign = sd[0] < 0
        if sign:
            sd = [d*-1 for d in sd]

        for i in xrange(len(sd)-1, 0, -1):
            if sd[i] < 0:
                sd[i] += 10
                sd[i-1] -= 1

        if sd[-2] > 5 or (sd[-2] == 5 and sd[-1]%2):
            sd[-3] += 1

        assert len(sd) >= 3
        sd.pop(-1)
        sd.pop(-1)

        return self.__class__((sign, sd, sf + 2 + len(sd) - 1))

    def __mul__(self, other):
        assert isinstance(other, SigFig)

        li,lf,ld = self.create_arithmetic_array()
        ri,rf,rd = other.create_arithmetic_array()

        #print li,lf,ld
        #print ri,rf,rd
        #print

        acc = defaultdict(int)
        def dump():
            keys = sorted(acc, reverse=True)
            print keys
            print map(acc.get, keys)
            print

        for ai,ad in zip(xrange(li, lf-1, -1), ld):
            for bi,bd in zip(xrange(ri, rf-1, -1), rd):
                acc[ai+bi] += ad*bd

        #dump()

        for i in sorted(set(acc) | set(i+1 for i in acc)):
            r, acc[i] = divmod(acc[i], 10)
            acc[i+1] += r

        #dump()

        while acc and acc[max(acc)] == 0:
            del acc[max(acc)]

        #dump()

        sign = acc[max(acc)] < 0
        if sign:
            for i in acc:
                acc[i] *= -1

        for i in sorted(acc):
            if acc[i] < 0:
                acc[i] += 10
                acc[i-1] -= 1

        sigfigs = min(self.sigfigs, other.sigfigs)
        lsp = max(acc) - sigfigs + 1

        if acc[lsp-1] > 5 or (acc[lsp-1] == 5 and acc[lsp-2]%2):
            acc[lsp] += 1
        for i in list(acc):
            if i<lsp:
                del acc[i]

        #dump()

        digits = [acc[i] for i in xrange(max(acc), min(acc)-1, -1)]
        #print digits

        return self.__class__((sign, digits, max(acc)))

    def create_arithmetic_array(self):
        return [self.most_significant_place, self.least_significant_place,
                [d*(-1 if self.sign else 1) for d in self.digits]]

    def __truediv__(self, other):
        assert isinstance(other, SigFig)


        li,lp = self.create_integer_representation(10)
        ri,rp = other.create_integer_representation(0)

        ai = li // ri

        acc = {}
        for i in xrange(1<<22):
            if 10**i > ai:
                break
            n = (ai if i==0 else ai//10**(i-1)) % 10
            if n:
                acc[i] = n

    __div__ = __truediv__

    def create_integer_representation(self, extra):
        shift = len(self.digits) - 1
        power = self.power - shift
        acc = 0
        for d in self.digits[::-1]:
            acc *= 10
            acc += d
        acc *= 10 ** extra
        power -= extra
        acc *= (-1 if self.sign else 1)
        return acc,power



class SigFigBuilder(object):

    def __init__(self):
        self.digit_powers =  defaultdict(int)
        self.sign = 0

    def set_digit(self, power, digit):
        self.digit_powers[power] = digit

    def add_digit(self, power, digit):
        self.digit_powers[power] += digit
        e, self.digit_powers[power] = divmod(self.digit_powers[power], 10)
        if e:
            self.add_digit(power+1, e)

    def clear_upper_zeros(self):
        while self.digit_powers and self.digit_powers[max(self.digit_powers)] == 0:
            del self.digit_powers[max(self.digit_powers)]

    def fix_sign(self):
        swap_sign = self.digit_powers[max(self.digit_powers)] < 0
        if swap_sign:
            self.sign = 0 if self.sign else 1
            for i in self.digit_powers:
                self.digit_powers[i] *= -1

    def carry_negative

def parse_string(bytes):
    lex = Lexer(bytes.strip())
    [pm, digs_pre_dot, dot, digs_post_dot, exp, exp_power
     ] = lex.pulls(r'[+-]', r'\d+', r'\.', r'\d+', r'[eE]', r'[+-]?\d+')
    if not lex.eof or (exp and not exp_power):
        raise LexicalError("bad sigfig literal %r" % (bytes,))

    sign = 1 if pm == '-' else 0
    power = int(exp_power) if exp_power else 0
    digs_pre_dot = map(int, digs_pre_dot)
    digs_post_dot = map(int, digs_post_dot)

    #remove insignificant trailing zeros
    if not dot:
        if set(digs_pre_dot) > set([0]):
            while digs_pre_dot and digs_pre_dot[-1] == 0:
                digs_pre_dot.pop(-1)
                power += 1

    if not digs_pre_dot:
        digs_pre_dot.append(0)

    #make scientific by shifting digits and power accordingly
    while len(digs_pre_dot) > 1:
        digs_post_dot.insert(0, digs_pre_dot.pop(-1))
        power += 1
    assert len(digs_pre_dot) == 1

    digits = digs_pre_dot + digs_post_dot

    #remove insignificant leading zeros
    while len(digits) > 1 and digits[0] == 0:
        digits.pop(0)
        power -= 1

    digits = tuple(digits)

    return sign, digits, power

__name__ == '__main__' and test_sigfigs(sigfig_class=SigFig)

print SigFig('90.8513423') + SigFig('12.6')
print SigFig('90.8513423') * SigFig('1.343e-8')
