'''Class SigFig for repressenting digits of a measurement that
   convey meaningful inormation (i.e. are significant).
'''

from __future__ import absolute_import

from decimal import Decimal
from collections import defaultdict

from hlab.lexing import Lexer, LexicalError
from hlab.bases import AutoRepr

valid_digits = tuple(range(10))

class SigFig(AutoRepr):

    def __init__(self, arg):
        if isinstance(arg, unicode):
            arg = str(arg)
        elif isinstance(arg, Decimal):
            sign, digits, exp = arg.as_tuple()
            arg = (sign, digits, exp - 1 + len(digits))
        if isinstance(arg, str):
            arg = parse_string(arg)
        elif isinstance(arg, (int,long)):
            arg = parse_string(str(arg))
        elif isinstance(arg, SigFig):
            arg = (arg.sign, arg.digits, arg.power)

        sign, digits, power = arg
        digits = tuple(digits)

        assert sign in (0,1)
        assert all(isinstance(digit, (int,long)) and digit in valid_digits
                   for digit in digits)
        assert isinstance(power, (int,long))

        self.sign = sign
        self.digits = digits
        self.power = power

    def repr_args(self):
        return [str(self)]

    #legacy method
    def as_scientific(self):
        return self

    def as_decimal(self):
        return Decimal((self.sign, self.digits, self.least_significant_place))

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

    def round_to_sigfigs(self, n):
        if n <= 0:
            raise ValueError("rounding %s to invalid number of sig figs %d" % (self, n))
        return self.round_at_index(n)

    def round_to_place(self, n):
        #print 'round to', n
        return self.round_at_index(self.power - n + 1)

    def round_at_index(self, i):
        #print 'round', self, self.sign, self.digits, self.power, 'at', i
        if i < 0:
            return self.__class__('0')

        if i >= len(self.digits):
            return self.shift_sigfigs(-(i - len(self.digits)))

        digits = list(self.digits)
        power = self.power

        if i==0:
            digits.insert(0, 0)
            power += 1
            i += 1

        round_dig = digits[i]

        if (round_dig > 5) or (round_dig == 5 and digits[i-1]%2):
            digits[i-1] += 1

        # carray
        for j in xrange(len(digits) - 1, 0, -1):
            if digits[j] == 10:
                digits[j] = 0
                digits[j-1] += 1

        digits = digits[:i:]

        if digits[0] == 10:
            digits[0] = 0
            digits.insert(0, 1)
            power += 1

        if digits[0] == 0:
            del digits[0]
            power -= 1

        if not digits:
            return SigFig('0')

        return self.__class__((self.sign, digits, power))

    def shift_sigfigs(self, sigfigs=1):
        if sigfigs == 0:
            return self

        if sigfigs > 0:
            digits = self.digits[:-sigfigs:]
            if not digits:
                return self.__class__(0)
        else:
            digits = self.digits + (0,) * -sigfigs

        return self.__class__((self.sign, digits, self.power))

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
            return self.get_exp_format_args()

        digs = map(str, self.digits)

        #deal with trailing zeros
        if sigfigs > power:
            insignificant_zeros = 0
            #prevet trailing decimal, i.e. 10.
            if power > 0 and len(digs) == power+1 and digs[-1] == '0':
                return self.get_exp_format_args()
        else:
            #check for significant trailing zeros
            if self.digits[-1::] == (0,):
                return self.get_exp_format_args()
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
            return self.get_exp_format_args()
        assert self.digits == (0,)
        if self.power == 0:
            return '0', None
        return '0.' + '0' * -self.power, None

    def get_exp_format_args(self):
        return ['%s%d%s' % ('-' if self.sign else '',
                            self.digits[0],
                            '.' + ''.join(map(str, self.digits[1::])) if
                            len(self.digits) > 1 else ''),
                self.power]

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__((1 if self.sign==0 else 0,
                               self.digits, self.power))

    def __nonzero__(self):
        return self.as_decimal() != 0

    def perform_binary_operation(self, other, func, rule):
        if not isinstance(other, (SigFig, Decimal, int, long)):
            return NotImplemented

        selfd = self.as_decimal()
        otherd = other.as_decimal() if isinstance(other, SigFig) else other
        valued = func(selfd, otherd)
        value = self.__class__(valued)

        sigfigs = (min(self.sigfigs, other.sigfigs)
                   if isinstance(other, SigFig) else
                   self.sigfigs)

        if rule=='mul':
            if valued == 0:
                value.digits = (0,) * sigfigs
            else:
                value = value.round_to_sigfigs(min(self.sigfigs, other.sigfigs)
                                               if isinstance(other, SigFig) else
                                               self.sigfigs)
        elif rule=='add':
            lsp = (max(self.least_significant_place, other.least_significant_place)
                   if isinstance(other, SigFig) else
                   self.least_significant_place)
            value = value.round_to_place(lsp)
#            if value.digits == (0,):
#                value.digits = (0,) * max(1, sigfigs)

        else:
            raise ValueError("bad sigfig rule %r" % (rule,))

        return value

    def __mul__(self, other):
        return self.perform_binary_operation(other, lambda a,b: a*b, 'mul')
    def __rmul__(self, other):
        return self.perform_binary_operation(other, lambda a,b: b*a, 'mul')
    def __div__(self, other):
        return self.perform_binary_operation(other, lambda a,b: a/b, 'mul')
    def __rdiv__(self, other):
        return self.perform_binary_operation(other, lambda a,b: b/a, 'mul')
    def __mod__(self, other):
        return self.perform_binary_operation(other, lambda a,b: a%b, 'mul')
    def __rmod__(self, other):
        return self.perform_binary_operation(other, lambda a,b: b%a, 'mul')
    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __add__(self, other):
        return self.perform_binary_operation(other, lambda a,b: a+b, 'add')
    def __radd__(self, other):
        return self.perform_binary_operation(other, lambda a,b: b+a, 'add')
    def __sub__(self, other):
        return self.perform_binary_operation(other, lambda a,b: a-b, 'add')
    def __rsub__(self, other):
        return self.perform_binary_operation(other, lambda a,b: b-a, 'add')

    def __gt__(self, other):
        return (self-other).as_decimal() > 0
    def __ge__(self, other):
        return (self-other).as_decimal() >= 0
    def __eq__(self, other):
        if not isinstance(other, (int,long,Decimal,SigFig)):
            return NotImplemented
        try:
            return not (abs((self-other).as_decimal()) > 0)
        except TypeError,e:
            return False
    def __ne__(self, other):
        return not (self == other)
    def __le__(self, other):
        return (self-other).as_decimal() <= 0
    def __lt__(self, other):
        return (self-other).as_decimal() < 0

    def __pow__(self, op):
        if not isinstance(op, (int,long)):
            raise TypeError("pow not implemented for non-integer type %r" % (op,))
        if op < 0:
            return 1 / self**-op
        if op == 0:
            return self / self
        return (self.as_decimal() ** op).round_to_sigfigs(self.sigfigs)

    def sqrt(self):
        return self.__class__(self.as_decimal().sqrt()).round_to_sigfigs(self.sigfigs)



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

#sf = SigFig('120')
#print sf.least_significant_place
#m = sf - 120
#print m
# print SigFig('0.00000') - 0

# sf = SigFig('22.534')
# print sf
# print sf.round_to_sigfigs(4)
# print SigFig(sf.as_decimal())
# print
# print SigFig('90.8513423') + SigFig('12.6')
# print SigFig('90.8513423') * SigFig('1.343e-8')

print SigFig('10') + SigFig('-5')
