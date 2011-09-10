'''Class SigFig for repressenting digits of a measurement that
   convey meaningful inormation (i.e. are significant).

   Current implementation is rather messy and have found we can
   get away simpler data structure to represent the same information.

   None the less, this implementation has been tested and found correct
   for a large class of cases.  Additionally, this file includes some
   simple and dirty testing to be ran as the code is evolved. Constant
   testing is very important as the SigFig class is the core component
   for the chemistry and physics calculations provided in this library.
'''


from decimal import Decimal

from hlab.lexing import Lexer, LexicalError
from hlab.bases import AutoRepr

class SigFig(AutoRepr):

    def __init__(self, arg):
        if isinstance(arg, unicode):
            arg = str(arg)
        if isinstance(arg, str):
            arg = parse_string(arg)
        elif isinstance(arg, Decimal):
            arg = convert_decimal(arg)
        elif isinstance(arg, (int,long)):
            arg = convert_integer(arg)
        elif isinstance(arg, SigFig):
            arg = (arg.sign, arg.digs_pre_dot, arg.dot, arg.digs_post_dot, arg.exp, arg.exp_power)
        (self.sign,
         self.digs_pre_dot,
         self.dot,
         self.digs_post_dot,
         self.exp,
         self.exp_power
         ) = arg

    def repr_args(self):
        return [str(self)]

    def as_tuple(self):
        return (self.sign, self.digs_pre_dot, self.dot, self.digs_post_dot, self.exp, self.exp_power)

    def __str__(self):
        base, exp = self.get_format_args()
        if exp is None:
            return base
        return '%se%s' % (base, exp)

    def get_format_args(self, min_exp_power=3):
        return self.as_scientific()._get_format_args(min_exp_power)

    def get_exp_format_args(self):
        return self.as_scientific()._get_exp_format_args()

    def _get_format_args(self, min_exp_power):
        if min_exp_power!=None and abs(self.exp_power) > max(self.sigfigs, min_exp_power):
            return self._get_exp_format_args()
        digs = map(str, self.digs_pre_dot + self.digs_post_dot)
        if self.sigfigs <= self.exp_power:
            #check for trailing zeros
            if self.digs_post_dot[-1::] == (0,):
                return self._get_exp_format_args()
            digs += ['0'] * (1 + self.exp_power - self.sigfigs)
        #prevet trailing decimal, i.e. 10.
        elif self.exp_power > 0 and len(digs) == self.exp_power+1 and digs[-1] == '0':
            return self._get_exp_format_args()
        else:
            if digs == ['0']:
                return '0', None
            if self.exp_power >= 0:
                if not (1+self.exp_power==len(digs) and digs[-1] != '0'):
                    digs.insert(1+self.exp_power, '.')
            else:
                digs = ['0.'] + ['0'] * (-1 - self.exp_power) + digs
        return ['%s%s' % ('-' if self.sign else '', ''.join(digs)),
                None]

    def _get_exp_format_args(self):
        return ['%s%d%s' % ('-' if self.sign else '',
                            (self.digs_pre_dot or [0])[0],
                            '.' + ''.join(map(str, self.digs_post_dot)) if
                            self.digs_post_dot else ''),
                self.exp_power]


    def as_decimal(self):
        x = self.as_scientific()
        return Decimal((x.sign,
                        x.digs_pre_dot + x.digs_post_dot,
                        x.exp_power - len(x.digs_post_dot)))

    def shift_sigfigs(self, sigfigs=1):
        s = self.as_scientific()
        if s is self:
            s = self.__class__(self)
        assert s is not self
        s.dot = s.exp = True
        if sigfigs > 0:
            s.digs_post_dot += (0,) * sigfigs
        elif -sigfigs > len(s.digs_post_dot):
            raise ValueError("bad shift")
        else:
            s.digs_post_dot = s.digs_post_dot[:sigfigs]
        return s

    @staticmethod
    def strip_leading_zeros(seq):
        i = iter(seq)
        for el in i:
            if el != 0:
                break
        else:
            return []
        return [el] + list(i)

    _sigfigs = None
    @property
    def sigfigs(self):
        '''cache this calculated property as accessed a lot!
        '''
        return self._calculate_sigfigs()
        #if self._sigfigs is None:
        #    self._sigfigs = self._calculate_sigfigs()
        #return self._sigfigs

    def _calculate_sigfigs(self):
        #handle tricky case of zero
        if list(set(self.digs_pre_dot) | set(self.digs_post_dot)) == [0]:
            if self.exp_power == 0:
                return max(1, len(self.digs_post_dot))
            return len(self.digs_pre_dot) + len(self.digs_post_dot)
        if self.dot:
            if set(self.digs_pre_dot) <= set([0]):
                return len(self.strip_leading_zeros(self.digs_post_dot))
            return (len(self.strip_leading_zeros(self.digs_pre_dot)) +
                    len(self.digs_post_dot))
        return len(self.strip_leading_zeros(self.digs_pre_dot[::-1]))

    @property
    def least_significant_place(self):
        if self.dot:
            return self.exp_power - len(self.digs_post_dot)
        n = len(self.digs_pre_dot)
        if tuple(set(self.digs_pre_dot)) == (0,):
            return 0
        while n and self.digs_pre_dot[n-1]==0:
            n -= 1
        return len(self.digs_pre_dot) - n + self.exp_power

    def as_scientific(self):
        if len(self.digs_pre_dot)==1 and self.digs_pre_dot != (0,):
            return self
        # strip trailing zeros
        trailing_zeros = 0
        digs_pre_dot = list(self.digs_pre_dot)
        digs_post_dot = list(self.digs_post_dot)
        if not self.dot:
            while digs_pre_dot and digs_pre_dot[-1] == 0:
                trailing_zeros += 1
                digs_pre_dot.pop(-1)
            while digs_post_dot and digs_post_dot[-1] == 0:
                digs_post_dot.pop(-1)
        #combine all digits
        exp_power = self.exp_power + trailing_zeros + len(digs_pre_dot) - 1
        digs = list(digs_pre_dot + digs_post_dot)
        #strip leading zeros
        while digs and digs[0] == 0:
            digs.pop(0)
            exp_power -= 1
        #handle zero
        if not digs:
            sigfigs = self.sigfigs
            return self.__class__((0, (0,), self.dot and sigfigs > 1,
                                   (0,) * sigfigs if self.dot else (),
                                   False, 0))
        return self.__class__((self.sign, tuple(digs[:1:]), len(digs) > 1,
                               tuple(digs[1:]), True, exp_power))

    def round_to_sigfigs(self, n):
        if n<1:
            raise ValueError("bad number of sigfigs to round to %r" % (n,))
        s = self.as_scientific()
        if list(s.digs_pre_dot) != [0]:
            n -= 1
        return self.round_scientific_to_index(n, s)

    def round_to_place(self, n):
        s = self.as_scientific()
        return self.round_scientific_to_index(max(-1, -(n - s.exp_power)), s)

    def round_scientific_to_index(self, index, s=None):
        if s is None:
            s = self.as_scientific()
        digs = list(s.digs_pre_dot + s.digs_post_dot)
        index += 1
        if len(digs) <= index:
            while len(digs) < index:
                digs.append(0)
        else:
            if index==0:
                digs.insert(0,0)
                s.exp_power += 1
                index += 1
            round_dig = digs[index]
            if (round_dig > 5) or (round_dig == 5 and digs[index-1]%2):
                digs[index-1] += 1
            del digs[index::]
            #carry
            for i in xrange(len(digs)-1, 0, -1):
                if digs[i]==10:
                    digs[i] = 0
                    digs[i-1] += 1
            if digs[0]==10:
                digs[0] = 0
                digs.insert(0, 1)
                s.exp_power += 1
        s.digs_pre_dot = (digs[0],)
        s.digs_post_dot = tuple(digs[1:])
        return s

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__((int(not self.sign), self.digs_pre_dot, self.dot,
                               self.digs_post_dot, self.exp, self.exp_power))

    def __nonzero__(self):
        return self.as_decimal() != 0

    def __hash__(self):
        return hash(self.as_decimal())

    def perform_binary_operation(self, other, func, rule):
        selfd = self.as_decimal()
        if not isinstance(other, (SigFig, Decimal, int, long)):
            return NotImplemented
        otherd = other.as_decimal() if isinstance(other, SigFig) else other
        valued = func(selfd, otherd)
        value = self.__class__(valued)
        sigfigs = (min(self.sigfigs, other.sigfigs)
                   if isinstance(other, SigFig) else
                   self.sigfigs)
        if rule=='mul':
            value = value.round_to_sigfigs(min(self.sigfigs, other.sigfigs)
                                           if isinstance(other, SigFig) else
                                           self.sigfigs)
            if valued == 0:
                value.digs_pre_dot = ()
        elif rule=='add':
            lsp = (max(self.least_significant_place, other.least_significant_place)
                   if isinstance(other, SigFig) else
                   self.least_significant_place)
            value = value.round_to_place(lsp)
            if value.digs_pre_dot == (0,) and value.digs_post_dot == ():
                value.digs_post_dot = (0,) * max(0, -lsp)
                value.digs_pre_dot = (0,) * max(0, sigfigs-max(0, -lsp))
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
        acc = self
        for i in xrange(op-1):
            acc = acc * self
        return acc

    def sqrt(self):
        return self.__class__(self.as_decimal().sqrt()).round_to_sigfigs(self.sigfigs)


def parse_string(bytes):
    lex = Lexer(bytes.strip())
    [pm, digs_pre_dot, dot, digs_post_dot, exp, exp_power
     ] = lex.pulls(r'[+-]', r'\d+', r'\.', r'\d+', r'[eE]', r'[+-]?\d+')
    if not lex.eof or (exp and not exp_power):
        raise LexicalError("bad sigfig literal %r" % (bytes,))
    if not digs_pre_dot:
        digs_pre_dot = '0'
    sign = 1 if pm == '-' else 0
    digs_pre_dot = tuple(map(int, digs_pre_dot))
    digs_post_dot = tuple(map(int, digs_post_dot))
    dot = bool(dot)
    exp = bool(exp)
    exp_power = int(exp_power) if exp_power else 0
    #print 'PARSE', bytes, (sign, digs_pre_dot, dot, digs_post_dot, exp, exp_power)
    return (sign, digs_pre_dot, dot, digs_post_dot, exp, exp_power)

def convert_decimal(d):
    sign, digits, exp = d.as_tuple()
    return (sign, digits, True, (), True, exp)

def convert_integer(i):
    #cheap trick, but works
    return parse_string(str(i))


__name__ == '__main__' and main()
