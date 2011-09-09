

raise RuntimeError("use sympy for this when ready")

from __future__ import division
from __future__ import absolute_import

from . import algebra as A

class SymbolicBase(A.DivAlgebraBase):

    pass

class Symbol(SymbolicBase):

    def __init__(self, print_form):
        self.print_form = print_form

class Add(SymbolicBase):

    def __init__(self, terms):
        self.terms = terms
