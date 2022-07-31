"""This module defines the alebraic operations on spherical tensors:

The Basics:

+v
-v
u + v
u - v

T/a

and the problem of multiplicative dispatch is handled (poorly) here:

u * v -->  u.dot(v)
a * v -->  u.dilate(v)

u ^ v --> u.cross(v)

u & v --> u.outer(v)

dot, cross, and outer are addressed in the spherical basis,
which is little different.

The projection operator is defined a little differently from Cartesian,
as it addmits an integer 'm' argument:

pure_state = tensor >> m

is tensor[j=n, q=0, m=m] times |n,0,m>.
"""

## \namespace geo.metric.wigner.eckart.ka Algebraic Operations of
# Spherical Tensors

import abc
import operator


## \f$ i \f$
I = (1j)  # sympy

from .one import cross, dot, outer

## Operator Decorator.
def check(op):
    """check(op) wraps a binary operation, op, that works on like-rank
    spherical tensors with error message diagonstics."""

    ## Fgt.
    def checked_op(self, other):
        """binary op, with error handling"""
        try:
            result = op(self, other)
        except TypeError as err:
            try:
                drank = self.L == other.L
            except AttributeError:
                from ....utils.exceptions import NonGeometricTypeError as err
                msg = "expected L={} spherical tensors, got: {}".format(
                    self.L, other)
            else:
                if not drank:  # raise
                    from ....utils.exceptions import (
                        NonCovariantOperation as err)
                    msg = "can't combine L={} an L={}".format(
                        self.L, other.L)
                else:
                    msg = "operand failed with type {}".format(type(other))
            raise err(msg)
        return result
    return checked_op


## Base class for algebraic operations on Spherical Tensors.
class SphericalTensor(object):
    """Algebraic Operation."""


    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def _mul(self, other):
        """Dispatch multiplication?"""

    @abc.abstractmethod
    def conjugate(self):
        """Complex Conjugate"""

    @abc.abstractmethod
    def tocart(self):
        """This indicates and MRO problem, in that it overides..."""

    ## +T
    def __pos__(self):
        return self._pack(*map(operator.pos, self))

    ## -T
    def __neg__(self):
        return self._pack(*map(operator.neg, self))

    ## T + T'
    # \throws exceptions.NonCovariantOperation
    # \throws exceptions.NonGeometricTypeError
    # \throws exceptions.TypeError
    @check
    def __add__(self, other):
        return self._pack(*map(operator.add, self, other))

    ## T - T'
    @check
    # \throws exceptions.NonCovariantOperation
    # \throws exceptions.NonGeometricTypeError
    # \throws exceptions.TypeError
    def __sub__(self, other):
        return self._pack(*map(operator.sub, self, other))

    def __rmul__(self, other):
        return self._pack(*[other*item for item in self])

    def __div__(self, other):
        return self.__rmul__(1./other)

    def __mul__(self, other):
        try:
            L = other.L
        except AttributeError:
            return self.__rmul__(other)
        return self._mul(other)

    def __invert__(self):
        return self.conjugate()

    C = property(__invert__)

    def __abs__(self):
        return self.dot(self) ** 0.5

    ## Let Cartesian logic handle pow(T, n)
    def __pow__(self, n):
        return pow(self.tocart(), n).spherical()

    ## Real Part
    @property
    def real(self):
        """Re(T)"""
        return self._pack(*[item.real for item in self])

    ## Imaginary Part
    @property
    def imag(self):
        """Im(T)"""
        return self._pack(*[item.imag for item in self])

    ## TODO: check conjugate
    ## Dot Product.
    # \param self \f$ {\bf A} \f$
    # \param other \f$ {\bf B} \f$
    # \returns \f$ \langle A|B\rangle =
    # \sum_{q=-1}^1{A_q B_q^*} \f$
    def dot(self, other):
        """Full contraction on all weights."""
        from . import Scalar as cls
        return cls(
            reduce(
                operator.add,
                (self[j, q, m] * other[j, q, m] for j, q, m in self.ijqm())))


class Scalar(SphericalTensor):
    """Scalar(w)"""

    ## Complex Conjugate:
    # \returns \f$w^* = w^* \f$
    def conjugate(self):
        """Complex conjugate."""
        return type(self)(self[0].conjugate())

    def _mul(self, other):
        return self[0] * other


## Spherical Vector Arithmetic/Algebra
class Vector(SphericalTensor):
    """Vector(m=0, m=1, m=-1)"""

    ## u & v --> u.outer(v)
    def __and__(self, other):
        return self.outer(other)

    ## Let Cartesian do the product.
    def _outercart(self, other):
        return (self.vector & other.vector).spherical()

    ## Complex Conjugate:
    # \returns \f$v^* = -v^*_{-1}{\bf \hat{e}}^+ + v^*_0{\bf \hat{e}}^0
    # -v^*_{1}{\bf \hat{e}}^ - \f$
    def conjugate(self):
        """Complex conjugate."""
        return type(self)(self[0].conjugate(),
                          -self[-1].conjugate(),
                          -self[1].conjugate())

    def _mul(self, other):
        if other.L == 1:  # forced polymorphism
            return self.dot(other)
        # note: this should pass it up __rmul__ (dilation).
        return NotImplemented

    ## ::cross product (purely spherical implementation).
    cross = cross

    ## ::dot product (purely spherical implementation).
    dot = dot

    ## ::outer product (purely spherical implementation).
    outer = outer
    
    ## u ^ v --> u.cross(v)
    def __xor__(self, other):
        return self.cross(other)

    ## Project out "m" values from highest weight
    # \param m An order integer
    # \returns Vector projection onto \f$ {\bf \hat{e}}^m \f$
    def __rshift__(self, other):
        try:
            m = int(other)
        except TypeError:  # external poly :-(
            result = NotImplemented
        else:
            args = [0, 0, 0]
            args[m] = self[m]
            result = type(self)(*args)
        return result


## Spherical Tensor Arithmetic/Algebra
class Tensor(SphericalTensor):
    """Tensor(scalar, [vector x 3], [tensor x 5])"""

    ## Transpose flips sign of ODD degrees.
    @property
    def T(self):
        """Transpose just flips sign of odd degrees."""
        return type(self)(
            *[[pow(-1, n*i) for i in l] for n, l in zip(range(3), self)]
        )

    ## Conjugation (spherical style)
    # \returns eckart.Tensor via
    # \f$ |l,m\rangle\rightarrow (-1)^l |l,-m\rangle\f$
    def conjugate(self):
        """a|l, m> --> (a^*)(-1)**l |l, -m>"""
        return type(self)(
            self[0, 0].conjugate(),
            [self[1, 0].conjugate(),
             -self[1, -1].conjugate(),
             -self[1, 1].conjugate()],
            [self[2, 0].conjugate(),
             -self[2, -1].conjugate(),
             self[2, -2].conjugate(),
             self[2, 2].conjugate(),
             -self[2, 1].conjugate(),
             self[2, 2].conjugate()])

    ## Project out j, m values
    def __rshift__(self, other):
        try:
            j, m = other
        except TypeError:
            raise ValueError("needs t and m")
        return self.basis()[j, 0, m] * self[j, 0, m]

    def _mul(self, other):
        if other.L == 1:
            print "Posterior not-implemented"
        elif other.L == 2:
            print "T*T' not implemented"
        else:
            print "rank {} not imp'd".format(other.L)


## Spherical Tensor Arithmetic/Algebra (rank 3)
class Three(SphericalTensor):
    """Math for rank 3"""

    def _mul(self, other):
        if other.L == 1:
            print "Posterior not-implemented"
        elif other.L == 2:
            print "T*T' not implemented"
        else:
            print "rank {} not imp'd".format(other.L)

    def conjugate(self):
        """NotImplemented"""


## Spherical Tensor Arithmetic/Algebra (rank 4)
class Four(SphericalTensor):
    """Math for rank 4"""

    def _mul(self, other):
        if other.L == 1:
            print "Posterior not-implemented"
        elif other.L == 2:
            print "T*T' not implemented"
        else:
            print "rank {} not imp'd".format(other.L)

    def conjugate(self):
        """NotImplemented"""
