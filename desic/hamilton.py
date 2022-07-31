"""The quaternion module is not for doing rotations. It is for doing the
algebra of "H"-- the quaternions, real or complexified (biquaternions).
"""
#pylint: disable=E1101,W0611

## \namespace geo::desic::hamilton
# (<a href="http://en.wikipedia.org/wiki/Biquaternion">Bi</a>)
# <a href="http://en.wikipedia.org/wiki/Quaternion">Quaternions</a>
# (__C__ \f$\otimes\f$) __H__

from ..utils import cauchy
from ..metric.euler import hamilton as versor
from ..metric.euler.hamilton import parity, Basis

## Mixin for operators on scalar and vector parts (prolly could be
# passed up the chain, but the composite nature of the components
# requires more thought and possibly baggage)
class FourOpsMixin(object):
    """element-wise operations on scalar and vector
    parts of an object: Quaternion or FourVector."""

    ## \f$ {\bf q} + {\bf p} = (q+p, \vec{q} + \vec{p}) \f$
    def __add__(self, other):
        return type(self)(self.scalar + other.scalar,
                          self.vector + other.vector)

    ## \f$ {\bf q} - {\bf p} = (q-p, \vec{q} - \vec{p}) \f$
    def __sub__(self, other):
        return type(self)(self.scalar - other.scalar,
                          self.vector - other.vector)

    ## dilation
    def __rmul__(self, other):
        from operator import mul
        from functools import partial
        return type(self)(*map(partial(mul, other), self.iter()))

    projection = cauchy.projection
    rejection = cauchy.rejection
    reflection = cauchy.reflection
    hat = cauchy.hat


## <a href="http://en.wikipedia.org/wiki/Quaternion">Quaternions</a>
## include addition, subtraction, and products: inner, outer, odd, even.
class Quaternion(versor.QuaternionMixin, FourOpsMixin):
    """Quaternion(Scalar, Vector) -anything in C X H:

    +, - are defined, as well as:
    constant * quaternion

    Additional Products are:
    inner
    outer
    odd
    even

    Extra methods:
    log
    exp
    spinor

    And for Biquaternions there's another type of conjugation:

    complex_conjugate()
    """

    ## This is the full Quaternion Space
    measure = 1

    __metaclass__ = versor.Lebesgue

    ## \returns \f$ \sqrt{q \cdot q^*} \f$
    def __abs__(self):
        return self.dot(self.complex_conjugate()).real ** 0.5

    ## Multiplicative inverse.
    def __invert__(self):
        return self.inv()

    ## Explicit Inverse
    # \returns \f$ {\bf \tilde{q}} = {\bf \bar{q}}/||{\bf q}||^2 \f$
    def inv(self):
        """inverse assumes unit norm"""
        return self.C / self.dot(self)

    ## Grassmann or Dilation by external polymorphism
    def __mul__(self, other):
        try:
            return super(Quaternion, self).__mul__(other.quaternion())
        except AttributeError:
            return self.__rmul__(other)

    ## \f$ q/p \rightarrow q\cdot p^{-1} \f$
    def __div__(self, other):
        return self * (other**(-1))

    ## Unit Quaternion.
    # \returns \f$ {\bf q} \rightarrow \frac{1}{||{\bf q}||} {\bf q} \f$
    def hat(self):
        """q.hat() --> q / ||q||**2"""
        return self / abs(self)

    ##  \f$ \frac{1}{2} ({\bf \bar{q} p} + {\bf \bar{p} q}) \f$
    def inner(self, other):
        """ Quaternion Inner Product: (~qp + ~pq)/2 """
        return ((self.C)*other.quaternion() +
                (other.quaternion().C)*self) / 2

    ## \f$ \frac{1}{2}(\bar{{\bf q}}{\bf p}-\bar{{\bf p}}{\bf q}) \f$
    def outer(self, other):
        """ Quaternion Outer Product: (~qp - ~pq)/2  --> (0; (~qp).vector) """
        return ((self.C)*other.quaternion() -
                (other.C).quaternion()*self) / 2

    ## \f$ \frac{1}{2}({\bf q}{\bf p}-{\bf p}{\bf q}) =
    ## (0; \vec{p} \times \vec{q})  \f$ \n
    def odd(self, other):
        """ Quaternion Odd Product: (qp - pq)/2 = (0, q.v X p.v) """
        return (self*other.quaternion() - other.quaternion()*self) / 2

    __and__ = outer  # is this the right thing to do?
    __xor__ = odd    # ditto.

    ## \f$ \frac{1}{2}({\bf q}{\bf p}+{\bf p}{\bf q})\f$ \n is the Even
    ## product.. the Grassman minus the Odd.
    def even(self, other):
        """ Quaternion Even Product  Quaternion: (qp + pq)/2 """
        return (self*other.quaternion() + other.quaternion()*self) / 2

    ## \f$ {\bf v} = {\bf q}/||{\bf q}|| \f$
    # \returns geo.metric.euler.hamilton.Versor
    def versor(self):
        """Scale into unit sphere- and recaste as a Versor"""
        return versor.Versor(*self.hat().iter())

    ## Hopf construction from geo.metric.euler.hamilton.hopf
    # \throws TypeError
    def AlibiMatrix(self):
        """Type Error"""
        raise TypeError("hopf assumes abs(q) = 1; not assured w/ {}".format(
            str(type(self).__name__)))

    ## Inverse of AlibiMatrix()
    def AliasMatrix(self):
        """Inverse of Alibi: this leads to an existential crisis
        regarding orthonormal coordinates in this package."""
        return self.AlibiMatrix().I

    ## \f$ [q, p] = qp - pq \f$
    def commutator(self, other):
        """q.commutator(p) --> q*p - p*q"""
        o = other.quaternion()
        return (self * o) - (o * self)

    ## Quaternion
    def quaternion(self):
        """self"""
        return self

    ## Quaternion Exponential (generalized Euler's Formula).
    # \returns \f$ \exp{{\bf q}} = (\cos{||\vec{q}||}, \hat{q}
    # \sin{||\vec{q}||}) \f$
    def exp(self):
        """q' = exp(q)"""
        from numpy import cos, sin, exp
        return (
            type(self)(
                type(self.scalar)(cos(abs(self.vector))),
                self.vector.hat()*sin(abs(self.vector))
                )*exp(self.scalar.w)
            )

    ## Quaternion Logarithm (Invers exp()).
    # \returns \f$ \ln{\bf q} = (\ln{||{\bf q}||}, \hat{q}
    # \cos^{-1}{\frac{q}{||{\bf q}||} })  \f$
    def log(self):
        """q' = log(q) is a lot like the log of a rotation matrix:

        (0; (angle/2) * axis) """
        from numpy import log, arccos
        return type(self)(
            type(self.scalar)(log(abs(self).w)),
            self.vector.hat()*arccos(self.scalar.w/abs(self).w)
            )

    ## 2j + 1 dimensional Spinor representation
    # \param \f$ j = \frac{1}{2} \f$ Spin of representation.
    # \returns \f$ (\sigma_0, -i\vec{\sigma})\cdot {\bf q} \f$ as numpy matrix
    def spinor(self, j=0.5):
        """xi = q.spinor(j=0.5)

        xi is a (2j+1) square matrix, Pauli matrix equivalent.
        """
        from .hilbert import pauli
        irr = pauli.Spin(j)
        sigma = 2*irr.J()
        sigma0 = irr.I()
        return (sigma0*self.scalar - 1j*sigma*self.vector).w

    ## \f$ {\bf q}^* = (q^*, \vec{q}^*) \f$
    def complex_conjugate(self):
        """(q, Q)* --> (q*, Q*)"""
        return type(self)(self.scalar.C, self.vector.C)

    ## \return \f$ \bar{q}^* \f$
    @property
    def H(self):
        """q.H -> q*.C (Hermitian conjugation)."""
        return self.complex_conjugate().C


## The Quaternion Basis
BASIS = tuple([item.quaternion() for item in versor.BASIS])

## The Unit Quaternions (not stuck on the unit-ball).
W, I, J, K = BASIS

## Synonyms.
QUATERNIONS = BASIS

