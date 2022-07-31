"""A Cayley Dickson Extension Algebra Factory:

For example:

r --> r + ir' = c     Pairs of Real Numbers make Complex Numbers
c --> c + jc' = q     Pairs of Real Numbers make Quaternions
q --> q + lq' = o     ...Octontions, ...etc...


The rules for extending have a sign dependence, so their are 2 classes:

Extension
Split
"""

## \namespace geo.desic.cayley_dickson
# <a href="http://en.wikipedia.org/wiki/Cayley-Dickson_construction">Cayley
# Dickson Construction</a>.
import abc
import functools
import operator
import itertools
import math

## Match involution with phase function.
_hash = {1: math.acos, -1: math.acosh}


## <a href="http://en.wikipedia.org/wiki/Associator">The associator</a>.
# \param x An element of the algebra
# \param y An element of the algebra
# \param z An element of the algebra
# \return \f$[x, y, z] = (xy)z - x(yz) \f$
def associator(x, y, z):
    """[x, y, z] => (xy)z - x(yz)

    is zero in Alternating Algebras."""
    return (x*y)*z - x*(y*z)


## Unary operator decorator
# \param method A unary operation
# \returns umethod wrapped to return type matching to that which it is bound.
def unary(method):
    """op(x) --> type(x)(op(x))"""

    @functools.wraps(method)
    def umethod(self):
        """call method and recast as type(self)"""
        return type(self)(method(self.a), method(self.b))

    return umethod


## Binary operator decorator
# \param method A binary operation
# \returns bmethod wrapped to return type matching to that which it is bound.
def binary(method):
    """op(x, y) --> type(x)(op(x, y))"""

    @functools.wraps(method)
    def bmethod(self, other):
        """call method and recast as type(self)"""
        return type(self)(method(self.a, other.a), method(self.b, other.b))

    return bmethod


## Cayley Dickson Extension Base Class
class CayleyDicksonConstructor(object):
    """Caley-Dickson Algebra base class."""

    ## The Algebra is not known at this level.
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def _gamma(self):
        """The involution."""

    ## Construct from 2 elements of root algebra
    def __init__(self, a, b=None):
        """z = Algebra(a [,b=0])"""
        ## 'real' part
        self.a = a
        ## imaginary part
        self.b = 0.*a if b is None else b

    ## +x
    __pos__ = unary(operator.pos)

    ## -x
    __neg__ = unary(operator.neg)

    ## Conjugation is the key
    # \returns \f$ \tilde{(a, b)} = (\tilde a, -b) \f$
    def conjugate(self):
        """Conjugation"""
        return type(self)(self.a.conjugate(), -self.b)

    ## Inverse, via methods that work for built-ins.
    # \return \f$ z^{-1} = \frac{\tilde z}{||z||^2} \f$
    def inv(self):
        """Inverse."""
        return self.conjugate() / abs(self)**2

    ## ~x --> x.inv()
    __invert__ = inv

    ## x + y
    __add__ = binary(operator.add)

    ## x - y
    __sub__ = binary(operator.sub)

    ## The float is the float of "a".
    def __float__(self):
        return float(self.a)

    ## \f$ ||z|| = \sqrt{z z^*} \f$
    def __abs__(self):
        return abs((self * self.conjugate()).a) ** 0.5

    ## A eval'able repr.
    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__,
                                   self.a, self.b)

    ## WETT Dilation, again.
    # \return \f$ c(a, b) = (ca, cb) \f$
    def __rmul__(self, c):
        return type(self)(c*self.a, c*self.b)

    ## Multiplication (with involution)
    # \param other Another element of the algebra.
    # \returns tuple \f$ (p, q)(r, s) = (pr-\gamma s^*, sp+qr^*) \f$
    def _mul(self, other):
        a = self.a*other.a - self._gamma*other.b.conjugate()*self.b
        b = other.b*self.a + self.b*other.a.conjugate()
        return a, b

    ## Type-check to distinguish group operator from dilation.
    def __mul__(self, other):
        if isinstance(other, type(self)):  # dilation vs group op
            a, b = self._mul(other)
            return type(self)(a, b)
        return self.__rmul__(other)

    ## Type-check to distinguish group operator from dilation.
    def __div__(self, other):
        if isinstance(other, type(self)):  # dilation vs group op
            return self * pow(other, -1)
        return self.__rmul__(1./other)

    ## Can you exponentiate a non-associative algebra?
    def __pow__(self, n):
        if n == -1:  # power branch
            return self.conjugate() / abs(self) ** 2
        elif n < 0:
            return pow(pow(self, -1), -n-1)
        elif n == 0:
            u = self.a + self.b
            return type(self)(pow(u, 0), 0*u)
        else:
            return reduce(operator.mul, itertools.repeat(self, n))

    ## Dimension of the Algebra
    def dim(self):
        """Dimension of Algebra."""
        n = 0
        r = self
        while True:
            try:
                r = r.a
            except AttributeError:
                n += isinstance(r, complex)
                break
            n += 1
            continue
        return 2 ** (n-1)

    ## Normalizer
    # \returns \f$ \frac{z}/||z|| \f$
    def hat(self):
        """normed element"""
        return self / abs(self)

    ## A phase factor,
    def phase(self):
        """A phase factor"""
        func = _hash[self._gamma]
        return (self.dim() / 2) * func(self.hat().a)


## TBD
def power(self, other):
    """Raise to a power."""
    theta = self.phase()
    ntheta = theta * other
    return pow(abs(self), other) * type(self)(
        math.cos(ntheta), math.sin(ntheta)
        )


## Extensions Algebra
class Extension(CayleyDicksonConstructor):
    """Extensions algebra. See __init__ for signature."""
    ## \f$ \gamma = +1 \f$ defines a nominal Extension.
    _gamma = 1


## Split Extensions Algebra
class Split(CayleyDicksonConstructor):
    """Split Extensions. See __init__ for signature."""

    ## \f$ \gamma = -1 \f$ defines a so-called Split Extension.
    _gamma = -1


## \f$ 1 \f$ (Complex)
l = Extension(1., 0.)

## \f$ i \f$ (Complex)
i = Extension(0., 1.)

## \f$ 0 \f$ (Complex)
O = 0. * l


## \f$ w \f$ (Hyper-complex)
W = Extension(l, O)

## \f$ i \f$ (Hyper-complex)
I = Extension(i, O)

## \f$ j \f$ (Hyper-complex)
J = Extension(O, l)

## \f$ k \f$ (Hyper-complex)
K = Extension(O, i)

## \f$ 0 \f$ (Hyper-complex)
N = 0 * W

## \f$ e_0 \f$, the Octonion identity
e0 = Extension(W, N)
## \f$ e_1 \f$ (Octonion)
e1 = Extension(I, N)
## \f$ e_2 \f$ (Octonion)
e2 = Extension(J, N)
## \f$ e_3 \f$ (Octonion)
e3 = Extension(K, N)
## \f$ e_4 \f$ (Octonion)
e4 = Extension(N, W)
## \f$ e_5 \f$ (Octonion)
e5 = Extension(N, I)
## \f$ e_6 \f$ (Octonion)
e6 = Extension(N, J)
## \f$ e_7 \f$ (Octonion)
e7 = Extension(N, K)
