"""The scalar module defines geometric objects that
are invariant under rotation. The class:

Scalar(w)

represents the rank-0 tensor object. It is NAN. Not a number, and a number
is not it.

"w" can be a number or a numpy array of any dimensions, in which case
Scalar(w) represents a scalar field or time-dependent scalar.
"""

#pylint: disable=C0301

## \namespace geo::metric::euclid.scalar The Rank-0 Cartesian Tensor: The
#Scalar is not _just_ a number.
import functools
import itertools
import operator
import sys

from . import euclid

__all__ = ('ZERO', 'ONE', 'Scalar')


## Sign
# \param s Something with a sign.
# \returns sgn(s)
def sgn(s):
    """Sign of argument, see pos."""
    return pos(s) - pos(-s)


## Positivity
# \returns \f$ s \ge 0 \f$ by hook or by crook.
def pos(s):
    """really is not negative."""
    try:
        return int(s.w >= 0)
    except TypeError:
        from numpy import array
        return array(map(pos, s))


## sad little scalar * quaternion--this would  be fixed by multimethods.
def _dilate_q(s, q):
    return type(q)(s*q.scalar, s*q.vector)


## Experimental -will need c-bindings to beat numpy
def swrap(func):
    """Experimental scalar maker decorator."""

    def method(self):
        """func as scalar?"""
        return type(self)(func(self.w))
    return method


## a decorator for
# <a href=
# "http://docs.python.org/reference/datamodel.html?highlight=__cmp__#object.__cmp__">
# __cmp__ operators</a> to try scalars, and then do numbers, works for
# singletons and np.ndarrays.
def _cmpdec(func):
    """decorator for cmp operators that need to know about 'w'"""

    @functools.wraps(func)
    def local_cmp(self, other):
        """docstring"""
        try:
            result = func(self.w, other.w)
        except AttributeError:
            result = func(self.w, other)
        return result
    return local_cmp


## A decorator: "dub" not "dub" -- the purpose is to help scalar operations
# with np.arrays--it seems to be impossible to cover every case, it
def _wnotw(func):
    """element-wise should decorate just fine, but it fails If you add
    an plain nd array on the right -- should that be allowed?--it is If you
    decorate with this.
    """

    @functools.wraps(func)
    def wfunc(self, right):
        """operator with Scalar checking"""
        try:
            result = (
                func(self.w, right) if euclid.rank(right) is None else  # guard
                func(self.w, right.w)
                )
        except AttributeError:
            from ...utils.exceptions import NonCovariantOperation
            from ...utils.exceptions import error_message
            raise NonCovariantOperation(error_message(func, self, right))
        return Scalar(result)
    return wfunc


## N dimensional
# <a href="http://en.wikipedia.org/wiki/Scalar_(physics)">Scalar</a>
# class transforms as \f$ s' = s \f$ and \f$ \hat{P}\, s = s \f$
class NScalar(euclid.Tensor_):
    """N-dimensional Scalar Base Class.

    N could be 3, 3+1, ...
    """

    ## Tensor rank
    rank = 0

    ## Pos: I forgot why I needed this
    pos = pos

    ## Sgn: IBID.
    sgn = sgn

    ## Scalars are deeply complicated in pure python: one would like them
    # to act like python numbers (float, complex), but the thing is, and
    # instance could also be a numpy array (of any shape)--and numpy's C
    # bindings are "deeper" than geo's: hence, it's a battle--and much is
    # won or lost in __getattr__.
    # \param attr An attribute that is not in (s)Scalar.__dict__
    def __getattr__(self, attr):
        """self.attr --> self.w.attr

        very complicated by numpy, and ipython-- is way out there too.

        Nevertheless, it is this method that allows you to use
        numpy and not-numpy interchangeably, though:
        numpy.array * Tensor can still cause trouble- as that turns
        the overload over to numpy."""
        try:
            # step one is ask 'w' for the attribite
            return getattr(self.w, attr)
        except AttributeError:
            # Ipython asks for this.
            if attr == '_ipython_display_':  # ipython leaky abstraction
                # someday we can add all kinds of fun.
                # http://nbviewer.ipython.org/github/ipython/ipython/blob/2.x/examples/Notebook/Custom%20Display%20Logic.ipynb
                pass
            else:
                print >> sys.stderr, (
                    "DEBUG %s failed on %s" % (type(self.w).__name__,
                                               attr))
            try:
                import numpy
                try:
                    f = getattr(numpy, attr)
                except AttributeError:
                    pass
                else:
                    print >> sys.stderr, "Debug: saved: ", attr
                    # To get here, you probably called numpy.cos (or similar)
                    # so we return what looks a method that returns
                    # numpy.cos(self.w), for example. There is a problem:
                    # If s is a singelton (array), this returns a python
                    # float (numpy array). Had the lambda cast it to a Scalar,
                    # the float would be promoted, but the array wouldn't.
                    return lambda: f(self.w)
                raise AttributeError(attr)
            except ImportError:
                # you don't have numpy.
                raise AttributeError(
                    "Couldn't ask numpy for: {}".format(attr))

    ## explicit __init__ is just to be nice (and it checks for nested
    # Scalars-- which should not happen).
    # \param w array_like
    def __init__(self, w):
        """Scalar(value)"""
        self.w = w.w if isinstance(w, NScalar) else w  # prevent nesting

    ## Multiplication ?
    def __call__(self, other):
        return self * other

    ## Float
    def __float__(self):
        try:
            return float(self.w)
        except TypeError:
            return self.w.astype(float)
            
    ## Complex
    def __complex__(self):
        try:
            return complex(self.w)
        except TypeError:
            return self.w.astype(complex)



    ## This is  a problem--it's not polymorphic enough- Scalars are a pain:
    # but "no scalars" is a failure--it's numpy- it takes control.
    def __div__(self, other):
        try:
            result = super(NScalar, self).__div__(other)
        except (TypeError, AttributeError):
            try:
                result = super(NScalar, self).__div__(other.w)
            except AttributeError:
                from ...utils.exceptions import (
                    UndefinedGeometricOperation, error_message
                    )
                raise UndefinedGeometricOperation(
                    error_message(type(self).__div__, self, other)
                    )
        return result

    ## number/scalar (use caution for array/scalar)
    def __rdiv__(self, other):
        #try:
        #    len(other)
        #except TypeError:
        #    print >> sys.stderr, 'Warning: numpy has jacked the interface!'
        return type(self)(other/self.w)


    @_wnotw
    def __sub__(self, other):
        return operator.sub(self, other)

    ## Send to "__radd__"...
    def __rsub__(self, other):
        return other + (-self)

    ## This is dicey-- no, not scalar addition, Tensor-space additions:
    # \param other A scalar or a vector
    # \returns Scalar or a quaternion (\f$ {\bf R} \oplus {\bf R^3} \f$)
    def __add__(self, other):
        from .euclid import rank
        if rank(other) == 1:
            # really tensor addition of R and R3 --> H.
            return self._v2q(other)
        # normal addition.
        return self._addition(other)

    ## This is for tensor space addition
    # \param vector_ A vector.Vector
    # \returns quaternion a full __H__ desic.hamilton.Quaternion
    def _v2q(self, vector_):
        from ...desic.hamilton import Quaternion
        return Quaternion(self, vector_)

    ## Normal Scalar Addition
    @_wnotw
    def _addition(self, other):
        return operator.add(self, other)

    ## Caution: a number + a scalar = scalar (is that OK?)
    def __radd__(self, other):
        return self + other

    ## pow is pretty regular
    def __pow__(self, other, m=None):
        if m is not None:  # guard on ternary power.
            raise NotImplementedError("Ternary pow is not implemented")
        try:
            result = (self.w)**other
        except TypeError:
            result = (self.w)**(other.w)
        return type(self)(result)

    ## Ironically, I am a geometric animal-- and I give you the
    # "I am not an animal' speech.
    def __mod__(self, other):
        #pylint: disable=W0613,R0201
        from ...utils.exceptions import ScalarAreNotNumbersError
        raise ScalarAreNotNumbersError

    ## See Scalar.__mod__().
    def __floordiv__(self, other):
        #pylint: disable=W0613,R0201
        from ...utils.exceptions import ScalarAreNotNumbersError
        raise ScalarAreNotNumbersError

    ## reflect pow is required, e.g:
    # \f$\phi = e^{i {\bf \vec{k}\cdot\vec{r}}} \f$ is a Scalar
    #  in the exponent (and phase is a scalar, too).
    def __rpow__(self, other):
        return type(self)(other**(self.w))

    ## <a href="http://docs.python.org/library/operator.html"> < </a>
    #  decorated with cmpdec() .
    @_cmpdec
    def __lt__(self, other):
        return operator.lt(self, other)

    ## <a href="http://docs.python.org/library/operator.html"> <= </a>
    #  decorated with cmpdec() .
    @_cmpdec
    def __le__(self, other):
        return operator.le(self, other)

    ## <a href="http://docs.python.org/library/operator.html"> > </a>
    #  decorated with cmpdec() .
    @_cmpdec
    def __gt__(self, other):
        return operator.gt(self, other)

    ## <a href="http://docs.python.org/library/operator.html"> >= </a>
    #  decorated with cmpdec() .
    @_cmpdec
    def __ge__(self, other):
        return operator.ge(self, other)

    ## Dot is rejected, though the reason is not fundamental.
    def dot(self, other):
        """s.dot(s') is a TypeError for ascetic reasons."""
        #pylint: disable=W0613,R0201
        raise TypeError('Use "*" for scalar multiplication')

    ## Trivial projection onto basis scalars
    # \yields self.
    def polyadics(self):
        """Needless-to-say: this is trivial."""
        yield self


## The 3D Scalar n Euclidean Space, which for architectural reasons,
# cannot be shared with Minkowski space (I mean, they _are_ different
# objects; a Euclidean 3 is not the same as a relativistic 3, no?).
class Scalar(NScalar):
    """s = Scalar(w) is a rank-0 tensor with one attribute:

    s.w

    which can be a singleton, array_like, or an iterator. You need Scalars
    because theyk now about Vector/Tensor operations, while singletons and
    np.ndarrays do not.


    ZERO
    ONE

    are module constants that are scalars."""

    ## The ranked meta class figures out the indices
    __metaclass__ = euclid.ranked

    ## A scalar times a number is straight dilation, tensor dilations
    # are added to this by the ranked metaclass (this has to be an
    # antipattern)-- again, for matrix components, right and left dilation
    # must be distinguished--as a scalar can represent a scalar field
    # with non-scalar _internal_ degrees of freedom.
    _dispatch_mul = {None: euclid.Tensor_.right_dilation,
                     1j: _dilate_q}

    ## What to do after a sandwiched Grassmann multiplication--while an
    # antipattern under the usual rule of OO programming, those concerns
    # are superseded by Clifford (geometric) algebra: that all 8
    # basis states can be represented by a biquaternion.
    # \param other A quaternion
    # \return scalar part of other.
    @staticmethod
    def _qt(other):
        """Here qt gets the scalar"""
        return other.scalar

    ## Promote to a Versor (un-normalized, see quaternion()).
    # \returns geo.metric.euler.hamilton.Versor
    def versor(self):
        """Promote a scalar to a versor (quaternion) by adding null vector.
        Now there's a problem: it is not normalized, so its not really a
        versor, but If you normalize it, that it's +/-1 are the only
        allowed scalar values.

        This gets used in order to do polymorphic rotations on geometric
        objects....cleanup is TBD
        """
        from .vector import Vector  # create array_like on the spot
        from ..euler.hamilton import Versor
        return Versor(self, Vector(*itertools.repeat(0*self.w, 3)))

    ## \f$ s \rightarrow (s; {\bf\vec{0}}) \f$
    # \returns geo.desic.euler.hamilton.Quaternion
    def quaternion(self):
        """Convert to a quaternion (see versor.__doc__)"""
        return self.versor().quaternion()

    ## \f$ s \rightarrow \frac{1}{3}s \delta_{ij} \f$
    def tensor(self):
        """Upgrade to |0, 0> Cartesian Tensor, with the same Trace"""
        from .tensor import DELTA
        return DELTA * (self.w/3.)

    ## Multiplicative Inverse (is this the best way?).
    # \returns \f$ s^{-1} \f$
    @property
    def I(self):
        """Inverse (mult.)"""
        return type(self)(self.w ** (-1))

    ## Spherical (Fundamental, complex) Representation:
    # \returns \f$T\f$ as a geo.metric.wigner.eckart.Vector
    def spherical(self):
        """Convert from Cartesian Scalar to SO(3) Trivial Scalar."""
        from ..wigner import eckart
        return eckart.Scalar.fromcart(self)

    ## Project Irreducible Representations (Trivially)
    # \returns scalar \f$ \langle 0,0|w|0,0\rangle \f$
    def irrep(self):
        """irrep(j=0, m=0, q=0) --> weight j, seniority index q."""
        return self


## Scalar Unit
ONE = Scalar.basis()[0]


## Scalar Null
ZERO = Scalar.null()


## The Scalar Basis, is for now, NOT iterable-- but that behavior needs to
# be considered.
BASIS = ONE
