"""Minkowski Space. This module defines:

FourScalar
FourVector

with the standard metric:

+1  0   0   0
0  -1   0   0
0   0  -1   0
0   0   0  -1
"""
#pylint: disable=E1101

## \namespace geo.desic.minkowski.minkowski
#  <a href="http://en.wikipedia.org/wiki/Four-vector">4-Vectors</a> and
#  4-Scalars

from ...metric.euclid import vector
from ...metric.euclid import scalar
from ...metric.euler import hamilton as versor

#from .. import hamilton as quaternion


## Helper
# \param t time-like component
# \param x space-like x
# \param y space-like y
# \param z space-like z
# \retval FourVector \f$ (t, x\hat{x}+ y\hat{y}+ z\hat{z}) \f$
def txyz2vector(t, x, y, z):
    """t, x, y, z, --> Four Vector"""
    ## Todo check If this can be an inherited classmethod?
    return FourVector(scalar.Scalar(t), vector.Vector(x, y, z))


## <a href="http://en.wikipedia.org/wiki/Scalar_(physics)#Scalars_in_relativity_theory">Four Scalars</a>.
class FourScalar(scalar.NScalar):
    """A Poincare invariant scalar:

    FourScalar(w)

    w array_like number.
    """
    ## prevent rank interaction (?) Need to decide what the ZOO is.
    __metaclass__ = type

    ## Four scalars have a 'w'.
    components = ('w',)


## <a href="http://en.wikipedia.org/wiki/Four-vector">Four Vector</a>.
class FourVector(versor.ExE3, quaternion.FourOpsMixin):
    """FourVector(scalar, vector)

    scalar: a 3-scalar
    vector: a 3-vector
    """
    components = ("_scalar", "_vector")

    ## Why do different? b/c classmethods don't call super; TODO: find
    # and abc trick to nicen this code.
    fromtxyz = classmethod(lambda cls, *a: cls.from_primitives(*a))

    @classmethod
    def fromspinor(cls, m):
        """Not Implemented"""
        raise RuntimeError("not V & V'd'")

    @classmethod
    def fromarray(cls, array_):
        """Straight unpacking of t, x, y, z"""
        return cls.fromtxyz(*array_)

    ## Time like component
    @property
    def time(self):
        """time-like component"""
        return self.scalar

    ## Space like component
    @property
    def space(self):
        """space-like component"""
        return self.vector

    ## t is time, and it is a 3-scalar
    t = time

    ## s is space, and it is a 3-vector
    s = space

    ## Minkowski space is pseudo normed, does any physics need this?
    def __abs__(self):
        try:
            return super(FourVector, self).__abs__()
        except ValueError:
            pass
        # this is ugly code: (1) Import deep into the method
        from cmath import sqrt
        # (2) redo super's code.
        normed2 = self.dot(self)
        # force recast as a FourScalar
        return type(normed2)(sqrt(normed2))

    ## We need this, because super will say 'w=' for the time component,
    # and that is just confusing.
    def __str__(self):
        return "({}; {})".format(self.time.w, self.space)

    ## 4-dot product:\n
    #  \f$ \rho = v_{\mu}u_{\mu} \equiv g_{\mu\nu}v_{\mu}u_{\nu} \f$
    def dot(self, other):
        """minkowski dot product"""
        return FourScalar(self.t * other.t - self.s.dot(other.s))

    ## "*" overload is dot product
    __mul__ = dot

    ## \f$ v^2 = v_{\mu}v_{\mu} \f$
    def __pow__(self, r):
        if r == 2:  # raise / power branch
            return self*self
        raise ValueError("2 or don't do it")

    ## 4-outer product:\n
    # \f$ \Lambda_{\mu\nu} = v'_{\mu}v_{\nu} \f$
    # \param self A FourVector (implicit, of course)
    # \param other A FourVector (explicit, of course)
    # \retval lorentz A lorentz.Lorentz object
    def __and__(self, other):
        from .lorentz import Lorentz
        return Lorentz(self.t * other.t,
                       self.t * other.s,
                       self.s * other.t,
                       self.s & other.s)

    ## \f$ v'_{\mu} =\Lambda_{\mu\nu}v_{\nu} \f$
    # \param self A FourVector object
    # \param lam A lorentz.Lorentz object
    # \retval lam*self Another FourVector object.
    def __rmul__(self, lam):
        """
        |t'|    |TT    TS| |t|    |TT*t + TS*s|
        |  | =  |        | | | =  |           |
        |s'|    |ST    SS| |s|    |ST*t + SS*s|
        """
        return type(self)(lam.tt*self.t + lam.ts*self.s,
                          lam.st*self.t + lam.ss*self.s)


## The Poincare Unit.
ONE = FourScalar(1)


def test():
    """test"""
    v = txyz2vector(3., 4., 12., -13.)
    return v
