"""Gibbs was the 1st to complexify vectors, and since geo's R3 vectors
went out and got complexified: here we are.

Complexification is harder than it looks. The root problem is that inner
product becomes:

<a, b> = a * b.C = conj(<b, a>)

Inner:  (1, 3)(3, 1) --> (1, 1)
         row    col

Outer:  (3, 1)(1, 3) --> (3, 3)
         col    row

Of course, this entire mess could have been avoided by defining the
inner product as:

<a|b> = 1/4 * [(a + b)**2 - (a - b)**2]

and then use the x.norm()**2 method inplace of x**2

"""
## \namespace geo.metric.euclid.gibbs Complexification of
# \f$ \mathbb R^3,\  \mathbb C^3\f$.
from . import vector

#__all__('Gibbs'

## A vector.Vector living in \f$ \mathbb C^3\f$
class Gibbs(vector.Vector):
    """vector = Gibbs(x, y, z)

    is a vector on C3, with inner product:

    v'*v= x'x* + y'y* + z'z*

    Mixing them with R3 vector may cause trouble, as they have different
    dot product."""

    ## Bra prepares an argument for dot/outer product according to the
    # module's inner product:\n
    # \f$ \langle u, v \rangle \equiv {\bf\vec{u}\cdot\vec{v}^*} \f$
    # \param v \f$ v \f$
    # \returns \f$ v \f$
    def bra(self):
        """bra = ket.C"""
        return self.C

    ## Send back to \f$ \mathbb R^3\f$ (complex or not).
    @property
    def vector(self):
        """as a Vector in R3"""
        return vector.Vector(self.x, self.y, self.z)

    ## Override vector.Vector.rejection for complex space.
    # \param self \f$ {\bf \vec{a}}\f$
    # \param other \f$ {\bf \vec{b}}\f$
    # \return \f$ {\bf -(\vec{a}\times\hat{b}^*)\times\hat{b}}\f$
    def rejection(self, other):
        """a.rejection(b) --> -(a ^ b.C) ^ b / ||b||**2"""
        v = other.hat()
        return -(self ^ v.C) ^ v

    ## Spinor map is for Real vectors only
    # \throws HomomorphismError
    def spinor(self, alpha=0):
        """Only real vectors can be mapped to spinors."""
        from ...utils.exceptions import HomomorphismError
        raise HomomorphismError("C3 cannot be mapped to C2, safely")    
    
## Complex Rejection Operator
# \param a \f$ {\bf \vec{a}} \f$
# \returns \f$ -[({\bf \hat{a}^*} \times)][(\bf \hat{a} \times)] \f$
def Rej(a):
    """R_v(a) = [-(v* X)*(v X)]a """  # check signs
    v = a.hat()
    return -v.cross().C * v.cross()


# pylint: disable=W0611
# These functions are complex/real agnostic; it is their arguments that
# implement the correct behavior while on the inside (via dot(), which
# get correct inner product from bra())
from .vector import Proj, Ref
