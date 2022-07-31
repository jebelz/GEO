"""An experimental module for geometric (Clifford) algebra on R3,
implemented with complex combinations of Parity (Even, Odd)
singlet/triplet pairs.

Overloads are:

-------------------
(A & B)           AB   Geometric Product
 A * B                 dot product
(A ^ B)                Wedge Product

~A                     Reverse
~A                     Hodge Dual

abs(A)                 norm

A >> B                 Right contraction
B << A                 Left contraction


There are no explicit blades.


Needless to say, geo should be a package of JUST geometric algebra, but
alas it is not. It's the old fashion vector/tensor motis operandi. Much of the
trouble in reconciling pauli with quaternion with spinors is exactly because
it is not adherent to G3.

So this module is really just for fun.
"""
## \namespace geo.desic.clifford
#  <a href="http://en.wikipedia.org/wiki/Geometric_algebra">Geometric
#  Algebra</a>, __CL(R3)__.

import collections
import operator
import itertools

## Need a Vector and Scalar base
from .. import Vector, Scalar

from ..utils.trig import Re, Im, conj


## Singlet projection function
singlet = operator.attrgetter("singlet")


## Triplet projection function
triplet = operator.attrgetter("triplet")


## See _grade() for implementation
Grader = collections.namedtuple("Grader",
                                "multiplicity complexity factor")

## Hash table to get grade
GRADES = {0: Grader(singlet, Re, 1),
          1: Grader(triplet, Re, 1),
          2: Grader(triplet, Im, 1j),
          3: Grader(singlet, Im, 1j)}


## \param blade A MultiVector
#  \param n     the grade of the blade
#  \returns MultiVector that is purely grade n.
def _grade(blade, n):
    """Get pure grade-n portion of a multi-vector"""
    try:
        # Get method, function, and conversion from hash-table.
        method, func, factor = GRADES[n]
    except KeyError:
        result = None
    else:
        result = factor * func(method(blade))
    return result


## Geometric Algebra's Multi Vector
class MultiVector(object):
    """MultiVector(_singlet, _triplet)

    _singlet is the scalar + i * pseudo-scalar
    _triplet is the vector + i * axial-vector
    """

    components = ('_singlet', '_triplet')

    grades = (0, 1, 2, 3)

    ## Get pure grade component.
    # \param n (grade)
    # \returns MultiVector with pure grade
    def grade(self, n):
        """pure grade-n part multi-vector of multi-vector"""
        cls = type(self)
        if n in (0, 3):  # grade a multivector
            value = _grade(self, n)
            return cls(value, Vector(0, 0, 0))
        elif n in (1, 2):
            value = _grade(self, n)
            return cls(Scalar(0), value)
        else:
            return cls(Scalar(0), Vector(0, 0, 0))

    ## Construct from blades
    @classmethod
    def fromblades(cls, *args):
        """add 4 blades to make a multi-vector"""
        return cls(args[0] + args[3],
                   args[1] + args[2])

    ## Construct from dictionary
    @classmethod
    def fromdict(cls, dict_):
        """make a multi-vector from a dictionary of blades"""
        return cls(dict_[0] + 1j * dict_[3],
                   dict_[1] + 1j * dict_[2])

    ## Complex singlet and triplet (2 x (1 + 3)) = 8 dimensions
    #  \param singlet_ scalar + (1j) pseudo scalar
    #  \param triplet_ vector + (1j) axial vector
    def __init__(self, singlet_, triplet_):
        ## \f$ s = \rho + i \pi \f$ \n
        #  scalar + \f$ i \f$ pseudo-scalar (trivector)
        self._singlet = singlet_
        ## \f$ \vec{T} = \vec{V} + i \vec{A} \f$ \n
        #   vector + \f$ i \f$ axial vector  (bivector)
        self._triplet = triplet_

    ## Complex representation of grade 0, grade 3.
    @property
    def singlet(self):
        """singlet"""
        return self._singlet

    ## Complex representation of grade 1, grade 2.
    @property
    def triplet(self):
        """triplet"""
        return self._triplet

    ## iter(multivector) --> singlet, triplet
    def __iter__(self):
        return iter((self.singlet, self.triplet))

    ## {singlet: triplet}
    def __str__(self):
        return "{%s ; %s}" % (str(self.singlet.w), str(self.triplet))

    ## 4
    # \returns 4 \f$ 4 \f$
    def __len__(self):
        return 4

    ## grade().
    # \param grade grade \f$s\f$
    # \returns grade   \f$ \langle A \rangle_s \f$
    def __getitem__(self, grade):
        return self.grade(grade)

    ## any non-zero stuff
    # \returns bool If non-trivial
    def __nonzero__(self):
        for item in self:
            if item:  # short-curcit any implementation
                return True
        return False

    ## Negation
    # \returns multivector Additive inverse: \f$ {\bf -a} \f$
    def __neg__(self):
        return type(self)(*map(operator.neg, self))

    ## Addition
    # \param self (implicit) \f$ {\bf a} \f$
    # \param other (explicit) \f$ {\bf b} \f$
    # \returns multivector Component-wise addition \f$ {\bf a+b} \f$
    def __add__(self, other):
        return type(self)(*[operator.add(a, b) for a, b in zip(self, other)])

    ## Subtraction
    # \param self (implicit) \f$ {\bf a} \f$
    # \param other (explicit) \f$ {\bf b} \f$
    # \returns multivector Component-wise subtraction  \f$ {\bf a-b} \f$
    def __sub__(self, other):
        return self + (-other)

    ## Antisymmetric Product (cross like) \n
    # via _mix() with operator.add
    # \param self (implicit) \f$ {\bf a} \f$
    # \param other (explicit) \f$ {\bf b} \f$
    # \returns multivector wedge product: \f$ {\bf a \vee b} \f$
    def __xor__(self, other):
        return reduce(operator.add,
                      self._mix(other, operator.add),
                      type(self)(Scalar(0), Vector(0, 0, 0)))

    ## Symmetric Product (dot like) \n
    # via _mix() with operator.sub
    # \param self (implicit) \f$ {\bf a} \f$
    # \param other (explicit) \f$ {\bf b} \f$
    # \returns multivector  \f$ {\bf a \cdot b} \f$
    def __mul__(self, other):
        return reduce(operator.add,
                      self._mix(other, operator.sub),
                      type(self)(Scalar(0), Vector(0, 0, 0)))

    ## Private mixer wrt a bilinear operator
    # \param other MultiVector
    # \param func a binary operation
    # \returns list_ a list of wedge-products of grades
    # \f$ f(u, v; o) = \sum_{i, j \forall i, j \in (0, 1, 2, 3) |
    # o(i, j)}{(v_iu_j)_{o(i, j)}} \f$
    def _mix(self, other, func):
        """mix wrt to func, according to clifford algebra"""
        return [(self[i] & other[j])[func(i, j)]
                for i in (0, 1, 2, 3)
                for j in (0, 1, 2, 3) if func(i, j) in self.grades and
                self[i] and other[j]]

    ## Involution ala Cayley-Dickson: \n
    # \returns multivector \f$ (S + \vec{T}) \rightarrow (S^ * - \vec{T}^*) \f$
    @property
    def involution(self):
        """ 0, 1, 2, 3 --> 1, -1, 1, -1"""
        return type(self)(conj(self.singlet), -conj(self.triplet))

    ## __THE__ Geometric Product
    # \param self (implicit) \f$ {\bf a} = (s, t) \f$
    # \param other (explicit) \f$ {\bf a'} = (s', t') \f$
    # \returns multivector
    # \f$ (ss'+\vec{T}\cdot\vec{T}', s\vec{T}'+s'\vec{T}+i \vec{T}
    # \times i\vec{T}') \f$ \n
    # It is  bound to "&" (a&b) since you can't interpret "ab" as such.
    def __and__(self, other):
        """A&B is the geometric product"""
        return type(self)(
            (self.singlet * other.singlet) +
            (self.triplet * other.triplet),
            (self.singlet * other.triplet +
             self.triplet * other.singlet) +
            ((self.triplet) ^ 1j*(other.triplet))
            )

    ## Dilation (Inversion)
    #  \param a A number
    #  \returns multivector \f$ {\bf a}/a \f$
    def __div__(self, a):
        return type(self)(*[item/a for item in self])

    ## Dilation
    #  \param a A number
    #  \returns multivector \f$ {\bf a}/a^{-1} \f$
    def __rmul__(self, a):
        return self/(1./a)

    ## Cross Product
    #  \param self (implicit) \f$ {\bf a} \f$
    #  \param other (explicit) \f$ {\bf b} \f$
    #  \returns multivector
    #  \f$ (-{\bf \Pi}A)\wedge B \rightarrow \vec {A} \times \vec{B} \f$
    def cross(self, other):
        """cross productor for bi-vectors"""
        return (-PI & self) ^ other

    __repr__ = __str__

    def __eq__(self, other):
        return all(left == right for left, right in zip(self, other))

    def __ne__(self, other):
        return not self == other

    ## \returns multivector
    @property
    def real(self):
        """real part"""
        return type(self)(*map(Re, self))

    ## \returns multivector
    @property
    def imag(self):
        """imaginary (so-called) part"""
        return type(self)(*map(Im, self))

    ## Conjugate
    # \returns MultiVector: \f$ (\pi, \vec{A}, \vec{V}, 1) \f$
    def conjugate(self):
        """conjugation"""
        return type(self)(*map(conj, self))

    ## Dual
    #  \returns multivector  \f$ {\bf \tilde{a}} \rightarrow {\bf a \Pi} \f$
    def dual(self):
        """dual"""
        return self & (-PI)

    ## Note: this is Not a MultiVector- which contrasts with geo's
    #  implementation.
    def __abs__(self):
        return abs(
            (reduce(
                operator.add,
                [(l*r).w for l, r in itertools.izip(self,
                                                    itertools.imap(conj,
                                                                   self))])))

    ## Inversion:
    # \returns \f$ {\bf a}^{-1}\rightarrow{\bf \tilde{a}}/|{\bf \tilde{a}}|\f$
    def __invert__(self):
        return self.dual() / abs(self)**2

    ## A << B \n (maybe)
    # \param self (implicit) \f$ {\bf a} \f$
    # \param other (explicit) \f$ {\bf b} \f$
    # \returns \f$ -\tilde{(a \wedge \tilde{b})} \f$
    def __lshift__(self, other):
        return -(self ^ (-other.dual())).dual()

    ## A >> B \n
    # \param self (implicit) \f$ {\bf a} \f$
    # \param other (explicit) \f$ {\bf b} \f$
    # \return \f$ -\tilde{(\tilde{a} \wedge b)} \f$
    def __rshift__(self, other):
        return -(-self.dual() ^ other).dual()


## The Pseudoscalar (tri-vector) or \f$ R^3 \f$.
PI = MultiVector(Scalar(1j), Vector(0, 0, 0))

## 8 Dimensional Basis:\n
# \f$ (1, e_1, e_2, e_3, e_1e_2, e_2e_3, e_3e_1, \Pi) \f$
BASIS = collections.namedtuple("MultiVectorBasis",
                               "rho x y z i j k pi")(
                                   *itertools.starmap(
                                       MultiVector,
                                       ((Scalar(1), Vector(0, 0, 0)),
                                        (Scalar(0), Vector(1, 0, 0)),
                                        (Scalar(0), Vector(0, 1, 0)),
                                        (Scalar(0), Vector(0, 0, 1)),
                                        (Scalar(0), Vector(1j, 0, 0)),
                                        (Scalar(0), Vector(0, 1j, 0)),
                                        (Scalar(0), Vector(0, 0, 1j)),
                                        PI)))

'''
## Practice using it...
if __name__ == '__main__':
    print "BASIS: \n", "\n".join(map(str, BASIS))
    ## Named Basis Elements
    ONE, X, Y, Z, I, J, K, PI = BASIS
    e1, e2, e3 = X, Y, Z
    e12 = K
    e23 = I
    e31 = J

    e123 = PI
    c = 0.5**0.5

    U = c*(X+Y)
    V = c*(X-Y)

    v = 3.*X + 4.*Y

    q = c*PI + c*K
    qbar = c*PI - c*K

    p = c*ONE + c*K
    pbar = c*ONE - c*K
'''
