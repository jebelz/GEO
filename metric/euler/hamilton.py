"""The versor module: These are unit quaternions, for doing rotations in R3.

Needless-to-say: the circle-ellipse problem is homomorphic to the
versor-quaternion problem.

These are *JUST* unit quaternions. You can't add/subtract them-- as they
would put them off-shell. They are JUST rotation reps:

Versor(scalar, vector) -- not that you would use the constructor.

You use:

>>>r = roll(theta)
>>>p = pitch(theta)
>>>y = yaw(theta)

and their products:

>>>yaw_pitch_and_then_roll = r * p * y

or something else from charts.py, or:

>>>Vector.versor(angle [,circumference=2pi])


>>>xyzw2versor(x, y, z, w)
>>>wxyz2versor(w, x, y, z)

as needed.

If you need to read them in from a file, then use:

Versor.fromfile(src)  #wxyz
Versor.fromelif(src)  #xyzw.
"""
## \namespace geo.metric.euler.hamilton
#  <a href="http://en.wikipedia.org/wiki/Versor">Versors</a>, or unit
#  quaternions.
import abc
import itertools
import math
import operator

from ...utils import arraylike #import ABCNumpy, enlister

from ..euclid import scalar
from ..euclid import vector
from ..euclid import euclid

#pylint: disable=R0904

__all__ = ('slerp', 'Versor', 'W', 'I', 'J', 'K', 'hopf')


##  <a href="http://en.wikipedia.org/wiki/Slerp">Spherical Linear Interpolation
#  </a>
# \param p0 Starting Rotation: \f$ p_0 \f$
# \param p1 Ending Rotation: \f$ p_1 \f$
# \param alpha Parameterization of geodesic: \f$ \alpha \f$
# \returns Versor \f$ q = p_0 [\bar{p}_0 p_1]^{\alpha} \f$
# @image html slerp.jpeg
# @image latex slerp.jpeg
def slerp(p0, p1, alpha):
    """q = slerp(p1, p2; alpha)

    Spherical Linear Interpolation:

    Start at:
    p0 (alpha=0) and move along hypersphere to
    p1 (alpha=1).

    q =  p0 * (p0 >> p1) ** alpha
    """
    # now that is operator overloading beauty
    return p0 * (p0 >> p1) ** alpha


## Direct matrix conversion (experimental-and more robust against
# numpy's array hi-jinx, but a complete fail wrt to geo's philosophy)
# \param w Versor.w
# \param u Versor.i
# \param j Versor.j
# \param k Versor.k
# \returns \f$ {\bf T} = 2\left[ \begin{array}{ccc}
# \frac{1}{2} - (j^2 + k^2) & ij+wk & ik-wj \\
# ij-wk & \frac{1}{2} -(i^2+k^2) & jk+wi \\
# ik + wi & jk - wi & \frac{1}{2} - (i^2+j^2) \end{array} \right] \f$
def hopf(w, i, j, k):
    """dcm = jopf(w, i, j, k)
    Tensor from a Versor's components'."""
    from ...metric.euclid import tensor
    return tensor.Tensor(
        1 - 2*(j**2 + k**2), 2*(i*j + w*k), 2*(i*k - w*j),
        2*(i*j - w*k), 1 - 2*(i**2 + k**2), 2*(j*k + w*i),
        2*(i*k + w*j), 2*(j*k - w*i), 1 - 2*(i**2 + j**2))


## A class to compose a Scalar and a Vector--for who knows what?
class ExE3(arraylike.ABCNumpy):
    """Composer Mixin for arraylike classes that are composed
    of a Scalar and a Vector."""

    __metaclass__ = abc.ABCMeta

    ## At least: \f$ {\bf R}_0 \oplus \bigotimes_{i=1}^{i=3}{\bf R}_i \f$
    components = ('_scalar', '_vector')

    ## Helper from components
    # \param w real component
    # \param x i-component
    # \param y j-component
    # \param z k-component
    # \returns Versor
    @classmethod
    def from_primitives(cls, w, x, y, z):
        """(w, x, y, z)-->Versor(w, (x, y, z))"""
        from ..euclid.scalar import Scalar
        from ..euclid.vector import Vector
        return cls(Scalar(w), Vector(x, y, z))

    ## Construct from 2 x 2 complex matrix
    # \param m from SU(2).
    # \returns Versor
    @classmethod
    def fromspinor(cls, m):
        """vector = Vector.fromspinor(m)


        m is a 2 x 2 complex matrix
        """
        w = (m[0, 0] + m[1, 1]) / 2.
        z = (-m[0, 0] + m[1, 1]) / 2j

        x = -(m[0, 1] + m[1, 0]) / 2j
        y = (-m[0, 1] + m[1, 0]) / 2.
        #pylint: disable=E1101
        return cls.fromwxyz(w, x, y, z)

    ## There is no map from a generic matrix to a Quaternion.
    # \throws ValueError
    @staticmethod
    def frommatrix(m):
        """This operation is not defined."""
        raise ValueError('from: {}? Try fromarray or fromspinor.'.format(m))

    ## \f$ {\bf q} \equiv (q; \vec{q}) \f$ \n Takes a euclid.Scalar and a
    #  euclid.Vector
    def __init__(self, scalar_, vector_):
        """Versor(scalar, vector):

        scalar --> sin(theta/2) as a Scalar instance
        vector --> cos(theta/2)*unit_vector as a Vector instance.

        Likewise, you can pull out:
        (w, x, y, z) If needed.
        """
        ## euclid.Scalar part
        self._scalar = scalar_
        ## euclid.Vector part
        self._vector = vector_

    ##   "{" + "; ".join(map(str, self.iter())) + "}"
    def __str__(self):
        return "{%s; %s}" % tuple(map(str, self.iter()))

    ## Inverse of eval.
    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__,
                                   *map(repr, self.iter()))

    ## Straight component negation.
    def __neg__(self):
        return type(self)(*itertools.imap(operator.neg, self.iter()))

    ## TBD hooks.
    @property
    def scalar(self):
        """scalar getter"""
        return self._scalar

    ## TODO: unify .scalar or .scalar() for scalar part
    @scalar.setter
    def scalar(self, scalar_):
        """scalar setter"""
        self._scalar = scalar_

    ## TODO same for vector.
    @property
    def vector(self):
        """vector getter"""
        return self._vector

    @vector.setter
    def vector(self, vector_):
        """vector setter"""
        self._vector = vector_

    ## Vector's 'x'-component.
    @property
    def x(self):
        """i-->vector.x"""
        return self.vector.x

    ## Vector's 'y'-component.
    @property
    def y(self):
        """j-->vector.y"""
        return self.vector.y

    ## Vector's 'z'-component.
    @property
    def z(self):
        """k-->vector.z"""
        return self.vector.z

    ## Norm: definition of metric is not here (see: dot() ).
    # \f$ ||{\bf q}|| \equiv \sqrt{{\bf{q \cdot q}}} \f$
    def __abs__(self):
        """||q|| = q.dot(q)**(1/2)"""
        #pylint: disable=E1101
        return (self.dot(self))**0.5


## General Quaternion Mixin
class QuaternionMixin(ExE3, euclid.LinearMap):
    """This mixin has ops for Versors and Quaternions."""

    __metaclass__ = abc.ABCMeta

    ## External polymorphism masquerading as math
    rank = 1j

    ## Why so different? b/c minkowski space hath no double ewe.
    fromwxyz = classmethod(lambda cls, *a: cls.from_primitives(*a))

    ## _IBID_
    fromxyzw = classmethod(lambda cls, *a: cls.from_primitives(a[-1], *a[:-1]))

    ## There is no unique map from a generic matrix to a Quaternion.
    # \throws ValueError
    @staticmethod
    def fromarray(a):
        """Illdefined Problem."""
        raise ValueError('from: {}? try fromwxyz or fromxyzw.'.format(a))

    ## \f$ ||{\bf q\cdot p}|| \equiv qp + \vec{q}\cdot\vec{p} \f$ \n
    # The Quaternion dot product
    def dot(self, other):
        """dot(q, q') = --> s*s' + v*v'"""
        return self.scalar * other.scalar + self.vector * other.vector

    ## \f$ {\bf \bar{q}} \rightarrow (q; -\vec{q}) \f$
    def conjugate(self):
        """(s; v).conjugate()--> (s; -v)"""
        return type(self)(self.scalar, -self.vector)

    C = property(conjugate)

    @property
    def T(self):
        """The Transpose in trivial"""
        return self

    @property
    def H(self):
        """H-->Hermitian Conjugate"""
        return self.T.C  # todo: pushup.

    ## grassmann() product (w/o type checking)
    def __mul__(self, versor_):
        """Grassmann product:

        qq' = (s; v)(s'; v')

        ... = (s*s' - v*v'; sv' - s'v + v X v')
        """
        return self.grassmann(versor_)

    ## q2 = q3/q1 -> q3(v) == q2(q1(v))
    # \throw Exception This is not well thought out.
    def __div__(self, other):
        """return (~other) * self"""
        raise Exception('ORDER error')

    ## SLERP (slerp()) from W or other to self.
    # \param alpha An exponent
    # \param other (optional) starting point will invoke slerp.
    # \return \f$ q^{\alpha} = (\cos{\frac{\theta}{2}},
    # \sin{\frac{\theta}{2}} {\bf \hat{n}})^{\alpha} =
    # (\cos{\frac{\alpha\theta}{2}},
    # \sin{\frac{\alpha\theta}{2}}  {\bf \hat{n}}) \f$
    def __pow__(self, alpha, other=None):
        """q' = pow(self, alpha, other=W):

        == slerp(other, self, alpha). Thus, traveling along the unit
        hypersphere:

        [p]--------[q']----------------------[q]-->
         0........<alpha>.....................1......

        q' is 100*alpha percent along the path from p to q.

        With:

        >>>self ** alpha  # other defaults to 'W', the indentity,

        is the standard 2-argument power function.
        see slerp.__doc__."""
        from ...utils.trig import arctan2, sin, cos
        if other is not None:  # kwd
            return slerp(self, other, alpha)
        f = abs(self).w ** alpha
        sinth = abs(self.vector).w  # but how?
        theta = arctan2(sinth, self.w)
        rat = sin(alpha*theta)/(sinth)
        return type(self)(scalar.Scalar(cos(alpha*theta)) * f,
                          self.vector * rat * f)

    ## str(q) --> {str(s); str(v)}
    def __str__(self):
        return "{"+str(self.w)+"; "+str(self.vector)+"}"

    ## Grassmann (antisymmetric) product
    # Is the antisymmetric product on \f$ {\bf H} \f$.
    # \param other A Versor, Quaternion, Vector, or Scalar, Tensor, YPR, etc.
    # \returns \f$ {\bf q}{\bf p} = (q; \vec{q})(p; \vec{p}) =
    # (qp-\vec{q}\cdot\vec{p}; q\vec{p} + p\vec{q}+\vec{q}\times\vec{p})\f$
    def grassmann(self, other):
        """Grassmann product with ANOTHER versor or quaternion"""
        try:
            return (type(self) | type(other))(
                self.scalar * other.scalar - self.vector * other.vector,
                self.scalar * other.vector + self.vector * other.scalar +
                (self.vector ^ other.vector)
            )
        except AttributeError:
            return self.grassmann(other.versor())  # external poly, -1 point.

    ## What to do after a sandwiched Grassmann multiplication
    @staticmethod
    def _qt(other):
        """_qt(versor) --> versor, again compare for polymorphism"""
        return other

    ## \throws Versor.DivisionError
    def __rdiv__(self, other):
        raise self.DivisionError

    ## \returns geo.utils.exceptions.OverloadNotImplementedError
    @property
    def DivisionError(self):
        """figure out the Division Error to throw"""
        from ...utils.exceptions import OverloadNotImplementedError
        return OverloadNotImplementedError(
            "/quaternion => * ~quaternion overload not implemented"
            )

    ## Inverse is inv()
    # \return \f$ {\bf \tilde{q}} = {\bf \bar{q}} \f$
    def __invert__(self):
        """Inverse."""
        return self.inv()

    ## Implicit inverse
    # \returns \f$ {\bf \tilde{q}} = {\bf \bar{q}}/||{\bf q}||^2 \f$
    def inv(self):
        """inverse assumes unit norm"""
        return self.C

    ## \returns inv()
    @property
    def I(self):
        """I --> inv()"""
        return self.inv()

    ## To a 4 component list
    def wxyz(self):
        """q --> [w, x, y, z]"""
        return self.scalar.tolist() + self.vector.tolist()

    ## To a 4 component list
    def xyzw(self):
        """q --> [x, y, z, w]"""
        return self.vector.tolist() + self.scalar.tolist()

    ## Forwarding property
    @property
    def w(self):
        """w-->scalar.w"""
        return self.scalar.w

    ## Non Clifford identity
    i = ExE3.x
    j = ExE3.y
    k = ExE3.z

    ## TODO: fix append to be more numpy like
    def append(self, other):
        return type(self)(self.scalar.append(other.scalar),
                          self.vector.append(other.vector))

    ## Todo: likewise keyword on rank 2 +
    # \kwd extend = True for extended list (versus appended)
    # \returns list of components or list of lists of components.
    def tolist(self, extend=True):
        """tolist(extend=True):

        [w, x, y, z] If extend else:
        [list(scalar), list(vector)]
        """
        result = super(QuaternionMixin, self).tolist()
        if extend:  # kwd
            result = sum([item.tolist() for item in result], [])
        return result

    ## Rep as SO(4).
    # \returns 4x4 matrix, or list of them.
    @arraylike.enlister
    def fourXfour(self):
        """Convert to 4 x 4 real matrix rep."""
        from . import van_elfrinkhof
        return sum(a*b for a, b in zip(self.tolist(extend=True),
                                       van_elfrinkhof.BASIS))

    @classmethod
    def _from4x4(cls, so4):
        from . import van_elfrinkhof
        return cls.fromwxyz(van_elfrinkhof.project_S(so4),
                            van_elfrinkhof.project_X(so4),
                            van_elfrinkhof.project_Y(so4),
                            van_elfrinkhof.project_Z(so4))

    ## Construct from SO(4) rep.
    # \param *args Any number of real 4 x 4 matrices.
    # \returns cls instance.
    @classmethod
    def from4x4(cls, *args):
        """Construct from 4 x 4 real matrices."""
        for count, value in enumerate(map(cls._from4x4, args)):
            if not count:  # init
                result = value
            else:
                result = result.append(value)
        return result

    ## Angle between 2 rotations.
    # \param other Versor
    # \returns \f$ \theta = \cos^{-1}{2({\bf +q \cdot +p})^2 -1 } \f$
    def __or__(self, other):
        """needs work.. pos is irrelavent b/c **2"""
        from ...utils.trig import arccos
        return arccos(2*(+(self).dot(+other)).w**2 - 1)

    ## R = Q >> P take Q to P:  P = QR
    # \param \f$ a \f$ (implicit as self)
    # \param \f$ b \f$ (explicit as other)
    # \returns \f$ c = \bar{a}b \f$ solves \f$ a = bc \f$
    def __rshift__(self, other):
        """Right Divisor (of argument!)

        If:
        R = Q >> P

        then R takes Q to P:

        P == Q * R

        Nuemonic:
        Q * (Q >> P) == P"""
        return (~self) * other

    ## Left Divisor:
    # \param \f$ a \f$ (implicit as self)
    # \param \f$ b \f$ (explicit as other)
    # \returns \f$ c = \bar{a}b \f$ solves \f$ a = cb \f$
    def __lshift__(self, other):
        """Left Divisor (of argument!)

        If:
        Q = B << A

        Then:
        A = Q * B

        Nuemonic:
        A == (A << B) * B"""
        return (~self) >> (~other)

    ## Ensure \f$ w > 0 \f$
    # \returns \f$ \mathop{\mathrm{sgn}}{q}{\bf q} \f$
    def __pos__(self):
        from functools import partial
        from operator import mul
        return self.fromwxyz(*map(partial(mul, self.scalar.sgn()),
                                  self.tolist()))


    ## \f$ {\bf q}(\vec{v}) \rightarrow \vec{v}' \f$ with \n
    #  \f$ (0, \vec{v}') = {\bf \tilde{q}}(0; \vec{v}){\bf q} \f$ \n
    #  is an alias transformation by similarity transform using grassmann()
    #  multiplication (of the versor inverse)\n
    #  __call__ calls this If Alias is a base class
    def alibi_transform(self, other):
        """Alibi transformation"""
        try:
            f = other._qt  # If other has this, it can be sandiwched.
        except AttributeError:
            return self.tensor().alibi_transform(other)  # catch-all deal.
        else:
            return f(self * other.quaternion() * self.C) # ORDER

    ## \returns \f$ \sqrt(q) \equiv q^{\frac{1}{2}} \f$
    def sqrt(self):
        """q' = sqrt(q)"""
        return self ** 0.5

    ## This should always == dcm
    def tensor(self):
        """Convert to a DCM."""
        from ..euclid.tensor import DELTA
        s, v = self.tolist(False)
        return (DELTA * (2*s**2-1) +
                2 * ((v & v) - v.dual()*s))

    ## Of course, a DCM is not a tensor, so maybe you don't want to call
    # it that.
    matrix = tensor

    ## another alias.
    dcm = matrix


## Possible Misnomer: The metaclass allows Versor and Quaternion inter-type
# operations to chose the correct target class.
class Lebesgue(abc.ABCMeta):
    """Experimental pattern, probably an ANTIPATTERN-- but you gotta try
    new things."""

    ## A dictionary of classes made by this metaclass
    cast = {}

    ## Make class, and memoize it Lebesque.cast
    # \sideeffect Updates Lebesque.cast
    def __new__(mcs, *args, **kwargs):
        cls = abc.ABCMeta.__new__(mcs, *args, **kwargs)
        mcs.cast[cls.measure] = cls
        return cls

    ## Bitwise or of cls.measure choses correct space (class) for the
    # resulting operation:
    # \param other A Versor, Quaternion, Vector, or Scalar.
    # \returns class from ::Lebesgue.cast
    # \throws TypeError If other doesn't work out
    # \sideeffect May update ::Lebesque.cast via geo.desic.hamilton import
    def __or__(cls, other):
        """Chose target class:
        V V -> V
        V Q -> Q
        Q V -> Q
        Q Q -> Q

        hence: or.
        """
        key = cls.measure | other.measure
        result = Lebesgue.cast.get(key)
        # See If look up failed
        if result is None:  # external polymorphism
            # If key is 1, then we need to load quaternions
            if key == 1:  # raise
                from geo.desic import hamilton
                result = Lebesgue.cast.get(key)
            else:
                raise TypeError(
                    "Could not cast {} for Quaternion multiplication.".format(
                        type(other)))
            if key != 1:  # raise
                raise TypeError(
                    "Could not cast {} for Quaternion multiplication.".format(
                        type(other)))
            else:
                from geo.desic import hamilton
                result = Lebesgue.cast.get(key)
        return result
            
    ## Todo: split BiQuaternions off from Quaternions.
    def __and__(cls, other):
        """Choose from 3 options:

          V  Q  B
         --------
        V|V  Q  B
        Q|Q  Q  B
        B|B  B  B
        """
        return Lebesgue.case[max(cls.measure, other.measure)]


## Limited <a href="http://en.wikipedia.org/wiki/Versor">Versor</a>
#  class for alias transformations.
class Versor(QuaternionMixin):
    """Versors are unit quaternions. They represent rotations. Alias rotations,
    that is rotation of coordinates, not of vectors.

    You can't add them, you can't divide them. You can:

    *   --> Grassmann product
    ~   --> conjugate (inverse)
    ()  --> transform a vector argument to a representation in a new frame
    q**n  --> spherical linear interpolation (slerp)

    See __init__ for signature %s

    You can get components as:
    w, x, y, z, i, j, k, scalar, vector, roll, pitch, yaw

    You can get equivalent rotation matrices:

    q.AlibiMatrix()
    q.AliasMatrix()
    q.Matrix()    (this pick the correct one from above)

    Or tait bryan angles:

    YPR()
    """

    __metaclass__ = Lebesgue

    ## This is a sub-space of H on the unit hyper-sphere
    measure = 0

    def versor(self):
        """versor is self"""
        return self

    ##\f${\bf q}\rightarrow M=(2q^2-1)I+2(q\vec{q}\times+2\vec{q}\vec{q})\f$
    def AlibiMatrix(self):
        """equivalent matrix for alibi rotation"""
        return self.AliasMatrix().T

    ##\f${\bf q}\rightarrow M=
    #  [(2q^2-1)I+2(q\vec{q}\times+2\vec{q}\vec{q})]^T\f$
    def AliasMatrix(self):
        """equivalent matrix for alias rotation"""
        return hopf(*self.tolist(extend=True))

    ## Default for mul is alibi.
    dcm = AlibiMatrix

    ## AliasMatrix()'s yaw
    @property
    def yaw(self):
        """Yaw angle (YPR ordering)"""
        return self.dcm().yaw

    ## AliasMatrix()'s pitch
    @property
    def pitch(self):
        """Pitch angle (YPR ordering)"""
        return self.dcm().pitch

    ## AliasMatrix()'s roll
    @property
    def roll(self):
        """Roll angle (YPR ordering)"""
        return self.dcm().roll

    ## as a YPR instance
    def ypr(self):
        """YPR instance equivalent"""
        return self.tensor().ypr()

    ## as a RPY instance
    def rpy(self):
        """RPY instance equivalent"""
        return self.ypr().rpy()

    ## A triplet of (x, y, z) in the rotated frame.
    def new_basis(self):
        """map(self, (x, y, z))"""
        print "new_basis is depreavgted- frame"
        return self.frame()

    ## Compute the look angles by transforming the bore-site and getting is
    #  Vector.Polar polar (elevation) and azimuth angle.
    def look_angle(self, boresite=vector.Z):
        """q.look_angle([boresite=vector.Z])

        get a euclid.LookAngle tuple."""
        polar_ = self(boresite)
        return polar_.theta, polar_.phi

    ## \f$ q^2 - \vec{q}^2  = \cos{\theta} \f$
    def angle_of_rotation(self, other=None):
        """q.angle() is the magnitude of the rotation"""
        from ...utils import trig
        if other is None:  # kwd
            value = self.scalar**2 - self.vector**2
        else:
            raise DeprecationWarning
        return trig.arccos(value.w)

    ## Degree version of angle().
    # \param other
    # \returns angle (degrees) of rotatuib
    def dangle_of_rotation(self, other=None):
        """Degrees angle of rotation."""
        return 180 * self.angle_of_rotation(other=other) / math.pi

    ## \returns Vector Unit vector of rotation axis
    def axis_of_rotation(self):
        """Axis of rotation."""
        return self.vector.hat()

    ## Change to a quaternion with LOD violation: but Versors and
    # Quaternion are intimate, so it's OK.
    def quaternion(self):
        """Change to a quaternion"""
        try:
            return self.__metaclass__.cast[1](self.scalar, self.vector)
        except KeyError:
            ## whoops, this is a work around
            from ...desic import hamilton
            return hamilton.Quaternion(self.scalar, self.vector)

    ## \returns Quaternion \f$[a, b] = ab - ba \f$
    def commutator(self, other):
        """[a, b] = ab - ba"""
        return self.quaternion().commutator(other)

    ## Non-scalar cast trace:
    # \f$ 1 + 2 \cos{\theta} \f$
    def trace(self):
        """Trace."""
        from ...utils.trig import cos
        return 1 + 2 * cos(self.angle())

    def euler(self):
        """As Euler angles."""
        return self.tensor().euler()

    def spinor(self):
        """see Quaternion.spinor()"""
        return self.quaternion().spinor()

    ## Construct from a numpy binary file (w, x, y, z)
    # \param args Arguments to np.fromfile
    # \kwd kwargs Keyword arguments to np.fromfile
    # \returns Versor (one versor)
    @classmethod
    def fromfile(cls, *args, **kwargs):
        """fromfile(*args, **kwargs) --> np.fromfile (with reshaping)"""
        from numpy import fromfile
        return cls.fromwxyz(*fromfile(*args, **kwargs).reshape(-1, 4))

    ## Contruct from a binary file with order x, y, z, w.
    # \param args Arguments to np.fromfile
    # \kwd kwargs Keyword arguments to np.fromfile
    @classmethod
    def fromelif(cls, *args, **kwargs):
        """See fromfile.__doc__. This take x, y, z, w order."""
        from numpy import fromfile
        return cls.fromxyzw(*fromfile(*args, **kwargs).reshape(-1, 4))


## (w, x, y, z) --> Versor
wxyz2versor = Versor.fromwxyz


## (x, y, z, w) --> Versor
xyzw2versor = Versor.fromxyzw


## The Real Quaternion Basis (depracting and moving to class).
W = Versor(scalar.ONE, vector.NULL)

## \f$ i \f$
I = Versor(scalar.ZERO, vector.Vector(1, 0, 0))

## \f$ j \f$
J = Versor(scalar.ZERO, vector.Vector(0, 1, 0))

## \f$ k \f$
K = Versor(scalar.ZERO, vector.Vector(0, 0, 1))

from collections import namedtuple
Basis = namedtuple('_Basis', 'w i j k')

## The basis of __H__: (w, i, j, k) with:\n
# \f$ i^2 = j^2 = k^2 = -w^2 = -1 \f$ \n
# \f$ ij = k,\  jk = i,\  ki = j \f$
BASIS = Basis(W, I, J, K)


## Parity Operator: \f$ R_+ \oplus R_-^3 \f$
# \param q A quaternion
# \returns P(q) (which is qbar)
def parity(q):
    """q' = P(q) is the parity operator."""
    from ..euclid import parity as P
    return type(q)(*map(P, q.iter()))
