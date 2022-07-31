"""The Real Vector(x, y, z) class lives here, along with some convinience
functions:

dot, cross, outer
scalar_triple_product
vector_triple_product

Proj, Rej, Ref (projection, rejection, and reflections),

and the constants:

NULL = Vector(0, 0, 0)
X, Y, Z
BASIS = X, Y, Z
"""
#!pylint: disable=C0301, R0904, R0201

## \namespace geo::metric::euclid.vector The Cartesian Vector and/or
# Vector Field.
import fractions
import functools
import sys

from ...utils import exceptions
from ...utils.trig import cos, sin, TWO_PI, CIRCUMFERENCE
from ...utils.arraylike import enlister

from . import euclid

from .kant import Polar


__all__ = ('scalar_triple_product', 'vector_triple_product',
           'scalar_quadruple_product', 'vector_quadruple_product',
           'Proj', 'Rej', 'Ref',
           'Vector', 'NULL', 'BASIS', 'X', 'Y', 'Z')

power_table = {-2: lambda v: 1/v**2,
               -1: lambda v: ~v,
               1: lambda v: v,
               2: lambda v: v*v,
               0.5: lambda v: v.sqrt()}


## This decorator makes a method that gets self.__module__.<name>,
# and calls it; why? because we don't know If the vector lives in
# \f$\mathbb R^3\f$ or \f$\mathbb C^3\f$ until runtime.
# \param name The name of the module function
# \returns method calling module.<name>(self)
def load_local_function(name):
    """method = load_local_function(name)"""

    # E.G: Vector.projector calls vector.Proj, while
    # Gibbs.projector calls gibbs.Proj
    @functools.wraps(getattr(sys.modules[__name__], name))
    def method(self):
        """Get the module function with name name"""
        return getattr(sys.modules[self.__module__], name)(self)

    return method


## Explicit Inner Product
# \param u Vector \f$ {\bf \vec{u}} \f$
# \param v Vector \f$ {\bf \vec{v}} \f$
# \returns scalar.Scalar
# \f$ s = {\bf \vec{u}\cdot\vec{v}}
# \equiv u_xv_x +  u_yv_y + u_zv_z\f$
# @image html dot.png
# @image latex dot.png
def dot(u, v):
    """Scalar product. The hard way:
    return euclid.ZOO[0](reduce(operator.add,
                                itertools.starmap(
                                    operator.mul,
                                    itertools.izip(u.iter(), v.iter()))))
    """
    # no C3 heroics here.
    return euclid.ZOO[0](u.x*v.x + u.y*v.y + u.z*v.z)


## Explicit Outer product: \f$ \vec{u}\vec{v} \f$
# \param u \f$ u_i \f$
# \param v \f$ v_j \f$
# \returns tensor.Tensor \f$ T_{ij} = u_iv_j \f$
# \throws TypeError on non-vector input
# \throws AttributeError on un Dx'd exception.
def outer(u, v):
    """outer (dyadic) product, conjugates 2nd arg.

    return euclid.ZOO[2](*itertools.starmap(operator.mul,
                                            itertools.product(u.iter(),
                                                              v.iter())))
    """
    try:
        # Get 9 dyads- fail If input doesn't have components: x, y, z.
        inputs = (u.x*v.x, u.x*v.y, u.x*v.z,
                  u.y*v.x, u.y*v.y, u.y*v.z,
                  u.z*v.x, u.z*v.y, u.z*v.z)
    except AttributeError as err:
        pass
    else:
        try:
            # Get tensor class: fails If u or v where quaternions.
            cls = euclid.ZOO[u.rank + v.rank]
        except TypeError as err:
            pass
        else:
            # RETURN: make the instance of Tensor with 9 inputs.
            return cls(*inputs)
    # Dx the exeception and raise it with an explanation.
    for count, item in enumerate([u, v]):
        if euclid.rank(item) != 1:  # raise
            msg = "Expected a Vector, got {} in parameter #{}"
            raise TypeError(msg.format(type(item), count))
    # Raise an un-Dx'd error.
    raise err


## \f$ {\bf [\vec{a}, \vec{b}, \vec{c} ]} \equiv {\bf \vec{a} \cdot
#  (\vec{b} \times \vec{c}})  \f$
# \param a Vector \f$ {\bf \vec{a}} \f$
# \param b Vector \f$ {\bf \vec{b}} \f$
# \param c Vector \f$ {\bf \vec{c}} \f$
# \returns (pseudo) scalar.Scalar \f$ \pi = \epsilon_{ijk}a_ib_jc_k \f$
# @image html scalar_triple_product.jpg
# @image latex scalar_triple_product.jpg
def scalar_triple_product(a, b, c):
    """scalar =  scalar_triple_product(a, b, c)
    a, b, c      Vector inputs:
    a * (b ^ c)"""
    return dot(a, cross(b, c))


## \f${\bf  \vec{a} \times (\vec{b} \times \vec{c})} \f$
# \param a Vector \f$ {\bf \vec{a}} \f$
# \param b Vector \f$ {\bf \vec{b}} \f$
# \param c Vector \f$ {\bf \vec{c}} \f$
# \returns \f$ v_i=\epsilon_{ijk}\epsilon_{klm}a_jb_lc_m \f$ Vector
# @image html vector_triple_product.gif
# @image latex vector_triple_product.gif
def vector_triple_product(a, b, c):
    """vector =  vector_triple_product(a, b, c)
    a, b, c      Vector inputs -->
    a ^ (b  ^ c)"""
    return reduce(cross, (c, b, a))


## \f$  {\bf (\vec{a}\times \vec{b}) \cdot
# (\vec{c} \times \vec{d}) } =
# ({\bf\vec{a}\cdot\vec{c}})({\bf\vec{b}\cdot\vec{d}}) -
# ({\bf\vec{b}\cdot\vec{c}})({\bf\vec{a}\cdot\vec{d}}) \f$
# \param a Vector \f$ {\bf \vec{a}} \f$
# \param b Vector \f$ {\bf \vec{b}} \f$
# \param c Vector \f$ {\bf \vec{c}} \f$
# \param d Vector \f$ {\bf \vec{d}} \f$
# \returns  scalar.Scalar
# \f$ w = \epsilon_{ijk}\epsilon_{imn}a_jb_kc_md_n \f$
def scalar_quadruple_product(a, b, c, d):
    """scalar =  scalar_quadruple_product(a, v, w, r)
    a, b, c, d     Vector inputs -->
    (a ^ b) * (c ^ d)"""
    return dot(cross(a, b), cross(c, d))


## \f$  {\bf (\vec{a}\times \vec{b}) \times
# (\vec{c} \times \vec{d}) }  \f$
# \param a Vector \f$ {\bf \vec{a}} \f$
# \param b Vector \f$ {\bf \vec{b}} \f$
# \param c Vector \f$ {\bf \vec{c}} \f$
# \param d Vector \f$ {\bf \vec{d}} \f$
# \returns Vector
# \f$ v_i = \epsilon_{ijk}\epsilon_{jlm}\epsilon_{knp}a_lb_mc_nd_p \f$
def vector_quadruple_product(a, b, c, d):
    """vector =  vector_quadruple_product(a, b, c, d)
    a, b, c, d     Vector inputs -->
    (a ^ b) ^ (c ^ d)"""
    return cross(cross(a, b), cross(c, d))


## External poly work around, when a spacecurve do the calling
_spacecurve = lambda x, y: dot(x.vector, y)


## Cross product: \f$ {\bf \vec{v'} } =
# {\bf \vec{a} \times \vec{b}} \f$
# \param a \f$ a_i \f$
# \param b \f$ b_j \f$
# \returns Vector \f$ c_i = \epsilon_{ijk}a_jb_l \f$
# @image html cross.png
# @image latex cross.png
def cross(a, b):
    """explicit cross(a, b) product:

    c.i = epsilon.ijk * a.j * b.k"""
    return type(a)(a.y*b.z - a.z*b.y,
                   a.z*b.x - a.x*b.z,
                   a.x*b.y - a.y*b.x)


## Look up table for scalar product functions (per number of arguments)
_scalar_product_table = {2: dot,
                         3: scalar_triple_product,
                         4: scalar_quadruple_product}

## Look up table for scalar product functions (per number of arguments)
_vector_product_table = {2: cross,
                         3: vector_triple_product,
                         4: vector_quadruple_product}


## Make a scalar from N vectors
# \param args Vector objects
# \returns scalar.Scalar per ::_scalar_product_table
# \throws ValueError If it can't be done.
def scalar_product(*args):
    """scalar = scalar_product(a, b, ...)"""
    try:
        func = _scalar_product_table[len(args)]
    except IndexError:
        raise ValueError("Can't make scalar from {} vectors".format(
            len(args)))
    return func(*args)


## Make a vector from N vectors
# \param args Vector objects
# \returns Vector per ::_vector_product_table
# \throws ValueError If it can't be done.
def vector_product(*args):
    """vector = vector_product(a, b, ...)"""
    try:
        func = _vector_product_table[len(args)]
    except IndexError:
        raise ValueError("Can't make vector from {} vectors".format(
            len(args)))
    return func(*args)


## Explicit v * M (is faster then generalized Form)
# \param v Vector
# \param m tensor.Tensor
# \returns  Vector \f$ v_i = M_{ij} v_j \f$
# \throws TypeError when rank-2 arg is not.
# \throws AttributeError on undiagnosed exception.
def anterior_product(v, m):
    """v' = anterior_product(v, m)
          = v*m
    v        Vector
    m        Tensor (rank-2)"""
    try:
        # step-1 is ask "m" for its tensor components
        x = m.xx * v.x + m.yx * v.y + m.zx * v.z
        y = m.xy * v.x + m.yy * v.y + m.zy * v.z
        z = m.xz * v.x + m.yz * v.y + m.zz * v.z
    except AttributeError as err:
        if euclid.rank(v) != 1:  # raise
            err = TypeError("Expected a Vector, got {}".format(type(v)))
        elif euclid.rank(m) != 2:  # raise
            err = TypeError("Expected a rank-2, got {}".format(type(m)))
        else:
            pass
        raise err
    try:
        # step-2 make a vector out the results, which fails If v is a versor.
        result = type(v)(x, y, z)
    except TypeError:
        err = TypeError("Expected a Vector, got {}".format(type(v)))
    return result


## Left vs Right dilation is under development: it matters, For example:
# \f$ \sigma_2 {\bf \overrightarrow{\sigma}} \ne
# {\bf \overrightarrow{\sigma}} \sigma_2 \f$--
# that is, If your doing quantum mechanics, the vector's components may
# be matrix-like.
_vector_dilation = euclid.Tensor_.right_dilation


## Explicit dilation by a scalar.Scalar
# \param v Vector
# \param s scalar.Scalar
# \returns Vector \f$ v_i \rightarrow sv_i \f$
def _vector_times_scalar(v, s):
    """v' = _vector_times_scalar(v, s)

    v, v'  are Vector
    s       is a Scalar"""
    return _vector_dilation(v, s.w)


## \f$ {\bf {\rm Proj}_a} = {\bf\hat{a}\hat{a}} \f$ \n
# Projection Operator onto a
# \param  a A Vector
# \returns The geo.metric.euclid.tensor.Tensor that projects vectors
# onto \f$ {\bf \vec{a}} \f$
# @image html proj.jpeg
# @image latex proj.jpeg
def Proj(a):
    """P = Proj(a) is the projection operator onto a

    >>>v = P(b) == b >> a

    is the projection of b onto a."""
    return a.hat().dyad()


## \f$ {\bf {\rm Rej}_a} = -({\bf \hat{a} \otimes})^2 \f$
# Rejection onto a
# \param  a A Vector
# \returns The geo.metric.euclid.tensor.Tensor that rejects vectors
# onto \f$ {\bf \vec{a}} \f$
# @image html rej.jpeg
# @image latex rej.jpeg
def Rej(a):
    """R = Rej(a) is the rejection operator onto a

    >>>v = R(b) == b << a

    is the rejection of b onto a.
    """
    return -(a.hat().cross()**2)


## \f$ {\bf {\rm Ref}_a} = {\bf {\rm Proj}_a} - {\bf {\rm Rej}_a} \f$ \n
# Reflection about a's orthogonal plane
# \param  a A Vector
# \returns The geo.metric.euclid.tensor.Tensor that refects vectors
# about the plane orthogonal \f$ {\bf \vec{a}} \f$
def Ref(a):
    """R = Ref(a) is the reflection operator onto a

    >>>v = R(b) == b << a

    is the reflection of b onto a, which means it's been relfected in
    by the plane perp-to-a."""
    return Proj(a) - Rej(a)


## Dual Tensor
# \f$ {\bf \tilde{v}} \equiv {\bf \vec{v} \otimes} \f$
# \param vec \f$ v_i \f$
# \returns tensor.Tensor \f$ \epsilon_{ijk}v_k \f$
def dual(vec):
    u"""tensor = dual(vec).

    T_ij = \u2211_k[\u03b5_ijk * v_k]

    Inverts tensor.dual(t)."""
    #
    from .three import EPSILON
    return EPSILON * vec


## quotient--is not bound to Vector.__div__
# \param a Vector
# \param b Vector
# \returns geo.desic.hamilton.quaternion.Quaternion \n
# \f$ {\bf q} =
# \frac{1}{||\vec{b}||^2}(\vec{a}\cdot\vec{b}, \vec{a}\times\vec{b}) \f$
def quotient(a, b):
    """Returns a quaternion: a*b, a^b-- prolly junk"""
    from ...desic.hamilton import Quaternion
    scale = abs(a).w/abs(b).w
    ahat = a.hat()
    bhat = b.hat()
    return scale*Quaternion(ahat * bhat, ahat ^ bhat)


## Experimental: promote v to a q and pull the trigger
def _vpromo(v, q):
    """Implement: vector * quaternion by promoting vector to R(1, 3)."""
    # convert vector to quaternion and do quaternion multiplcation, c.f:
    # K = X * Y.versor(pi)
    return v.quaternion() * q


## <a href="http://en.wikipedia.org/wiki/Vector_(physics)">Vector</a> class
# transForms as \f$ v_i' = M_{ij}v_j \f$
# @image html vector.png
# @image latex vector.png
class Vector(euclid.Tensor_):
    u"""\u20D1v = Vector(x, y, z) is a vector with 3 attributes:

    v.x, v.y, v.z

    Vector operations are overloaded:

    "*"  does dilation by a scalar, dot product, matrix multiply For
         For rank 0, 1, 2  objects.

    For vector arguments, u:

    v*u --> v.dot(u)
    v ^ u --> v.cross(u)
    v & u --> v.outer(u)

    The methods cover all the manifestly covariant equation you can write down.

    abs(v)     --> a scalar
    abs(v).w  --> a regular number or array

    v.hat()   is the unit vector along v.
    v.versor(angle)  makes a versor that rotates around v by angle.
    v.polar()        makes a Polar object out of it (spherical coordinates)
    v.cylinder()     IBID, cylindrical coordinates
    ~v --> v.dual()  makes a matrix that does a cross product with v.
    v.right_quaternion() makes a right quaternion: q = (0, v)
    v.right_versor()    makes a right versor: q = (0, v.hat())

    v.spherical()    Spherical (fundamental SO(3)) representation

    v.angle(u)       Angle between v and u

    v.projection(u)  v >> u       Project v onto u
    v.rejection(u)   v << u       Reject v onto u
    v.reflection(u)  v | u        Reflect v wrt to u.


    v.spacecurve(t)  Paramatize SpaceCurve wrt to array 't'
    v.var()          makes the covariance matrix of an array_like vector.

    v.e(n)          Get n-th component 0-->v.x, 1-->v.y, 2-->v.z

    type(v).baseunit(n) n= 0,1,2 --> X, Y, Z

    and so on and so-Forth. It's ALL there.
    """
    from ...utils.exceptions import CliffordAlgebraError as DivisionError

    ## The ranked meta class figures out the indices
    __metaclass__ = euclid.ranked

    ## Vectors are rank 1
    rank = 1

    ## Sneaky an scandolous way to get quaternions to multiply vectors, see
    # The __R3__ Vector visits __H__, and this is one way to deal (isinstance
    # is another--but that is too easy).
    scalar = 0

    ## (continues) quaternion multiplication seemlessly-- this occurs b/c
    # versors and quaternions are different (the former is confined to the
    # unit hypersphere and can't be added and other things).
    measure = 1

    ## Reconstructor.
    # \brief Sometimes, what looks like a vector may really be a polar or
    # space curve, hence. this.
    # \relates Proj
    # \param vector_
    # \returns cls instance with vector_'s data.
    # \bug may need to call argument's vector property.
    @classmethod
    def fromvector(cls, vector_):
        """Construct from a vector (trivial)"""
        return cls(*vector_.iter())

    ## Construct from a 2 x 2 complex matrix
    # \returns \f$ \frac{-1}{2}(m_{01}+m_{10}){\bf \hat{x}} +
    #  \frac{-i}{2}(m_{01}-m_{10}){\bf \hat{y}} +
    #  \frac{-1}{2}(m_{00}-m_{11}){\bf \hat{z}} \f$
    @staticmethod
    def fromspinor(m):
        """Construct from a 2x2 matrix.

        Differs from quaternion contructor's vector part by 1j (c.f.,
        Clifford Algebra's' bi-vector vs. cross product)
        """
        from ..euler.hamilton import Versor
        return Versor.fromspinor(1j*m).vector

    ## \f$ \vec{v} = v_x \hat{e}_x + v_y \hat{e}_y + v_z \hat{e}_z \f$ \n
    # Explicit dependency injection as a linear combination of the
    # <a href="http://en.wikipedia.org/wiki/Standard_basis">standard basis</a>
    # \param x The x-component
    # \param y The y-component
    # \param z The z-component
    def __init__(self, x, y, z):
        ## \f$ v_x \equiv \vec{v} \cdot \hat{e}_0 \f$, x component
        self.x = x
        ## \f$ v_y \equiv \vec{v} \cdot \hat{e}_1 \f$, y component
        self.y = y
        ## \f$ v_z \equiv \vec{v} \cdot \hat{e}_2\f$, z component
        self.z = z

    ## R3 --> R
    # \param other Some Tensor_.
    # \return \f$ {\bf {\vec{V}(u)}} \equiv {\bf \vec{V}\cdot u}\f$
    def __call__(self, other):
        """V(u) = V * u,

        Result? 'u' decide."""
        return self * other

    ## Addition with non-cartesian front end, that is, other might
    # not be in cartessian coordinates.
    # \param (self) \f$ a_i \f$
    # \param other  \f$ b_i \f$
    # \returns \f$ c_i = a_i + b_i \f$
    # \throws geo.utils.exceptions.NonCovariantOperation
    # @image html vector_add.gif
    # @image latex vector_add.gif
    def __add__(self, other):
        """Vector addition, and only vector addition."""
        try:
            return super(Vector, self).__add__(other.vector)
        except AttributeError:
            from ...utils.exceptions import NonCovariantOperation
            raise NonCovariantOperation("Can't add %s to a Vector" %
                                        type(other).__name__)

    ## Subtraction passed off to addition.
    # \param (self) \f$ {\bf \vec{a}} \f$
    # \param other  \f$ {\bf \vec{b}} \f$
    # \returns \f${\bf\vec{a}}-{\bf\vec{b}}={\bf\vec{a}}+(-{\bf\vec{b}})\f$
    # \returns Vector difference.
    # @image html vector_sub.gif
    # @image latex vector_sub.gif
    def __sub__(self, other):
        return self + (-other)
        
    ## Define what you want str to look like.
    _format_str = "[{}]"

    ## str of components() (cleaned up....)
    def __str__(self):
        """See Vector._format_str For visual wrapper."""
        return self._format_str.format(", ".join(map(str, self.tolist())))

    ## Go from a covariant Vector to a contravariant vector
    # (for orthonormal coordinates --see geo.detic.christoffel)
    # \returns \f$ \frac{{\bf{\vec{v}}}}{||v||^2} \f$
    def __invert__(self):
        """Invert magnitude- for co/contra-variant in orthonormal"""
        return self / (self*self)

    ## This is a limited "pow"-- don't do silly exponents.
    # \param n An integer in (-2, 1, 2)
    # \returns It really depends
    # \throws UndefinedGeometricOperation
    # \throws ValueError
    def __pow__(self, n):
        try:
            return power_table[n](self)
        except KeyError:
            raise ValueError("Exponent violates: ||{}|| <= 2".format(n))


    ## Division of mostly done by super, but quotient() may
    # be invoked
    def __div__(self, other):
        """V/c is normal, V/V'?, see quotient."""
        if euclid.rank(other) == 1:  # external polymorphism.
            return quotient(self, other)
        return super(Vector, self).__div__(other)

    ## \f$ x/\vec{v} \f$ is a Vector.DivisionError
    # or call quaotient with external polymorphism :-(
    # \param other It doesn't matter
    # \throws Vector.DivisionError
    def __rdiv__(self, other):
        raise self.DivisionError("Use clifford.py for division by vectors.")

    ## Dot product, converts self to bra() and uses vector.dot().
    # \param other Another Vector \f$ |v\rangle \f$ (or self as default)
    # \returns \f$ \langle u| v \rangle \f$ 
    def dot(self, other=None):
        """scalar = u.dot(v) for vector v.

        u.dot() --> u.dot(u)."""
        return dot(self.bra(), self if other is None else other)  # kwd

    ## The hash table assigns multiplication based on the 2nd operands rank,
    # the space curve entry-- it's a problem with the architecture,
    # 1j is for quaternion sandwiches- it *was* a nice multimethod emulator,
    # but grew out of control- it is what it is.
    _dispatch_mul = {None: _vector_dilation,
                     0: _vector_times_scalar,
                     1: dot,  # Vector.dot, not vector.dot -a lil' messy.
                     2: anterior_product,
                     3: lambda self, other: other.jki * self,  # avoid inner
                     4: lambda self, other: other.jkli * self,  # IBID
                     'spacecurve': _spacecurve,
                     1j: _vpromo}

    ## \f$ c_{i} = \epsilon_{ijk}a_jb_k \f$ \n The (pseudo)Vector wedge,
    # cross product
    # \param other Another Vector (on None)
    # \returns Vector cross produce (or tensor.Tensor operator)
    def cross(self, other=None):
        """Cross product.
        All forms are equivalent:

        >>>w = u.cross(v)
        >>>w = u ^ v
        >>>w = cross(u, v)
        >>>w = u.wedge(v)  # full contraction with Levi-Civita.
        where w is a Vector--- there is no AxialVector class, and
        parity not considered."""
        return -self.dual() if other is None else cross(self, other)  # kwd

    ## u^v --> \f$ {\bf \vec{u} \times \ \vec{v}} \f$
    # \param other A Vector (Tensor)
    # \returns self^other cross() (super)
    # @image html wedge.png
    # @image latex wedge.png
    def __xor__(self, other):
        u"""(a ^ b)_i = \u2211_i\u2211_j(a_i * b_j * \u025b_ijk)"""
        try:  # external polymorphism...
            return cross(self, other)
        except AttributeError:
            return super(Vector, self).__xor__(other)

    ## Call dyad() for vector args, or pass is up to super.
    def __and__(self, other):
        """outer() product:
        (a & b)_ij = a_i * b_j"""
        try:  # external polymorphism (try vector, except tensor)
            return self.dyad(other)
        except AttributeError:
            return super(Vector, self).__and__(other)

    ## dyad coverts other to a bra() and then calls outer(): |self><other|
    # \param other=None  A Vector, or self.
    # \returns \f$ |u\rangle\langle v| \f$
    def dyad(self, other=None):
        """T = v.dyad(other=None) --> v.outer(other)

        If other is None, then use 'v':

        v.dyad() --> v & v --> v.outer(v).
        """
        return outer(self, (self if other is None else other).bra())  # kwd

    ## Define a rotation about \f$ \hat{v} \f$ \n, relative to kwarg:
    # circumference = \f$2\pi\f$
    # \param angle The angle to rotate round self
    # \param circumference=TWO_PI The units for the angle
    # \returns A versor.versor.Versor reppin' the rotation around self.
    # \throws ValueError for unknown circumference string.
    def versor(self, angle=TWO_PI/2, circumference=TWO_PI):
        """vector(angle, circumference=2*pi)

        return a unit quaternion (versor) that represents an
        alias rotation by angle about vector.hat().

        Circumference can also be a string from
        geo.utils.trig.CIRCUMFERENCE
        """
        from ..euler.hamilton import Versor
        try:
            f = TWO_PI/circumference
        except TypeError:
            try:
                f = TWO_PI/CIRCUMFERENCE[
                    str(circumference).lower().rstrip('s')]
            except (AttributeError, KeyError, TypeError):
                raise ValueError(
                    "Valid circumference names are: %s" %
                    str(CIRCUMFERENCE.keys()))

        # ZOO[0] is Scalar, _ones_like handles array cases.
        return Versor(euclid.ZOO[0](self._ones_like(cos(f*angle/2.))),
                      self.hat()*(sin(f*angle/2.)))

    ## Convert to a
    # <a href=
    # "http://en.wikipedia.org/wiki/Classical_Hamiltonian_quaternions#Right_versor">
    # right versor</a>
    # after normalization.
    # \returns quaternion() \f$ {\bf \vec{v}} \rightarrow (0, {\bf \hat{v}}) \f$
    def right_versor(self):
        """v.hat().quaternion()"""
        return self.hat().quaternion()

    ## Convert to a
    # <a href=
    # "http://en.wikipedia.org/wiki/Classical_Hamiltonian_quaternions#Right_quaternion">
    # right versor</a> (for transformation)\n
    # That is: as add a::ZERO Scalar part and don't normalize to unit
    # hyper-sphere-- but it is still a Versor, not a Quaternion
    # \returns versor.Versor
    # \f$ {\bf \vec{v}} \rightarrow (0, {\bf \vec{v}}) \f$
    def quaternion(self):
        """v.quaternion() ---> Versor(Scalar(0), v)"""
        from ..euler.hamilton import Versor
        return Versor(euclid.ZOO[0](self._ones_like(0.)), self)

    ## Need to track this down and deprecate it: grasmann needs to
    # call a vector's versor as a non-normalized versor. meanwhile.
    # quaternions need to use quaternions. it's a mess
    # \returns geo.desic.hamilton.Quaternion().
    def true_quaternion(self):
        """Deprecated?"""
        return self.quaternion().quaternion()

    ## 2 x 2 Complex Rep
    # \kwd j \f$ \frac{1}{2} \f$
    # \returns matrix \f$ \left[\begin{array}{cc}
    # z  & x - iy \\
    # x + iy & -z \end{array}\right ] (for 1/2).
    # \bug Need to align pauli and spinor method names
    def pauli(self, j=0.5):
        """pauli(j=1/2) is the spinoral rep of dimensions 2j+1."""
        # note the -1j factor-- that's from clifford algebra.
        return self.true_quaternion().spinor(j)/1j

    ## 2-component spinoral rep
    # \kwd alpha =0 is an optional phase factor.
    # \returns geo.desic.cartan.Spinor
    # @image html spinor.tiff
    # @image latex spinor.tiff
    def spinor(self, alpha=0):
        """2-component spinor = vector.spinor(alpha=0)

        That is, convert a vector into the equivalent spinor, with its
        little hidden variable (alpha--the phase of that extra DoF).

        The scalar component may differ from Vector.quaternion()-- that
        is TBD. (Spinors can be bijected to quaternions, so do they agree?)

        It think this takes   X --> W + I
        while quaterion() is  X --> I 
        """
        from geo.desic.hilbert.cartan import Spinor
        return Spinor.fromvector(self, alpha)

    ## dual() Tensor
    # \returns Tensor
    def dual(self):
        u"""T_ij = \u2211_k(v_k * \u025b_ijk)"""
        return dual(self.vector)

    ## A >> B == A.Proj(B)
    @euclid.wrank
    def __rshift__(self, other):
        if other.rank == 1:  # external polymorphism.
            return self.projection(other)
        return NotImplemented

    ## A << B == A.Rej(B)
    @euclid.wrank
    def __lshift__(self, other):
        if other.rank == self.rank:  # external polymorphism.
            return self.rejection(other)
        return NotImplemented

    ## A | B == A.Ref(B)
    @euclid.wrank
    def __or__(self, other):
        if other.rank == self.rank:  # external polymorphism.
            return self.reflection(other)
        return NotImplemented

    ## Projection operator from local module--we don't use super's b/c
    # we need to worry about Gibbs complex vectors calling this
    projector = load_local_function('Proj')

    ## Rejection operator from local module
    reflector = load_local_function('Ref')

    ## Rejection operator from local module
    rejector = load_local_function('Rej')

    ## Override euclid.Tensor_.rejection, for the Vector only op.
    # \param self \f$ {\bf \vec{a}}\f$
    # \param other \f$ {\bf \vec{b}}\f$
    # \return \f$ {\bf -(\vec{a}\times\hat{b})\times\hat{b}}\f$
    def rejection(self, other):
        """rejection self onto other."""
        v = other.hat()
        return -(self ^ v) ^ v

    ## Householder matrix
    # \returns \f$ I - 2 \frac{\vec{v}\vec{v}}{\vec{v}\cdot\vec{v}} \f$
    def householder(self):
        """todo: why is this -reflector()"""
        from .tensor import DELTA
        return DELTA - 2*self.dyad()/self.dot(self).w

    ## Gramian
    # \returns \f$ G(u, v) =| \vec{u} \times \ldots \ \vec{v} |^2 \f$
    def gramian(self, other):
        """The gramian of 2 vectors"""
        return abs(self ^ other)**2

    ## Transpose is self: there is no column or row-- that is nonsense,
    # there is only permutations of indices, and Sym(1) is a pretty nada
    # group.
    def transpose(self, index=0):
        """transpose is self"""
        if index:  # raise
            raise exceptions.IndicesError("Vector transpose is trivial.")
        return self

    ## Variance of a composite Vector.
    # \returns \f$\bar{(\vec{v}-\langle\bar{v}\rangle)
    # (\vec{v}-\langle\bar{v}\rangle)}\f$
    # For an iterable Vector, compute the variance tensor
    def var(self):
        """For an iterable vector, return a covariance matrix"""
        v = (self - self.mean())
        return (v & v).mean()

    ## Get likewise zeros to fill in matrices and quaternions\n-
    # this may not be the best way
    def _ones_like(self, constant=1):
        """private way to make zeros"""
        return constant + self.x*0

    ## Centralized random vector.
    # \returns \f$ \vec{\mu} = \vec{v} - \bar{v} \f$
    def centralized(self):
        """vector = vector - <vector>"""
        return self - self.mean()

    ## Self
    # \returns self
    @property
    def vector(self):
        """I am vector"""
        return self

    ## What to do after a sandwiched Grassmann multiplication, this
    # may be an explicit implementation of something that could by
    # implicitly polymorphic- this is what happens when you have
    # Clifford algebra w/o actually implementing it.
    @staticmethod
    def _qt(versor):
        """_qt: what is it? Here it takes a versor
        and gets the vector part."""
        return versor.vector

    ## Upgrayedd to a Space Curve.
    # \param t Optional time coordinate.
    # \returns frenet_serret.motion.SpaceCurve
    def spacecurve(self, t=None):
        """For array_like vectors, make a full fledged motion.SpaceCurve:

        spacecurve = vector.spacecurve([t=None])
        """
        from ..frenet_serret.motion import SpaceCurve
        return SpaceCurve(self, t=t)

    ## Convert to a numpy matrix?
    # \returns matrix/list
    @enlister
    def asmatrix(self):
        """as a numpy matrix."""
        from numpy import matrix
        return matrix(self.tolist()).T

    ## Spherical (Fundamental, complex) Representation:
    # \returns \f$T\f$ as a geo.metric.wigner.eckart.Vector
    def spherical(self):
        """As a spherical rep'd vector."""
        from ..wigner import eckart
        return eckart.Vector.fromcart(self)

    ## Project Irreducible Representations
    # \param m Azimuth Projection
    # \returns \f$ v^{(q=m)} \f$
    def irrep(self, m):
        """irrep(j=1, q=0) --> weight j, seniority index q."""
        return super(Vector, self).irrep(1, 0, m)

    ## Root Spinor (Quaternion rep-- differs from Cartan rep).
    # \returns \f$ q \f$ such that \f$ qq = \{0, {\bf \vec{v}}\} \f$
    def sqrt(self):
        """Quaternion (spinor) rep sqrt."""
        return self.quaternion() ** 0.5

    ## The is the conjugate()-- this exist to match biquaternions, which
    # have 2 conjugates.
    # \returns \f$ {\bf \vec{v}}^* \f$
    def complex_conjugate(self):
        """The complex_conjugate."""
        return self.conjugate()

    def test(self, *args):
        """A test method"""
        return lambda x: self, args

    ## Ket prepares an argument for dot product according to the module's
    # inner product:
    # \f$ \langle u, v \rangle \equiv {\bf\vec{u}\cdot\vec{v}} \f$
    # \param v \f$ v \f$
    # \returns \f$ v \f$
    def bra(self):
        """ket(v) --> v"""
        # This Vector is in R3, so there is NO conjugation.
        return self

    ## Complexify vectors inner product space:
    # \returns gibbs.Gibbs \f$ R^3 \rightarrow C^3 \f$
    def gibbs(self):
        """Same vector, different space."""
        from .gibbs import Gibbs
        return Gibbs(*self.iter())

    ## Make a scalar from N vectors
    scalar_product = scalar_product

    ## Make a vector from N vector
    vector_product = vector_product

    ## Convert representation to Spherical coordinates
    # \returns kant.polar.Polar
    def polar(self):
        """vector --> polar"""
        from .kant.spherical import Polar
        return Polar.fromvector(self)

    ## Convert representation to Cylindrical coordinates
    # \returns kant.cylindrical.Cylinder
    def cylinder(self):
        """vector --> cylinder"""
        from .kant.cylindrical import Cylinder
        return Cylinder.fromvector(self)
    
    ## Convert representation to Parabolic coordinates
    # \returns kant.parabolic.Parabola
    def parabola(self):
        """vector --> parabola"""
        from .kant.parabolic import Parabola
        return Parabola.fromvector(self)
    

    
## The NULL vector
NULL = Vector.null()


from collections import namedtuple
Basis = namedtuple('Basis', 'x y z')


## \f$ \hat{e}_x, \hat{e}_y, \hat{e}_z \f$ \n The
# <a href="http://en.wikipedia.org/wiki/Standard_basis">standard basis</a>,
BASIS = Basis(*Vector.basis())


## \f$ {\bf \hat{x}} \f$
X = BASIS[0]

## \f$ {\bf \hat{y}} \f$
Y = BASIS[1]

## \f$ {\bf \hat{z}} \f$
Z = BASIS[2]
