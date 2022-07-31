"""The Space Curve is composed of an array_like Vector and its 1-d time
array:

>>>SpaceCurve(vector, t=t)

SpaceCurves are composed of their 2 parts, as inheritance from vectors
is problematic.

Moreover, SCs do everything, and perhaps that's not right. It could be that
transformations only by call are defined."""
#pylint: disable=C0301,R0903,C0111

## \namespace geo::metric::frenet_serret::motion SpaceCurve's and their
# Tangent Frames

import functools
import itertools
import operator
import types

from ..euclid import tensor
from ..euclid import vector
from ..euler import charts

from .dxdt import deriv as D

## UP in the Body Frame
UP = vector.Z

## Here, DOWN __is__ minus UP.
DOWN = -UP


## tangent-left-up to tangent-right-down transformation
tlu2trd = charts.roll(180)

## tangent-right-down to tangent-left-up  transformation
trd2tlu = ~tlu2trd


## <a href="http://mathworld.wolfram.com/AngularVelocity.html">
# Angular Velocity</a>
# \param r Position
# \param v Velocity
# \returns \f$ \vec{\omega} = \vec{v}/r \f$ \n
def angular_velocity(r, v):
    """Angular Velocity"""
    return v / abs(r)


##  <a href="http://mathworld.wolfram.com/Speed.html"> Speed</a>
# \param v Velocity
# \returns \f$ |\vec{v}| \f$ \n
def speed(v):
    """Speed"""
    return abs(v)


##  <a href="http://mathworld.wolfram.com/TangentVector.html">Tangent Vector
#  </a>
# \param v Velocity
# \returns \f$ \hat{T} \equiv \hat{v} \f$ \n
def tangent(v):
    """Tangent"""
    return v.hat()


## Why? so later I can add better derivatives as keywords
def _vxa(v, a):
    return v ^ a


##  <a href="http://mathworld.wolfram.com/BinormalVector.html">Binormal
#  vector</a> as a hat() decorated velocity_X_acceleration()
# \param v Velocity
# \param a Acceleration
# \returns  \f$ \hat{B} = \vec{B}/|\vec{B}| \f$ \n
def binormal(v, a):
    """binormal"""
    return tangent(_vxa(v, a))


## <a href="http://mathworld.wolfram.com/NormalVector.html">Normal Vector
# </a>
# \param v Velocity
# \param a Acceleration
# \returns  \f$\hat{N}\equiv\hat{\dot{\hat{T}}}=\frac{d\hat{T}}{dt}/
# |\frac{d\hat{T}}{dt}|\f$\n
def normal(v, a):
    """Normal"""
    return binormal(v, a) ^ tangent(v)


##  <a href="http://mathworld.wolfram.com/Curvature.html">The Curvature</a>
# \param v Velocity
# \param a Acceleration
# \returns \f$ \kappa = \frac{|\vec{v} \times \vec{a}|}{|\vec{v}|^3} \f$ \n
def curvature(v, a):
    """curvature"""
    return abs(_vxa(v, a)) / speed(v)**3


##  <a href="http://mathworld.wolfram.com/Torsion.html">Torsion</a>
# \param v Velocity
# \param a Acceleration
# \param j jerk
# \returns \f$ \tau = [\vec{v}, \vec{a}, \dot{a}]/\sigma^2  \f$ \n
# @image html torusan.gif
def torsion(v, a, j):
    """Torsion"""
    return (v ^ (a ^ j)) / _vxa(v, a)**2


## <a href=" http://mathworld.wolfram.com/Centrode.html">The Centrode</a>
def centrode():
    """This code was forgotten."""


##  <a href="http://mathworld.wolfram.com/Curvature.html">The Curvature</a>
# \param v Velocity
# \param a Acceleration
# \param j Jerk
# \returns  \f$ \vec{C} = \tau \vec{T} + \kappa \vec{B} \f$ \n
def darboux(v, a, j):
    """The darboux, from v, a, j."""
    return torsion(v, a, j)*tangent(v) + curvature(v, a)*binormal(v, a)


## TNB
# \param v Velocity
# \param a Acceleration
# \return \f$ {bf \hat{x}T} + {bf \hat{y}N} + {bf \hat{z}B} \f$
def TNB(v, a):
    """TNB, from v,  a"""
    from ..euclid.tensor import ziprows
    return ziprows(tangent(v), normal(v, a), binormal(v, a))


## T'N'B'
# \param v Velocity
# \param a Acceleration
# \param j Jerk
# \return \f$ {bf \hat{x}T'} + {bf \hat{y}N'} + {bf \hat{z}B'} \f$
def TNBprime(v, a, j):
    """T'N'B', from v, a, j"""
    from ..euclid.tensor import ziprows
    print "CHECK!"
    d = darboux(v, a, j)
    return ziprows(*map(d.__xor__, TNB(v, a).rows()))


##  This decorator normalizes the output with abs(self):
# \param func \f$ f(x) \f$
# \returns \f$ f \rightarrow \frac{f(x)}{||x||} \f$
def normalized(func):
    """ func() --> func()/||self|| decorator"""
    @functools.wraps(func)
    def nfunc(self, *args, **kwargs):
        return func(self, *args, **kwargs)/abs(self)
    return nfunc


##  This decorator normalizes the output with abs(output)
# \param func \f$ \vec{f}(x) \f$
# \returns \f$ f \rightarrow \frac{\vec{f}(x)}{||f(x)||} \f$
def _hat(func):
    """func() --> func().hat() decorator"""
    @functools.wraps(func)
    def hfunc(self, *args, **kwargs):
        return func(self, *args, **kwargs).hat()
    return hfunc


##  This decorator returns the reciprocal of the function
# \param func \f$ f(x) \f$
# \returns \f$ f \rightarrow \frac{1}{f(x)} \f$
def radius_of(func):
    """func --> 1/func decorator"""
    @functools.wraps(func)
    def rfunc(self, *args, **kwargs):
        return func(self, *args, **kwargs)**(-1)
    return rfunc


#  This decorator returns the magnitude of the vector function
# \param func \f$ \vec{f}(x) \f$
# \returns \f$ \vec{f} \rightarrow ||f(x)|| \f$
def magnitude(func):
    """func --> |func| decorator"""
    @functools.wraps(func)
    def mfunc(self, *args, **kwargs):
        return abs(func(self, *args, **kwargs))
    return mfunc


## Comparision Decoration (Deprecated)
def comp(method):
    def comp_method(self, other):
        try:
            other.rank
        except AttributeError:
            result = self.vector
        else:
            try:
                operand = other.vector
                if isinstance(operand, types.MethodType):  # external poly
                    raise AttributeError
            except AttributeError:
                operand = other
            result = method(self, other)(self.vector, operand)
        if result.rank == 1:
            result = type(self)(result, t=self.t)
        return result
    return comp_method


## Direct method transform from operator module.
def comp2_0(func):
    def comp_method(self, other):
        try:
            other.rank
        except AttributeError:
            result = self.vector
        else:
            try:
                operand = other.vector
                if isinstance(operand, types.MethodType):  # external poly
                    raise AttributeError
            except AttributeError:
                operand = other
            result = func(self.vector, operand)
        if result.rank == 1:
            result = type(self)(result, t=self.t)
        return result
    return comp_method


## Direct method transform from operator module.
def comp2(func):
    def comp_method(self, other):
        result = func(self.vector, other)
        if result.rank == 1:
            result = type(self)(result, t=self.t)
        return result
    return comp_method


## The Space Curve is an array_like Vector, and its time coordinate.
class SpaceCurve(object):
    """v = SpaceCurve(vector, t=t)"""

    ## Need this for composition, but not inheritance?
    rank = 'spacecurve'  # special hack workaround

    ## A Vector (element-wise) with optional t parameter
    def __init__(self, vector_, t=None):
        self.vector = vector_
        if t is None:  # kwd
            from numpy import arange
            t = arange(len(vector_))
        self.t = t

    def __getattr__(self, attr):
        try:
            result = getattr(self.vector, attr)
        except AttributeError as err:
            raise err
        try:
            result = result.spacecurve(self.t)
        except Exception:
            pass
        return result

    def __abs__(self):
        return abs(self.vector)

    def __str__(self):
        return "\n".join(itertools.imap("{}\t: {}".format, self.t, self.vector))

    __repr__ = __str__

    def __pos__(self):
        return self

    ## \f$ -{\bf \vec{v}}(t) \f$
    def __neg__(self):
        return type(self)(-(self.vector), t=self.t)

    __add__ = comp2(operator.add)
    __sub__ = comp2(operator.sub)
    __mul__ = comp2(operator.mul)
    __xor__ = comp2(operator.xor)
    __rshift__ = comp2(operator.rshift)
    __lshift__ = comp2(operator.lshift)

    def __rrshift__(self, other):
        return (other >> self.vector).spacecurve(self.t)

    def __rlshift__(self, other):
        return (other << self.vector).spacecurve(self.t)

    __and__ = comp2(operator.and_)
    __or__ = comp2(operator.or_)

    def __pow__(self, n):
        return self.vector**n

    def __rmul__(self, a):
        return self / (1./a)

    __div__ = comp2(operator.div)

    ## Not a Liskov safe architecture-- you cannot have a non-array
    # space curve-- since it's a sampled curve, it has to be a sequence.
    def __getitem__(self, index):
        return type(self)(self.vector[index], t=self.t[index])

    def __setitem__(self, index, value):
        self.vector[index] = value.vector
        self.t[index] = value.t

    ## len() --> length of the attributes (and they have to be equal)
    def __len__(self):
        if len(self.vector) == len(self.t):  # raise
            return len(self.t)
        raise ValueError("There are {} vectors and {} time stamps.".format(
            len(self.vector), len(self.t)))

    ## \f$ {\bf \hat{v}}(t) \f$
    def hat(self):
        """Unitize vector part."""
        return type(self)(self.vector.hat(), t=self.t)

    def broadcast(self, func):
        """Broadcast func to components."""
        return type(self)(self.vector.broadcast(func), self.t)

    ## \f$ \vec{v} \equiv \frac{d\vec{r}}{dt} \f$ \n
    #  <a href="http://mathworld.wolfram.com/Velovelcity.html for details.">
    #  Velocity</a>
    def velocity(self):
        """1st derivative of self()"""
        return self.broadcast(D)

    ## \f$ \vec{a} \equiv \frac{d\vec{v}}{dt} = \frac{d^2\vec{r}}{dt^2} \f$ \n
    #  <a href="http://mathworld.wolfram.com/Acceleration.html for details.">
    #  Velocity</a>
    def acceleration(self):
        """2nd derivative of self()"""
        return self.velocity().velocity()

    ## \f$ \vec{j} \equiv \frac{d\vec{a}}{dt} = \frac{d^3\vec{r}}{dt^3} \f$ \n
    #  <a href="see http://en.wikipedia.org/wiki/Jerk_(physics)">Jerk</a>, as
    #  a derivative() decorated SpaceCurve.__call__
    def jerk(self, v=None, a=None, j=None):
        """3rd derivative of self()"""
        return self.acceleration().velocity()

    def speed(self):
        """Speed"""
        return speed(self.velocity())

    def tangent(self):
        """Tangent is the 'hat'-ed Velocity"""
        return tangent(self.velocity())

    def normal(self):
        """Normal is the cross product of the binormal and the Tangent"""
        return normal(self.velocity(), self.acceleration())

    def angular_velocity(self):
        """angular_velocity os the normalized velocity"""
        return angular_velocity(self.vector, self.velocity())

    ## \f$ \vec{\alpha} = \vec{\omega}' = \vec{a}/r \f$ \n
    # <a href="http://mathworld.wolfram.com/AngularAcceleration.html">
    # Angular Acceleration</a>
    @normalized
    def angular_acceleration(self):
        """Angular Acceleration is the normalized Acceleration"""
        return self.acceleration()

    def binormal(self):
        """Binormal is the 'hat'-ed velocity_X_acceleration"""
        return binormal(self.velocity(), self.acceleration())

    def torsion(self):
        """torsion: -- it's along story"""
        return torsion(self.velocity(), self.acceleration(), self.jerk())

    def curvature(self):
        """||V X A||/||V||**3"""
        return curvature(self.velocity(), self.acceleration())

    def darbaux(self):
        """The Darbaux."""
        return darboux(self.velocity(),
                       self.acceleration(),
                       self.jerk())

    ## \f$ \sigma = 1/\tau \f$ \n
    #  <a href="http://mathworld.wolfram.com/RadiusofTorsion.html">Radius of
    #  torsion</a> as a radius_of() decorated torsion().
    @radius_of
    def radius_of_torsion(self):
        """radius_of() torsion"""
        return self.torsion()

    ## \f$ \rho^2 = 1/ |v \times a |^2 \f$ \n
    #  <a href="http://mathworld.wolfram.com/RadiusofCurvature.html">
    #  Radius of curvature</a> as a radius_of() decorated curvature()
    @radius_of
    def radius_of_curvature(self):
        """radius_of() curvature"""
        return self.curvature()

    ## \f$ s = \int_{\gamma} ds = \int_{\gamma} |\dot{\vec{r}}| \f$\n
    #  <a href="http://mathworld.wolfram.com/ArcLength.html">Arc Length</a>
    #  (computed for fixed time steps).
    def arc_length(self, axis=None):
        """TODO: use scipy to integrate to make a nested function...."""
        return self.speed().cumsum(axis=axis)

    ## \f$ {\bf T} = \hat{x}\hat{T}+\hat{y}\hat{N} + \hat{z}\hat{B} \f$ \n
    #  <a href="http://mathworld.wolfram.com/Trihedron.html">The Trihedron
    #  Tensor</a>
    def TNB(self):
        """TNB"""
        return TNB(self.velocity(), self.acceleration())

    ## \f$ \kappa \f$, read-only curvature()
    @property
    def kappa(self):
        """kappa"""
        return self.curvature()

    ## \f$ \tau \f$, read-only torsion()
    @property
    def tau(self):
        """tau"""
        return self.torsion()

    ## \f$ \vec{T} \f$, read-only tangent()
    @property
    def T(self):
        """T"""
        return self.tangent()

    ## \f$ \vec{B} \f$, read-only normal()
    @property
    def N(self):
        """N"""
        return self.normal()

    ## \f$ \vec{B} \f$ , read-only binormal()
    @property
    def B(self):
        """B"""
        return self.binormal()

    @property
    def Tprime(self):
        """T'"""
        return self.kappa*self.N

    @property
    def Nprime(self):
        """N'"""
        return -self.kappa*self.T + self.tau*self.B

    @property
    def Bprime(self):
        """B'"""
        return -self.tau*self.B

    ## \f$ {\bf \vec{T}'} = {\bf \vec{\omega} \times \vec{T}} \f$ \n
    #  \f$ {\bf \vec{N}'} = {\bf \vec{\omega} \times \vec{N}} \f$ \n
    #  \f$ {\bf \vec{B}'} = {\bf \vec{\omega} \times \vec{B}} \f$ \n
    #  <a href="http://en.wikipedia.org/wiki/Frenet-Serret_formulas">
    #  Frenet-Serret</a> formulae.
    def TNBprime(self):
        """THE T'N'B'"""
        return TNBprime(self.velocity(),
                        self.acceleration(),
                        self.jerk())

    ## Matplotlib plotter
    def plot3d(self, f=0):
        """Requires matplotlib."""
        import pylab
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d import axes3d
        fig = plt.figure(f)
        ax = Axes3D(fig)
        ax.plot(self.x, self.y, self.z)
        pylab.xlabel("x")
        pylab.ylabel("y")
        pylab.show()
        return ax

    ## return space curve
    def spacecurve(self):
        """self"""
        return self

    ##  Compute tangent plane triplet from SpaceCurve \f${\bf \vec{r}}\f$:\n
    #  \f$ \bf{ \hat{x} } \propto {\bf \hat{r'}} \f$ \n
    #  \f$ {\bf \hat{z}} \propto {\bf \hat{z}} -
    #  {\bf (\hat{z}\cdot\hat{x})\hat{x}}  \f$ \n
    #  \f$ \bf{ \hat{y} }= {\bf \hat{z} \times \hat{x} } \f$ \n
    #  (where "prime" is differentiation via prime())\n
    #  keyword z defines body z coordinate.
    def tangent_frame_triplet(self, z=DOWN):
        """i, j, k = r.tangent_frame_triplet([z=DOWN])

        i, j, k for a right handed orthonormal triplet, with:

        i   parallel to r.velocity()
        j   is k ^ i (cross product)
        k   is keyword z's normalized vector rejection of i
        """
        i = self.tangent().vector
        k = (z.rejection(i)).hat()
        j = (k ^ i)

        return i, j, k

    ## Tangent, Level, Up' to Level \n: Rotation connecting
    #  tangent-to-ellipsoid to tangent-to-motion frame \n computed from
    #  itertools.product:\n take all dot product combinations (as numbers,
    #  not vector.Scalar()) and make a euclid.Matrix().
    def tlu2level(self):
        """Tangent Left Up to Level."""
        return tensor.Tensor(
            *[(e_body*e_level).w for e_body, e_level in itertools.product(
                self.tangent_frame_triplet(z=UP),
                vector.BASIS)]).versor()

    ## invert tlu2level()
    def level2tlu(self):
        """Level to Tangent Left Up."""
        return ~(self.tlu2level())

    ## compose level2tlu() and tlu2rd()
    def level2trd(self):
        """Level to Tangent Right Down."""
        #pylint: disable=E1102
        return self.level2tlu.compose(tlu2trd())

    ## Compute level frame to body frame-- with keyword defined system
    #  (TLU or TRD)
    def level2body(self, imu, method='level2tlu'):
        """Level to Body frame."""
        return (operator.methodcaller(method)(self)).compose(imu)

    ## To Be Debugged
    def level2body_affine(self, imu, method='level2tlu'):
        """Level to body affine Transform."""
        from ..euclid.affine import Affine
        print "order not debugged"
        R = self.level2body(imu, method=method)
        T = -(~R)(self)
        return Affine(R, T)
