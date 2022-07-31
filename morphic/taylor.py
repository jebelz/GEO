"""The quadric (with no z terms) and it's properties.


This packages needs vectors with an x, y, z:

>>>v = Vector(x, y, z)

With methods:

>>>v.hat()

and operator overloads:

>>>-v  # negation
>>>v ^ u  # cross product
>>>v << u  # vector rejection = I - projection


For curvatures It needs tensors

>>>t = Tensor(xx, xy, ..., zz)

with a constructor:

>>>t.ziprows(v1, v2, v3)

Transposition:

>>t' = t.T   


"""
## \namespace geo.morphic.taylor Quadratic Taylor Expansions of DEMS
import abc
import collections

import numpy as np
from numpy import linalg

from . import gauss
from .monge import FTYPE, CTYPE

## In Progress: how to propagate NULLs..
DEFAULT = -100  # np.nan

## hypot(1, x)**2 is faster than 1+x**2, and faster if 1 is array_like x.
HYP2 = lambda dummy: np.hypot(1., dummy)**2


## A collection of derivatives.
Derivatives = collections.namedtuple("Derivatives",
                                     ("xx", "xy", "yy", "x", "y"))


## Wrapping decorator-- pretty basic
def wrap_gauss(func):
    def run_gauss(self):
        return func(*self.geomorphics())
    return run_gauss


## Methods for both singleton and array-like Quadric object
class Qbase(object):
    """Base class--- to be ABC's"""

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def shape_operator(self):
        """Abstract method is defined for Quadric and QuadricArray
        differently"""

    @staticmethod
    @abc.abstractmethod
    def _make_symmetric_matrix(E, F, G):
        """Abstract method is defined for Quadric and QuadricArray
        differently"""

    ## See Quadric.z and
    @abc.abstractproperty
    def z(self):
        """Abstract property is defined for Quadric and QuadricArray
        differently"""

    ## str(p) --> "p[0]x**2+p[1]xy+p[2]y**2+p[3]x+p[4]y+p[5]"
    def __str__(self):
        return (
            str(self[0]) + "\t x**2 +\n" + str(self[1]) + "\t xy +\n" +
            str(self[2]) + "\t y**2 +\n" + str(self[3]) + "\t x +\n" +
            str(self[4]) + "\t y +\n" + str(self[5])
            )

    ## @returns \f$\frac{\partial f}{\partial x}{\bf \hat{x}} +
    #  i\frac{\partial f}{\partial y}{\bf hat{y}}\f$
    def __complex__(self):
        return complex(self.x, self.y)

    ## \param x x-coordinate
    #  \param y y-coordinate
    #  \returns z z(x, y)
    def __call__(self, x, y):
        return (
            self[0] * x**2 +
            self[1] * x*y +
            self[2] * y**2 +
            self[3] * x +
            self[4] * y
            )

    ## \returns  \f$ (\frac{\partial f}{\partial x})_{\{(x,y)=(0,0)\}} =
    # p_3 \f$
    @property
    def x(self):
        """dQ/dx at (0,0)"""
        return self[3]

    ## \returns  \f$ (\frac{\partial f}{\partial y})_{\{(x,y)=(0,0)\}} =
    #  p_4 \f$
    @property
    def y(self):
        """dQ/dy at (0,0)"""
        return self[4]

    ## \returns  \f$ (\frac{\partial^2 f}{\partial x^2})_{\{(x,y)=(0,0)\}} =
    #  2p_0 \f$
    @property
    def xx(self):
        """d2Q/dx2 at (0,0)"""
        return 2*self[0]

    ## \returns  \f$ (\frac{\partial^2 f}{\partial y^2})_{\{(x,y)=(0,0)\}} =
    #  2p_2 \f$
    @property
    def yy(self):
        """d2Q/dy2 at (0,0)"""
        return 2*self[2]

    ## \returns \f$(\frac{\partial^2 f}{\partial x\partial y})_{\{(x,y)=
    #  (0,0)\}}=p_1\f$
    @property
    def xy(self):
        """d2Q/dxdy at (0,0)"""
        return self[1]

    ## \returns  \f$ \frac{\partial^2 f}{\partial y \partial x} =
    #  \frac{\partial^2 f}{\partial x \partial y} \f$
    yx = xy

    ## @returns  \f$ E = 1+f_x^2 \f$
    @property
    def E(self):
        """E = 1 + (dQ/dx)**2 at (0, 0)"""
        return HYP2(self.x)

    ## @returns  \f$ FJJ = f_x f_y \f$
    @property
    def F(self):
        """F = 1 + (d2Q/dxdy) at (0, 0)"""
        return self.x*self.y

    ## @returns  \f$ G = 1+f_y^2 \f$
    @property
    def G(self):
        """G = 1 + (dQ/dy)**2 at (0, 0)"""
        return HYP2(self.y)

    ## @returns  \f$ L' = f_{xx} = 2p_0 \f$
    @property
    def L1(self):
        """L' = d2Q/dx2 at (0,0)"""
        return self.xx

    ## @returns  \f$ M' = f_{xy} = p_1 \f$
    @property
    def M1(self):
        """M' = d2Q/dxdy at (0,0)"""
        return self.xy

    ## @returns  \f$ N' = f_{yy} = 2p_2 \f$
    @property
    def N1(self):
        """L' = d2Q/dy2 at (0,0)"""
        return self.yy

    ## @returns  \f$ N \equiv L'/||\nabla F|| \f$ see L1()
    @property
    def L(self):
        """L = L'/||gradient||"""
        return self.L1/self.slope3D()

    ## @returns  \f$ M \equiv M'/||\nabla F|| \f$ see M1()
    @property
    def M(self):
        """M = M'/||gradient||"""
        return self.M1/self.slope3D()

    ## @returns  \f$ N \equiv N'/||\nabla F|| \f$  see N1()
    @property
    def N(self):
        """N = N'/||gradient||"""
        return self.N1/self.slope3D()

    ## The First Fundamental Form
    ## @returns
    ## \f${\bf I}=
    ## \left[\begin{array}{cc}  E & F \\ F & G \end{array}\right]\f$\n
    ## NB: first_fundamental_form() computation is class-dependent
    @property
    def I(self):
        """I --> first fundamental form"""
        return self.first_fundamental_form()

    ## The Second Fundamental Form
    ## \returns  \f$ {\bf II} =
    ## \left[\begin{array}{cc} L&M\\ M&N \end{array}\right]\f$\n
    ## NB: second_fundamental_form() computation is class-dependent
    @property
    def II(self):
        """II --> second fundamental form"""
        return self.second_fundamental_form()

    ## The shape operator
    ## \returns \f$ {\bf S} = {\bf I}^{-1}{\bf II} \f$ \n
    ## NB: shape_operator() computation is class-dependent
    @property
    def S(self):
        """S = I**(-1)*II"""
        return self.shape_operator()

    ## TBD: comment
    def slope(self):
        """complex slope-- not used in production since it's computed
        before the Quadric is"""
        return CTYPE(self[3] + 1j * self[4])

    ## TBD: comment
    def diagonal(self):
        """ The diagonal term as xx + i*yy"""
        return CTYPE(self[0] + 1j * self[2])

    ## TBD: comment
    def cross(self):
        """ The xy cross term"""
        return FTYPE(self[1])

    def _grad(self):
        return self.x, self.y, self.z

    ## 3D gradient Vector.
    ## \returns \f$ \frac{\partial f}{\partial x}{\bf \hat{x}} +
    ## \frac{\partial f}{\partial y}{\bf \hat{y}} +
    ## \frac{\partial f}{\partial z}{\bf \hat{z}} \f$
    ## (three dimensional)
    def grad(self):
        """vector gradient, as a tuple"""
        from .. import Vector
        return Vector(*self._grad())

    ## \f$ \hat{N} = -\vec{\nabla}f(x, y, z) \f$
    def normal(self):
        """Normal."""
        return -self.grad().hat()

    ## \f$ \vec{s} = -\frac{\partial z}{\partial x}{\bf \hat{x}} -
    # \frac{\partial z}{\partial y}{\bf \hat{y}} \f$
    def steepest(self):
        """Steepest Descent."""
        from .. import Vector
        slope = self.slope()
        t = Vector(-slope.real, -slope.imag, 0.)
        try:
            return t.hat()
        except ZeroDivisionError:
            return Vector(0, 0, 1)  #UserWarning("UmbilicError")

    ## \f$ \vec{T} = Rej(\hat{N})(\vec{s}) \f$
    def tangent(self):
        """Rejection of steepest on normal"""
        return (self.steepest() << self.normal()).hat()

    ## \f$ \hat{B} = \hat{N} \times \hat{T} \f$
    def binormal(self):
        """Binormal"""
        return self.normal() ^ self.tangent()

    ## The Darboux Frame:\n
    # \f$ \hat{x}\hat{T} +
    # \hat{y}\hat{B} +
    # \hat{z}\hat{n} \f$
    def Darboux(self):
        """The Darboux Frame Transformation:

        D(X) == Tangent
        D(Y) == Binormal
        D(Z) == Normal
        """
        from .. import ziprows
        return ziprows(self.tangent(), self.binormal(), self.normal())

    ## \f$ \psi = \tan^{-1}{\frac{s_y}{s_x}} \f$
    def aspect(self):
        """Aspect angle (degrees)"""
        from ..utils.trig import arctand2
        s = self.slope()
        return arctand2(s.imag, s.real)

    ## \f$ \delta = \cos^{-1}{\hat{n}\cdot\hat{z}} \f$
    def dip(self):
        """Dip angle (degrees)"""
        from ..utils.trig import arccosd
        return arccosd(self.normal().z)

    ## root quadrature sum of gradient components
    #  \returns \f$ \sqrt{1+f_x^2+f_y^2} \f$
    def slope3D(self):
        """slope3D is ||gradient||"""
        return reduce(np.hypot, self._grad())

    ## The First Fundamental Form
    #  \returns
    # \f$ {\bf I}=\left[\begin{array}{cc}  E & F \\ F & G \end{array}\right]\f$
    #  \n Breaking it down saves 0.02 seconds on 1 processors (so don't).
    def first_fundamental_form(self):
        """I"""
        return self._make_symmetric_matrix(self.E, self.F, self.G)

    ## No Short Cuts II.
    def second_fundamental_form_pedagogical(self):
        """II"""
        return self._make_symmetric_matrix(self.L, self.M, self.N)

    ## The Second Fundamental Form
    ## \returns
    ## \f${\bf II}=\left[\begin{array}{cc} L & M \\ M & N \end{array}\right]\f$
    ## norm is factored out for speed, and applied via its __rdiv__ method.
    def second_fundamental_form_fast(self):
        """II"""
        return self._make_symmetric_matrix(
            *map(self.slope3D().__rdiv__, (self.L1, self.M1, self.N1))
            )

    ## Set the form here-- various computation styles have different speeds.
    second_fundamental_form = second_fundamental_form_fast

    ## Convert to geo.morphics.gauss.Derivatives data.
    def geomorphics(self):
        """Convert to Derivatives."""
        from .gauss import Derivatives
        return Derivatives(*map(self.__getattribute__, Derivatives._fields))

    plan = wrap_gauss(gauss.plan)
    profile = wrap_gauss(gauss.profile)
    streamline = wrap_gauss(gauss.streamline)
    tangential = wrap_gauss(gauss.tangential)
    total = wrap_gauss(gauss.total)
    gaussian = wrap_gauss(gauss.gaussian)
    mean = wrap_gauss(gauss.mean)


## The Single Quadric: Class represents ONE and only ONE quadric.
class Quadric(Qbase):
    """Quadric(p0,p1,p2,p3,p4,p5=0) (no non-trivial z terms-- so it's really a
    Quadratic"""
    ## With
    ## \f$ z = f(x,y) = p_0 x^2 + p_1 xy + p_2 y^2 + p_3 x + p_4 y + p_5 \f$,
    ## define:\n
    ## \f${\bf P}\equiv \left(
    ## \begin{array}{c} p_0\\p_1\\p_2\\p_3\\p_4\\p_5 \end{array}\right)   \f$
    ## as a matrix
    def __init__(self, p0, p1, p2, p3, p4, p5=0.):
        """Quadric(p0,p1,p2,p3,p4,p5)
        -->
        p0 x**2           (E)
        p1 xy             (2F)
        p2 y**2           (G)
        p3 x
        p4 y
        p5
        """
        ## a matrix of parameters
        self.p = np.matrix([p0, p1, p2, p3, p4, p5]).T

    ## Get an entry from the parameter matrix: p[i] = \f$p_i\f$
    def __getitem__(self, index):
        return self.p[index, 0]

    ## \f$ f(x, y) \f$ \n Evaluate like a function...
    def __call__(self, x, y):
        return super(Quadric, self).__call__(x, y) + self[5]

    ## \f$\frac{\partial f}{\partial z} \equiv -1\f$ (Paraboloid...)
    @property
    def z(self):
        """Singleton z gradient is by definition -1"""
        return -1.

    ## The shape operator\n \f$ {\bf S} = {\bf I}^{-1}{\bf II} \f$.
    def shape_operator(self):
        """S = I**(-1)*II"""
        return self.invI()*self.second_fundamental_form()

    ## \f$ I^{-1} \f$ from I().
    def invI(self):
        """inverse of I"""
        return self.I.I

    ## straight outta numpy: linalg.eig
    def eig(self):
        """shape operator eigenvalues"""
        return linalg.eig(self.S)

    ## sort eig() results into correct order
    ## \returns tuple \f$ \lambda_1, \lambda_2,
    ## {\bf {\hat{v}_1}}, {\bf {\hat{v}_2}} \f$
    def order_eigens(self):
        """order the eigen values of shape operator"""
        (lam1, lam2), matrix_ = self.eig()
        v1, v2 = map(np.ravel, map(np.array, matrix_.T))
        # swap order if needed
        if abs(lam1) < abs(lam2):
            v1, v2 = v2, v1
            lam1, lam2 = lam2, lam1

        return (lam1, lam2), (v1, v2)

    ## Principle Curvature "Vectors" in (x, y, z)?
    def kappa(self):
        """Principal Curvatures"""
        from itertools import starmap, izip
        from operator import mul
        from .. import Vector
        return map(self.Darboux().T,
                   [Vector(item[0], item[1], 0) for item in
                    starmap(mul, izip(*self.order_eigens()))])

    ## Method is polymorphic (in spite of being static), and cannot be a
    # function: (2, 2) matrix
    #  \param xx \f$ f_{xx} \f$
    #  \param xy \f$ f_{xy} = f_{yx}\f$
    #  \param yy \f$ f_{xx} \f$
    #  @returns \f$\left(\begin{array}{cc}
    #  f_{xx} &  f_{xy} \\
    #  f_{yx} &  f_{yy} \end{array}\right).\f$
    @staticmethod
    def _make_symmetric_matrix(xx, xy, yy):
        """_make_symmetric_matrix(xx,xy,yy) --> matrix([[xx,xy],[xy,yy]])"""
        return np.matrix([[xx, xy],
                          [xy, yy]])
