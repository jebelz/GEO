"""Monge Patch

patch = MongePatch(x, y, z [, weight=ONE])

patch.quadric() is the surface fit at (x=0, y=0).
"""

## \namespace geo.morphic.monge Monge Patch Sampling
import abc

import numpy as np

from ..utils.misc import INV_6X6_FAST

## Real Data Type
FTYPE = np.float64

## Complex Data Type (Precision should match FTYPE):
CTYPE = np.complex128

## Number of Degrees of Freedom
DOF = 6

## Floating point 1
ONE = FTYPE(1)

## Complex 1j.
I = CTYPE(1j)

## local look up for fast square function
SQUARE = np.square

## using __slots__ is a little faster.
MOMENTS = ('w', 'x', 'y', 'z', 'x2', 'y2', 'xy', 'y3', 'xy2', 'x2y', 'x3',
         'y4', 'xy3', 'x4', 'zy', 'zx', 'zy2', 'zxy', 'zx2', 'x2y2', 'x3y')

## Empty 6 x 1 matrix (float64).
SIX_BY_ONE = np.matrix(np.empty(DOF).astype(np.float64))

## A numpy 6 x 6 matrix to use over and over again; type is preserved.
SIX_BY_SIX = np.outer(*list((SIX_BY_ONE,)*2))


## MCS Makes MongePatch, and fills in computed attributes (never do
# work in the constructor, and this pattern avoids that).
class mongepatch(type):
    """Metaclass makes the MongePatch class, and its instances.

    The purpose of this class is to create the MongePatch instance
    with both its coordinates and their MOMENTS. The client does not
    need to know about MOMENTS, but the instances need them--and they
    need to compute them once as instance attributes. As it is an
    antipattern to do work in a constructor, and setting instance
    attributes outside of the constructor is well, an ultimate python sin,
    you do it in the class's metaclass's __call__ method.

    The metaclass also creates static class MOMENTS attributes--as
    creating instance attributes outside of __init__ where there are
    no deafults, is not a good idea."""
    
    ## Create MongePatch class with static momments set to None
    def __new__(mcs, *args, **kwargs):
        cls = super(mongepatch, mcs).__new__(mcs, *args, **kwargs)
        for attr in MOMENTS:
            setattr(cls, attr, None)
        return cls

    ## Create MongePatch instance with instance moments set.
    # \param x
    # \param y
    # \param z
    # \param weight ONE
    # \returns MongePatch instance, fully constructed.
    def __call__(mcs, x, y, z, weight=ONE):
        return mcs + super(mongepatch, mcs).__call__(x, y, z, weight=weight)

    ## Add moments instance attributes to an new instance (and return it).
    # \param mcs (implicit) The metaclass
    # \param self An unfinished MongePatch instance
    # \returns self as a FULLY CONSTRUCTED MongePatch instance.
    # \sideeffect Sets all the instance attributes in ::MOMENTS.
    def __add__(mcs, self):
        ## \f$ x^2 \f$, \n x**2 is 5% faster than x*x, but np.square is best
        self.x2 = SQUARE(self.x)
        ## \f$ y^2 \f$
        self.y2 = SQUARE(self.y)
        ## \f$ xy \f$
        self.xy = self.x*self.y
        ## \f$ y^3 \f$
        self.y3 = self.y2*self.y
        ## \f$ xy^2 \f$
        self.xy2 = self.x*self.y2
        ## \f$ x^2y \f$
        self.x2y = self.x2*self.y
        ## \f$ x^3 \f$
        self.x3 = self.x2*self.x
        ## \f$ y^4 \f$
        self.y4 = self.y3*self.y
        ## \f$ xy^3 \f$
        self.xy3 = self.x*self.y3
        ## \f$ x^2y^2 \f$
        self.x2y2 = SQUARE(self.xy)
        ## \f$ x^3y \f$
        self.x3y = self.x3*self.y
        ## \f$ x^4 \f$
        self.x4 = SQUARE(self.x2)
        ## \f$ zy \f$
        self.zy = self.z*self.y
        ## \f$ zx \f$
        self.zx = self.z*self.x
        ## \f$ zy^2 \f$
        self.zy2 = self.z*self.y2
        ## \f$ zxy \f$
        self.zxy = self.z*self.xy
        ## \f$ zx^2 \f$
        self.zx2 = self.z*self.x2
        return self
    

## The MongePatch is a bunch of points representing a suface,
class MongePatch(object):
    """patch = MongePatch(x, y, z[ ,weight=1])

    holds arrays of x, y, z points (with weights). A local fit of
    a quadratic is:

    q = patch.quadric()
    
    Curvatures can be derived from the quadratic's coefficients.
    The slopes *are* the cooefficents.
    """

    __metaclass__ = mongepatch
    
    ## Construct from a vector of points.
    # \param vector_ A vector
    # \param weight =ONE
    # \returns MongePatch
    @classmethod
    def fromvector(cls, vector_, weight=ONE):
        return cls(vector_.x, vector_.y, vector_.z, weight=weight)

    ## (x, y, z, [weight= ::ONE ]) are Cartesian.
    # \param x x coordinate of postings
    # \param y y coordinate of postings
    # \param z z coordinate at postings
    # \param weight array of weights for z
    def __init__(self, x, y, z, weight=ONE):
        ## x-postings for the surface
        self.x = x
        ## y-postings for the surface
        self.y = y
        ## z-posting for the surface
        self.z = z
        ## weights of the z coordinate.
        self.w = weight

    ## Number of points in the surface (same as numpy len)
    def __len__(self):
        return len(self.z)

    ## numpy.size of arrays.
    @property
    def size(self):
        """forward property to z"""
        return self.z.size

    ## See vandermonde().
    # \returns vandermonde matrix (array) if called by np.matrix (array).
    def __array__(self):
        """see vandermonde, return type is late-binding determined"""
        return self.vandermonde()

    ## Call vandermonde() to make:\n
    ## | X4()    X3Y()     X2Y2()    X3()    X2Y()   X2()  |\n
    ## | X3Y()  X2Y2() XY3()  X2Y()  XY2()  XY()  |\n
    ## | X2Y2()  XY3() Y4()  XY2()  Y3()  Y2()  |\n
    ## | X3()  X2Y() XY2()  X2()  XY()  X()  |\n
    ## | X2Y()  XY2() Y3()  XY()  Y2()  Y()  |\n
    ## | X2()  XY()  Y2()  X()  Y()  ONE()  |, \n
    # \returns matrix
    ##\f$ (\sum_{i=1}^m{\bf{Q}_i{\bf Q}_i^T}) = \left(\begin{array}{cccccc}
    ## s(x^4) & s(x^3y) & s(x^2y^2) & s(x^3) & s(x^2y) & s(x^2) \\ s(x^3y) &
    ## s(x^2y^2) & s(xy^3) & s(x^2y) & s(xy^2) & s(xy) \\ s(x^2y^2) & s(xy^3) &
    ## s(y^4) & s(xy^2) & s(y^3) & s(y^2) \\ s(x^3) & s(x^2y) & s(xy^2) &
    ## s(x^2) & s(xy) & s(x) \\ s(x^2y) & s(xy^2) & s(y^3) & s(xy) & s(y^2) &
    ## s(y)\\s(x^2) & s(xy) & s(y^2) & s(x) & s(y) & s(1)\end{array}\right)\f$.
    def vandermonde(self, m=SIX_BY_SIX):
        """construct symmetric 6x6 matrix from 15 unique elements (in order)

        -- kwarg is not to be used... it's a python
        trick to put an existing matrix in the namespace (for speed, and it
        makes a 1000 sec difference on 1 apple CPU.
        """
        # Use multiple assignment (where required) to speed it up (that is
        # compute the properties one and only once).
        m[0, 0] = self.X4
        m[0, 1] = m[1, 0] = self.X3Y
        m[0, 2] = m[2, 0] = m[1, 1] = self.X2Y2
        m[0, 3] = m[3, 0] = self.X3
        m[0, 4] = m[4, 0] = m[3, 1] = m[1, 3] = self.X2Y
        m[0, 5] = m[5, 0] = m[3, 3] = self.X2

        m[1, 2] = m[2, 1] = self.XY3
        m[1, 4] = m[4, 1] = m[2, 3] = m[3, 2] = self.XY2
        m[1, 5] = m[5, 1] = m[3, 4] = m[4, 3] = self.XY

        m[2, 2] = self.Y4
        m[2, 4] = m[4, 2] = self.Y3
        m[2, 5] = m[5, 2] = m[4, 4] = self.Y2

        m[3, 5] = m[5, 3] = self.X

        m[4, 5] = m[5, 4] = self.Y
        m[5, 5] = self.ONE

        # cast to package precision
        return m.astype(FTYPE)

    ## Solve for Quadric() parameters via.
    # \returns matrix
    #\f$(\sum_{i=1}^m{\bf{Q}_i{\bf Q}_i^T}){\bf P}=
    #\sum_{i=1}^m{z_i{\bf Q}_i}\f$ via:\n\n
    #\f${\bf P}=
    #(\sum_{i=1}^m{\bf{Q}_i{\bf Q}_i^T})^{-1} \sum_{i=1}^m{z_i{\bf Q}_i}\f$\n
    # using __array__() and row() methods.
    def __invert__(self):
        """~surface solves for a Quadric surface such that:

        q = Quadric(*~surface)

        is the quadric fitting z = Q(x=0, y=0)
        """
        from numpy import linalg
        return self.row()*(linalg.inv(self.vandermonde()))

    ##  Transpose row().
    # \returns matrix
    # \f$ \sum_{i=1}^m{z_i{\bf Q}_i}
    # \left(\begin{array}{c} s(zx^2)\\s(zxy)\\s(zy^2)\\s(zx)\\s(zy)\\s(z)
    # \end{array} \right) \f$ \n
    # not used in processing.
    def col(self):
        """row.T"""
        return self.row().T

    ## Make a row matrix from: [ZX2(), ZXY(), ZY2(), ZX(), ZY(), Z()].
    # \returns
    # that is: \n \f$  (\sum_{i=1}^m{z_i{\bf Q}_i})^T \f$ \n
    def row(self, m=SIX_BY_ONE):
        """matrix of regressands-- kwarg is not o be used... it's a python
        trip to put an existing matrix in the namespace (for speed, and it
        makes a 240 sec difference on 1 apple CPU"""
        ## Put the computed elements right into a pre-existing array
        m[:] = (self.ZX2,
                self.ZXY,
                self.ZY2,
                self.ZX,
                self.ZY,
                self.Z)
        return m.astype(FTYPE)

    ## Convert solution (__invert__()) to a Quadric object.
    def quadric(self):
        """invert and make a quadric"""
        from .taylor import Quadric
        p = (~self).T
        ## this is 8 time faster than *-magic  unpacking matrix.tolist()
        ## but it doesn't seem to matter... (0 seconds)-not understood...
        return Quadric(p[0, 0], p[1, 0], p[2, 0], p[3, 0], p[4, 0], p[5, 0])

    ## Return quadric() and size.
    # \returns tuple  quadric(), ::Surface.size
    def __call__(self):
        return self.quadric(), self.size

    ##   \f$ \sum_{i=1}^{i=m}{z_i x_i^2 w_i} \f$
    @property
    def ZX2(self):
        """sum of zx2 """
        return (self.w * self.zx2).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{z_i x_i y_i w_i} \f$
    @property
    def ZXY(self):
        """sum of zxy """
        return (self.w * self.zxy).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{z_i y_i^2 w_i} \f$
    @property
    def ZY2(self):
        """sum of zy2 """
        return (self.w * self.zy2).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{z_i x_i w_i} \f$
    @property
    def ZX(self):
        """sum of zx """
        return (self.w * self.zx).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{z_i y_i w_i} \f$
    @property
    def ZY(self):
        """sum of zy """
        return (self.w * self.zy).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{z_i w_i} \f$
    @property
    def Z(self):
        """sum of z """
        return (self.w * self.z).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x^3_i w_i} \f$
    @property
    def X3(self):
        """sum of x3"""
        return (self.w * self.x3).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x^2_i w_i} \f$
    @property
    def X2(self):
        """sum of x2"""
        return (self.w * self.x2).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x_i w_i} \f$
    @property
    def X(self):
        """sum of x"""
        return (self.w * self.x).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{y^3_i w_i} \f$
    @property
    def Y3(self):
        """sum of y3"""
        return (self.w * self.y3).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{y^2_i w_i} \f$
    @property
    def Y2(self):
        """sum of y2"""
        return (self.w * self.y2).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{y_i w_i} \f$
    @property
    def Y(self):
        """sum of y"""
        return (self.w * self.y).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x^4_i w_i} \f$
    @property
    def X4(self):
        """sum of x4"""
        return (self.w * self.x4).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x^3_i y_i w_i} \f$
    @property
    def X3Y(self):
        """sum of x3y"""
        return (self.w * self.x3y).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x^2_i y_i^2 w_i} \f$
    @property
    def X2Y2(self):
        """sum of x2y2"""
        return (self.w * self.x2y2).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x_i y^3_i w_i} \f$
    @property
    def XY3(self):
        """sum of xy3"""
        return (self.w * self.xy3).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{y^4_i w_i} \f$
    @property
    def Y4(self):
        """sum of y4"""
        return (self.w * self.y4).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x_i y_i^2 w_i} \f$
    @property
    def XY2(self):
        """sum of xy2"""
        return (self.w * self.xy2).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x^2_i y_i w_i} \f$
    @property
    def X2Y(self):
        """sum of x2y"""
        return (self.w * self.x2y).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{x_i y_i w_i} \f$
    @property
    def XY(self):
        """sum of xy"""
        return (self.w * self.xy).ravel().sum()

    ##   \f$ \sum_{i=1}^{i=m}{w_i} \f$
    @property
    def ONE(self):
        """weighted sum of 1"""
        try:
            # always OK in operational mode, not-so in development mode
            len(self.w)
        except TypeError:
            # If weight=1 was a python scalar, then you use this branch:
            result = FTYPE(self.size)
        else:
            result = self.w.sum()
        return result
