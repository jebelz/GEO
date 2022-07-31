"""Pauli matrices (and higher dimensional irreducible reps) as
components to the geo.Vector angular momentum operator.
The SU2 ABC makes finite dimensional irreducible reps of SU(2) via the
Irr metaclass. The instances of the concrete classes are normalized
eigenstates.

They will be mixed into pure states, which will in turn me mixed into
density matrices (mixed states of the same irr-rep)-- hence allowing
a complete description of polarization via ADJ-- the adjoint representation.


"""
#pylint: disable=R0921

## \namespace geo.desic.hilbert.pauli
# <a href="http://en.wikipedia.org/wiki/Algebra_of_physical_space">
# The Algebra of Physical Space</a> with Spinors

import abc
import operator
import functools

import numpy as np

from ...utils import exceptions
from ...metric.wigner import casimir


## local complex matrix maker.
_matrix = functools.partial(np.matrix, dtype=complex)


## Base Class for Irreducible Representations
class SU2(object):
    """SU2(m)  (degree m)"""

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def j(self):
        """ 'j' is the spin on the rep."""

    ## Instances are \f$J_z\f$ eigenvectors.
    def __init__(self, m):
        if m not in type(self).span():  # raise
            raise exceptions.NilpotencyError("Invalid order: %s" % str(m))
        ## \f$ J_z|j, m\rangle = m|j, m\rangle \f$
        self.m = m

    ## \f$ 2j +1 \f$
    @classmethod
    def dim(cls):
        """Dimension of class."""
        return casimir.dim(cls.j)

    ## \f$ J^2|j, m\rangle = j(j+1)|j, m \rangle \f$
    @classmethod
    def casimir(cls):
        """Casimir invariant."""
        return casimir.invariant(cls.j)

    @classmethod
    def J2(cls):
        """J-squared operator."""
        return (cls.J()**2).real

    ## \f$ \lambda \in [-j, -j+1, \cdots, j-1, j] \f$
    @classmethod
    def span(cls):
        """Array range of allowed 'm' values."""
        return np.linspace(cls.j, -cls.j, cls.dim())

    @classmethod
    def I(cls):
        """Identity operator."""
        return _matrix(np.diag([1.]*(cls.dim())))

    ## \f$ J^{\pm} = \delta_{m(m'\pm1)} \sqrt{j(j+1) - mm'} \f$
    @classmethod
    def ladder(cls, d):
        """ladder(d) for d = +/- 1 is the raising (lowering) operator."""
        l = np.diag([1.]*(cls.dim()-1), d)
        casimir_ = cls.casimir()
        for i, m in enumerate(cls.span()):
            for j, dummy in enumerate(cls.span()):
                if l[i, j]:
                    l[i, j] *= np.sqrt(casimir_ - m*(m-d))
        return _matrix(l)

    ## \f$ \left[ \begin{array}{ccccc}
    # j & 0 & \ldots & \ldots & 0 \\
    # 0 & j-1 &  0 & \vdots &  \\
    # \vdots &  & \ddots &   & \vdots  \\
    # \vdots &  &   &   -(j-1) & \vdots \\
    #  0  &\ldots  & \ldots  &   0 &  -j  \\
    # \end{array} \right] \f$
    @classmethod
    def Jz(cls):
        """Jz operator (diagonal)"""
        return _matrix(np.diag(cls.span()))

    ## \f$ J_x = \frac{1}{\sqrt{2}}(J^+ + J^-)\f$
    @classmethod
    def Jx(cls):
        """Jx operator."""
        return (cls.ladder(1) + cls.ladder(-1))/2

    ## \f$ J_x = \frac{-i}{\sqrt{2}}(J^+ +  J^-)\f$
    @classmethod
    def Jy(cls):
        """Jy operator."""
        return (cls.ladder(1) - cls.ladder(-1))/2j

    ##\f${\bf\vec{J}}=J_x{\bf \hat{x}}+J_y{\bf \hat{y}}+J_y{\bf\hat{z}}\f$
    @classmethod
    def J(cls):
        """Vector 'J' operator."""
        from ... import Vector
        return Vector(cls.Jx(), cls.Jy(), cls.Jz())

    ## \f$ |j, m\rangle \f$
    def __array__(self):
        return _matrix(
            [0]*int(self.j - self.m) + [1] + [0]*int(self.j + self.m)
            ).T

    def __rmul__(self, other):
        return other * self.__array__()

    def __str__(self):
        return str(np.matrix(self))


## Irreducible Representation Factory
class Spin(type):
    """irr = Spin(j)"""

    ## \f$j\f$ defines irreducible rep.
    # \return SU_2<j> A new class just for this rep.
    def __new__(mcs, j):
        if not casimir.isrep(j):
            raise exceptions.CasimirInvariantError("CasimirError")
        return type('SU2_%i' % casimir.dim(j), (SU2,), {'j': j})


## Trivial Rep.
SCALAR = Spin(0)

## Fundamental Representation
SPINOR = Spin(0.5)

## \f$ {\bf \vec{\sigma}} = {\sigma_1}{\bf  \hat{x}} + {\sigma_2}
# {\bf  \hat{y}} + {\sigma_3}{\bf  \hat{z}}  \f$
sigma = 2*SPINOR.J()

## \f$\sigma_1=\left(\begin{array}{cc}  0 & 1 \\ 1 & 0 \end{array} \right) \f$
sigma1 = sigma.x

## \f$\sigma_2=\left(\begin{array}{cc}  0 & -i \\ i & 0 \end{array} \right) \f$
sigma2 = sigma.y

## \f$\sigma_3=\left(\begin{array}{cc}  1 & 0 \\ 0 & -1 \end{array} \right) \f$
sigma3 = sigma.z

## \f$ \sigma_0 = -i \sigma_1 \sigma_2 \sigma_3 \f$
sigma0 = (-1j*reduce(operator.mul, sigma.iter()))


## \f$ | \frac{1}{2}, +\frac{1}{2} \rangle \f$ (Up)
UP = SPINOR(0.5)

## \f$ | \frac{1}{2}, -\frac{1}{2} \rangle \f$ (Down)
DOWN = SPINOR(-0.5)


## Vector Representation
VECTOR = Spin(1)

## 5-dimensional Rep.
TENSOR = Spin(2)
