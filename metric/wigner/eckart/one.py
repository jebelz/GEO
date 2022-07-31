"""Pure rank-1 Conversion Formulae:

This is the defining (fundamental) representation of SO(3) versus the
Adjoint (Real, Cartesian) representation.

Weyl Module  dimW Irreps (j; q)         
Tableau                                 
=====================================================================
[0]            3     3   (1; 0)                 v = v


To select the stength of m-eigenstates from a Cartesian vector:

m    function
--------------
+1    V11
+0    V10
-1    V1m1


To get the Cartesian components from a spherical vector:
tensor:


e    function
--------------
x    Vx
y    Vy
z    Vz

Spherical Basis Vectors in Cartesian form are:

q    CONSTANT
------------------------
+1    EPLUS      BASIS[1]
0     EZERO      BASIS[0]
-1    EMINUS     BASIS[-1]

The definition of spherical vector operations are:

op       function       result
--------------------------------
a * b    dot(a, b)      spherical scalar
a ^ b    cross(a, b)    spherical vector
a & b    outer(a, b)    spherical tensor.
"""
## \namespace geo.metric.wigner.eckart.one Fundamental and Adjoint
# representation of SO(3).

import operator
import functools

from geo.metric.euclid.vector import Vector as Cart

## \f$ \sqrt{2} \f$
ROOT2 = 2**0.5


## \f$ \sqrt{3} \f$
ROOT3 = 3**0.5


## \f$ \sqrt{6} \f$
ROOT6 = ROOT2 * ROOT3


## The Scalar Norm codifies the great scalar problem; that is: tradition
# says the  scalar is normalized by 1/3, while racah.py indicates that
# 1/sqrt(3) or even -1/sqrt(3) is correct-- on cannot make the outer product
# self-consitent without addressing this; this way, racah is satisfied,
# and DELTA has a norm for sqrt(3)-- which makes sense as the inchoherent sum
# of 3 normalized states.
SCALAR_NORM = ROOT3


## ## Project m=0 from a Cartesian Vector.
# \param v A Cartesian Vector
# \returns \f$ v^0 = v_z \f$
def V00(v):
    """V00(v) --> |1, 0>

    for v a Cartesian vector."""
    return v.z


## Project m=1 from a Cartesian Vector.
# \param v A Cartesian Vector
# \returns \f$ v^+ = -\sqrt{\frac{1}{2}}(v_x - iv_y) \f$
def V11(v):
    """V11(v) --> |1, 1>

    for v a Cartesian vector."""
    return (-v.x + 1j * v.y) / ROOT2


## Project m=-1 from a Cartesian Vector.
# \param v A Cartesian Vector
# \returns \f$ v^+ = +\sqrt{\frac{1}{2}}(v_x - iv_y) \f$
def V1m1(v):
    """V11(v) --> |1, -1>

    for v a Cartesian vector."""
    return (v.x + 1j * v.y) / ROOT2 
    

## Extract Cartesian 'x' component
# \param v A Spherical Vector
# \returns \f$ v_x = \frac{1}{\sqrt{2}}[-v_+ + v_-] \f$
def Vx(v):
    """v.x = Vx(v) For a Spherical Vector, v."""
    return (-v[1] + v[-1]) / ROOT2


## Extract Cartesian 'y' component
# \param v A Spherical Vector
# \returns \f$ v_y = i\frac{1}{\sqrt{2}}[v_+ + v_-] \f$
def Vy(v):
    """v.y = Vy(v) For a Spherical Vector, v."""
    return 1j * (-v[1] - v[-1]) / ROOT2


## Extract Cartesian 'z' component
# \param v A Spherical Vector
# \returns \f$ v_z = v^0 \f$
def Vz(v):
    """v.z = Vz(v) For a Spherical Vector, v."""
    return v[0]


## \f$ {\bf \hat{e}}^+ = \frac{1}{\sqrt{2}}[{\bf \hat{x}} - {\bf \hat{y}}]\f$\n
# The overall choice of a minus sign is a result of the Condon-Shortley
# phase convention.
EPLUS = Cart(-1, -1j, 0) / ROOT2


## \f$ {\bf \hat{e}}^- = \frac{1}{\sqrt{2}}[{\bf \hat{x}} + {\bf \hat{y}}]\f$
EMINUS = Cart(1, -1j, 0)/ ROOT2
    

## \f$ {\bf \hat{e}}^0 = {\bf \hat{Z}} \f$
EZERO = Cart(0, 0, 1)


## Spherical BASIS[m] --> \f$ {\bf \hat{e}}^{(q=m)} \f$
BASIS = [EZERO, EPLUS, EMINUS]


## A dictionary of  "m" state projecting functions.
PROJECTIONS = {1: V11, 0: V00, -1: V1m1}



## Spherical form of dot product
# \param u Spherical Vector
# \param v Spherical Vector
# \returns \f$ \begin{array}{cccc}
#  & J=0 & J=1 & J=2 \\
# M=2&- &- &-  \\
# M=1&- &- & - \\
# M=0& -u^+v_- + u^0v_0 - u^-v_+&- &-  \\
# M=-1&- &- &- \\
# M=-2&- &- &- \\
# \end{array} \f$
def rdot(u, v):
    """s = dot(u, v)

    This implementation does not take complex conjugates, rather like
    angular momentum addition, it ensures m + m' = 0 = M. Thus, m=+/-1
    basis states are NULL."""
    from . import Scalar
    return Scalar(-u[1]*v[-1] + u[0]*v[0] - u[-1]*v[1])


## Complex inner product space's dot product
# \param u Spherical Vector
# \param v Spherical Vector
# \returns \f$ \begin{array}{cccc}
#  & J=0 & J=1 & J=2 \\
# M=2&- &- &-  \\
# M=1&- &- & - \\
# M=0& u_{-}+v_{-}^{*} + u^{0}v_{0}^{*} + u^{+}v_{-}^{*} + & -  & -  \\
# M=-1&- &- &- \\
# M=-2&- &- &- \\
# \end{array} \f$
def cdot(u, v):
    """s = dot(u, v) with complex conjugation"""
    from . import Scalar
    return Scalar(u[-1]*v[-1].conjugate() +
                  u[0]*v[0].conjugate() +
                  u[1]*v[1].conjugate())


## This is important: why choose a real dot product? It's not positive
# definite, but, it conserves J_z: azimuth quantum number, so that
# all the terms in a scalar have q + q' = 0.
dot = rdot


## Spherical form of cross product
# \param u Spherical Vector
# \param v Spherical Vector
# \returns \f$ \begin{array}{cccc}
#  & J=0 & J=1 & J=2 \\
# M=2&- &- & - \\
# M=1&- & -i(u^+v^0 - u^0v^+) & - \\
# M=0&- & -i(u^+v^- - u^-v^+) & - \\
# M=-1&- & +i(u^-v^0 - u^0v^+) &- \\
# M=-2&- &- &- \\
# \end{array} \f$
def cross(a, b):
    """Cross Product: like dot, there is no complex conjugation. Note
    that the cross product introduces a phase of (1j) when compared to
    a Lie bracket [i, j] = ij-ji = EPSILON.ijk used in the Cartesian
    basis-- they essentially different thanks to 1j."""
    from . import Vector
    return 1j* Vector((a[1]*b[-1] - a[-1]*b[1]),
                      (a[1]*b[0] - a[0]*b[1]),
                      (a[0]*b[-1] - a[-1]*b[0]))


## Spherical form of outer product
# \param u Spherical Vector
# \param v Spherical Vector
# \returns \f$ \begin{array}{cccc}
#  & J=0 & J=1 & J=2 \\
# M=2& - &- & u^+v^+\\
# M=1&- & \sqrt{\frac{1}{2}}(u^+v^0 - u^0v^+) &
#    \sqrt{\frac{1}{2}}(u^+v^0+u^0v^+) \\
# M=0& \sqrt{\frac{1}{3}}(u^+v^- - u^0v^0 + u^-v^+)&
# \sqrt{\frac{1}{2}}(u^+v^- - u^-v^+)&
# \sqrt{\frac{1}{6}}(u^+v^- + 2u^0v^0 + u^-v^+) \\
# M=-1&- & -\sqrt{\frac{1}{2}}(u^-v^0 - u^0v^+) &
#  \sqrt{\frac{1}{2}}(u^-v^0+u^0v^-) \\
# M=-2&- &- & u^-v^-\\
# \end{array} \f$
# This can be verified using racah.VECTOR
def outer(u, v):
    """T = outer(u, v) is a spherical outer prod."""
    from . import Tensor

    # why it is not -1/root(3) is a mystery.
    scalar_ = (u * v) / SCALAR_NORM

    # The 1j phase factor is explained in cross.__doc__
    vector_ = -1j * (u ^ v) / ROOT2

    # This part can be computed from racah.py
    tensor_ = [(u[1]*v[-1] + 2*u[0]*v[0] + u[-1]*v[1]) / ROOT6,
               (u[1]*v[0] + u[0]*v[1]) / ROOT2,
               u[1]*v[1],
               u[-1]*v[-1],
               (u[-1]*v[0] + u[0]*v[-1]) / ROOT2]

    return Tensor(scalar_[0], vector_, tensor_)




