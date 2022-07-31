"""van_elfrinkhof

Implements 4 x 4 real quaternion matrices (SO(4) Right) -- SO(4)_L

The bases are numpy arrays derived from the kron product of
pauli spin (1/2) matrices.

W, I, J, K = BASIS

are the SO(4)_R quaternions as 4 x 4 real numpy matricies.

W_L, I_L, J_L, K_L = BASIS_LEFT

are likewise for SO(4)_L.

Projections onto these and other bases are avialble, including
chiral projection on the light cone.
"""
## \namespace geo.metric.euler.van_elfrinkhof \f$ SO(4)_R) \f$
# representation of __H__.
from ...desic.hilbert import pauli

import numpy as np


## \f$ {\bf w}_R = \sigma_0 \otimes \sigma_0 \f$
W = np.kron(pauli.sigma0, pauli.sigma0).real

## \f$ {\bf i}_R = \sigma_1 \otimes i\sigma_2 \f$
I = np.kron(pauli.sigma3, 1j*pauli.sigma2).real

## \f$ {\bf j}_R = i\sigma_2 \otimes \sigma_0 \f$
J = np.kron(1j*pauli.sigma2, pauli.sigma0).real

## \f$ {\bf k}_R = i\sigma_1 \otimes i\sigma_2 \f$
K = np.kron(pauli.sigma1, 1j*pauli.sigma2).real


## \f$ ||q|| = det({\bf q})^{\frac{1}{4}} \f$
def norm(q):
    """Norm of a 3x3 quaternion."""
    from numpy.linalg import det
    return det(q)**(1./4.)


## \f$ \bar{q} = {\bf q}^T \f$
def conj(q):
    """Conjugate is Transpose."""
    return q.T


## Basis of unit Quaternions.
BASIS = (W, I, J, K)


## v -> (0, v) -> 4 x 4 matrix.
def v2q(vector_):
    """Construct 4x4 from a vector."""
    return vector_.quaternion().fourXfour()


## Alias transformation
# \param v Vector
# \param q SO(4) quaternion
# \returns v' Vector
def alias(q, v):
    """Alias transform"""
    return asvector(conj(q) * v2q(v) * q)


## Alias transformation
# \param v Vector
# \param q SO(4) quaternion
# \returns v' Vector
def alibi(q, v):
    """Alibi transform"""
    return asvector(q * v2q(v) * conj(q))


## Scalar part.
# \returns \f$ \frac{1}{4}[q_{00}+q_{11}+q_{22}+q_{33}] \f$
def project_S(q):
    """Project out w component."""
    return q.trace()[0, 0]/4.


## \f$ SO(4)_R \f$ X-projection
# \returns \f$ \frac{1}{4}[(q_{01}-q_{10})-(q_{23}-q_{32})] \f$
def project_X(q):
    """Project out right handed 'i' component."""
    return ((q[0, 1] - q[1, 0]) - (q[2, 3] - q[3, 2]))/4.


## \f$ SO(4)_L \f$ X-projection
# \returns \f$ \frac{1}{4}[(q_{01}-q_{10})+(q_{23}-q_{32})] \f$
def project_XL(q):
    """Project out left handed 'i' component."""
    return ((q[0, 1] - q[1, 0]) + (q[2, 3] - q[3, 2]))/4.


## \f$ SO(4)_R \f$ Y-projection
# \returns \f$ \frac{1}{4}[(q_{02}+q_{13})-(q_{20}+q_{31})] \f$
def project_Y(q):
    """Project out right handed 'j' component."""
    return ((q[0, 2] + q[1, 3]) - (q[2, 0] + q[3, 1]))/4.


## \f$ SO(4)_L \f$ Y-projection
# \returns \f$ \frac{1}{4}[(q_{02}-q_{13})-(q_{20}-q_{31})] \f$
def project_YL(q):
    """Project out left handed 'j' component."""
    return ((q[0, 2] - q[1, 3]) - (q[2, 0] - q[3, 1]))/4.


## \f$ SO(4)_R \f$ Z-projection
# \returns \f$ \frac{1}{4}[(q_{03}-q_{12})+(q_{21}-q_{30})] \f$
def project_Z(q):
    """Project out right handed 'k' component."""
    return ((q[0, 3] - q[1, 2]) + (q[2, 1] - q[3, 0]))/4.


## \f$ SO(4)_L \f$ Z-projection
# \returns \f$ \frac{1}{4}[(q_{03}+q_{12})-(q_{21}+q_{30})] \f$
def project_ZL(q):
    """Project out left handed 'k' component."""
    return ((q[0, 3] + q[1, 2]) - (q[2, 1] + q[3, 0]))/4.


## Chiral Projection: +x light cone
def project_Xv(q):
    """01-10"""
    return (project_X(q) + project_XL(q))/2


## Chiral Projection: +y light cone
def project_Yv(q):
    """02-20"""
    return (project_Y(q) + project_YL(q))/2


## Chiral Projection: +y light cone
def project_Zv(q):
    """03-30"""
    return (project_Z(q) + project_ZL(q))/2


## Chiral Projection: -x cross product
def project_Xa(q):
    """32-23"""
    return (project_X(q) - project_XL(q))/2


## Chiral Projection: -y cross product
def project_Ya(q):
    """13-31"""
    return (project_Y(q) - project_YL(q))/2


## Chiral Projection: -z cross product
def project_Za(q):
    """21-12"""
    return (project_Z(q) - project_ZL(q))/2


## Chiral Projection: +x light cone
def project_iXv(q):
    """01+10"""
    return project_Xv(q) + project_Xv(q.T)


## Chiral Projection: +y light cone
def project_iYv(q):
    """02+20"""
    return project_Yv(q) + project_Yv(q.T)


## Chiral Projection: +y light cone
def project_iZv(q):
    """03+30"""
    return project_Zv(q) + project_Zv(q.T)


## Chiral Projection: -x cross product
def project_iXa(q):
    """32+23"""
    return project_Xa(q) + project_Xa(q.T)


## Chiral Projection: -y cross product
def project_iYa(q):
    """13+31"""
    return project_Ya(q) + project_Ya(q.T)


## Chiral Projection: -z cross product
def project_iZa(q):
    """21+12"""
    return project_Za(q) + project_Za(q.T)


def asvector(q):
    """Extract vector part."""
    from .. import Vector
    return Vector(project_X(q), project_Y(q), project_Z(q))

def asquaternion(q):
    from .. import Scalar
    return Scalar(q.trace()) + asvector(q)



## Convert Right Handed Quaternion to Left Handed
def r2l(q):
    """Reverse Handed-ness."""
    r = 1*q
    r[1, 2] *= -1
    r[1, 3] *= -1
    r[2, 1] *= -1
    r[2, 3] *= -1
    r[3, 1] *= -1
    r[3, 2] *= -1
    return r

## Left Handed Basis.
BASIS_LEFT = map(r2l, BASIS)

W_L, I_L, J_L, K_L = BASIS_LEFT


project_SL = project_S


'''


def associate(q):
    s = project_S(q)
    x = project_X(q)
    y = project_Y(q)
    z = project_Z(q)
    xl = project_XL(q)
    yl = project_YL(q)
    zl = project_ZL(q)

    return matrix([
        [s, -x, -y, -z],
        [-xl, s - 2*(q[0, 0] + q[1, 1]), z - 2*(q[2, 1] - q[3, 0]),
   , 0],
        [-yl, 0,  s - 2*(a[0, 0] + a[2, 2]), 0, 0],
        [-yl, 0, 0, s - 2*(a[0, 0] + a[3, 3])]
    ]
    )'''


## The associative matrix of a 4d rotation is the outer product of 2
# 4-element things.
def svd(self):
    """in progress..."""
    from numpy import linalg
    u, sigma, v = linalg.svd(self.asmatrix())
    return [(s_i*u_i, v_i) for u_i, s_i, v_i in
            zip(self.frommatrix(u).cols(),
                sigma,
                self.frommatrix(v).rows())
            if s_i]  # filter



from .hamilton import Basis
BASIS_RIGHT = _Basis(W, I, J, K)
BASIS_LEFT = _Basis(W_L, I_L, J_L, K_L)
BASIS = BASIS_RIGHT
