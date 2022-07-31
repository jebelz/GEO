"""Generate SO3 Rotation Lie Algebra Style"""
## \namespace geo.metric.SO3 Rotation Generators.

from .. import BASIS, Vector

## The Angular Momentum Vector generates rotations
# \f$ {\bf \hat{J}} \f$
J = Vector(*[(1j)*(e.dual()) for e in BASIS])


## \f$ {\bf R} = e^{ -i\theta{ \bf\hat{n}\cdot \hat{J}}} \f$
def R(v, theta):
    """Construct a rotation matrix from exp(i theta J*n)"""
    return (-(1j)*theta*(J*v.hat()).w).exp()


## \f$ A^k_{ij} = [J_j, J_j] \cdot J_k \f$ \n
# The Adjoint Representation of SO(3).
def adjoint():
    """Fascinating."""
    #pylint: disable=E1103
    return list(J.imag.iter())
