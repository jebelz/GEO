"""Composite object in Minkowski space. There is no nod to Einstein b/c
the module is not worthy. That would require a full GR implementation."""
## \namespace geo.desic.minkowski <a
# href="http://en.wikipedia.org/wiki/Minkowski_space">Minkowski space</a>
# @image html cone.png
# @image latex cone.png


## <a href="http://en.wikipedia.org/wiki/Speed_of_light">The ultimate scalar
#  </a>
c = 299792458.


## \f$ u_{\mu} = \gamma(c, {\bf \vec{v}}) \f$
#  \param v a geo.metric.euclid.vector.Vector velocity (SI)
#  \returns u a minkowski.FourVector four-velocity
def u(v):
    """v (SI) --> u, the 4-velocity in natural units"""
    from .minkowski import FourVector, ONE
    return gamma(v) * FourVector(ONE, beta(v))


## \f$ {\bf\vec{\beta}} = {\bf\vec{v}}/c \f$
#  \param v a geo.metric.euclid.vector.Vector velocity (SI)
#  \returns beta a geo.metric.euclid.vector.Vector velocity (natural units)
def beta(v):
    """S.I.->Natural Units for 3-velocity"""
    return v / c


## \f$ \gamma = \frac{1}{\sqrt{1-|{\bf \vec{\beta}}|^2}} \f$
#  \param beta a geo.metric.euclid.vector.Vector velocity (natural units)
#  \returns gamma a geo.metric.euclid.scalar.Scalar Lorentz factor
def gamma(beta_):
    """Lorentz factor (natural units)"""
    return (1. - beta_**2) ** (-0.5)


## Boost Helper
# \param v a geo.metric.euclid.vector.Vector velocity (SI)
# \returns \f$ \Lambda({\bf \vec{v}}) = \left[
# \begin{array}{cccc}
# \gamma & -\gamma\beta n_x &  -\gamma\beta n_y &  -\gamma\beta n_z \\
# -\gamma\beta n_x & 1+(\gamma-1)n_x^2 & (\gamma-1)n_xn_y &(\gamma-1)n_xn_z\\
# -\gamma\beta n_y & (\gamma-1)n_yn_x & 1+(\gamma-1)n_y^2 &(\gamma-1)n_yn_z\\
# -\gamma\beta n_z & (\gamma-1)n_zn_x & (\gamma-1)n_zn_y &1+(\gamma-1)n_z^2
# \end{array} \right]
def boost(v):
    """SI 3-velocity to Lorentz transform:

    M = boost(v); v = vx + vy + vz
    """
    from ...metric.euclid.tensor import DELTA
    from .lorentz import Lorentz

    b = beta(v)
    g = gamma(b)
    g_b = -g*b
    return Lorentz(g, g_b, g_b, DELTA + (g-1)*(b & b)/(b*b))



## A pure spacial rotation in __M4__:
# \param tensor_ A rotation matrix masquerading as a rank-2 Euclidean tensor.
# \return \f$ \Lambda \f$ A lorentz.Lorentz transform that is pure spacial
# rotation.
def rotation(tensor_):
    """Compute (1, -1, -1, -1) eta-metric around a rotation:

    (T or R3 x R3 --> Lambda on M4 x M4)"""
    from ...metric.euclid.scalar import ONE
    from ...metric.euclid.vector import NULL
    from .lorentz import Lorentz
    return Lorentz(ONE, NULL, NULL, tensor_)
