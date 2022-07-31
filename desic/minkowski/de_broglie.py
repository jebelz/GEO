"""A module for relativistic Doppler/aberration corrections"""
## \namespace geo.desic.minkowski.de_broglie 4-wave-vectors


## Get your radar's 4-wave vector.
# \param nu Radar's center frequency (cps)
# \param n Radar's pointing direction (metric.euclid.vector.Vector)
# \param m [=0] is the effective mass per photon.
# \retval k \f$k_{\mu} = (2\pi\nu, \frac{\nu}{c}{\bf \hat{n}}) \f$
def k(nu, n, m=0.):
    """SI radar frequency, antenna pointing unit vector --> Four k-vector"""
    from ..minkowski import c
    from .minkowski import FourVector

    # normalize
    n = n.hat()

    # get wave 3-vector
    k_ijk = (nu/c) * n

    # Compute Frequency (with effective mass)
    k_0 = k_ijk**2 + m**2

    return FourVector(k_0, k_ijk)
