"""Dirac Matrices in Minkowski Space"""

#pylint: disable=C0103,E1101

## \namespace geo.desic.hilbert.dirac \f$ \gamma^{\mu} \f$ and
# <a href="http://en.wikipedia.org/wiki/Spacetime_algebra">
# The Mother Ship Algebra</a>.
from numpy import kron

from ..minkowski.minkowski import FourVector

from . import pauli


## Spineducible Representation.
class Spin(pauli.Spin):
    """Unfinished meta class extrapolates pauli to dirac"""


## \f$ \gamma^{\mu}_{(2j+1)} \f$
def Gamma(j):
    """Gamma(j) for 2(2j+1) square matrix rep."""
    from operator import mul
    paul = pauli.Spin(j)
    sig = paul.J() * 2 #(2*j+1) # or just 2?
    sig0 = (-1j*reduce(mul, sig.iter()))

    gam0 = kron(pauli.sigma0, sig0)
    gam = FourVector.fromwxyz(gam0, *[kron(1j*pauli.sigma2, s) for s in
                                      sig.iter()])

    return gam


## Fundamental Rep of Dirac Matrices.
gamma = Gamma(0.5)

##  \f$ \gamma^0 = \sigma_3 \otimes \sigma_0 \f$
gamma0 = gamma.t

##  \f$ \gamma^{\mu} = (\gamma^0, i\sigma_2 \otimes {\bf \vec{\sigma}}) \f$
gamma1 = gamma.x
##  \f$ \gamma^{\mu} = (\gamma^0, i\sigma_2 \otimes {\bf \vec{\sigma}}) \f$
gamma2 = gamma.y
##  \f$ \gamma^{\mu} = (\gamma^0, i\sigma_2 \otimes {\bf \vec{\sigma}}) \f$
gamma3 = gamma.z
## \f$ \gamma^5 = i \gamma^0 \gamma^1 \gamma^2 \gamma^3 \f$
gamma5 = (1j*gamma0*gamma1*gamma2*gamma3)
