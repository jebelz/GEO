"""TBD UTM"""
## \namespace geo.detic.mercator.utm TBD
# <a href="http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system">UTM</a>.

from __future__ import division

from ...utils import trig


## TBD Zone
class Zone(object):
    """Zone"""


## TBD Grid
class Grid(object):
    """Grid"""


## Point Scale Factor
k0 = 0.9996

## \f$N_0\f$ vs sgn(lat)
N0 = {1: 0, -1: 10000}


## \param lat \f$ \phi \f$
#  \param lon \f$ \lambda \f$
#  \returns
#  \f$ \frac{k_0}{\sqrt{1-(\sin{\lambda}\cos{\varphi})^2}} \f$
def point_scale(lat, lon):
    """The point-scale function"""
    return k0/trig.sqrt(1 - (trig.sind(lon)*trig.cosd(lat))**2)


## \param lat \f$ \phi \f$
#  \param lon \f$ \lambda \f$
#  \returns
#  \f$ \tan^{-1}{\cos{\lambda}\sin{\varphi}} \f$
def convergence(lat, lon):
    """The convergence function"""
    return trig.arctan(trig.tand(lon)(trig.sind(lat)))
