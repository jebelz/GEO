"""Peg points and predefined coordinate relations.

There is the PegPoint (aka, Peg):

peg = Peg(lat, lon, hdg)

This defines a point on an (as yet to be defined ellipsoid) at which to
attach a tangent coordinate system.

The heading is defined here as clockwise with respect to North = 0 degrees.

The steps to rotate from ECEF to an LTP are defined here too. They are:

ECEF --> ENU at (lat=0, lon=0)
>>>f = ecef2enu00

(0, 0) to (lat, lon)
>>>g = enu002enu

And then finally, a rotation from ENU to any heading:
>>>h = enu2ltp

Together, they are:

fgh(v) = (f * g * h)(v)

"""
#pylint: disable=E1101

## \namespace geo.detic.origins More like PegPoints

import collections

from ..metric.euclid.tensor import zipcols, Tensor
from ..metric.euler import charts

__all__ = ('NORTH', 'EAST', 'bearing', 'Graticule', 'PegPoint', 'Peg')


## North Heading
NORTH = 0.

## The East Heading Direction
EAST = NORTH + 90


## This takes ECEF to NWU (and visa versa); NWU is the local tangent plane
# coordinate system with zero heading:\n
# \f$ \left[\begin{array}{ccc} 0&0&1\\0&-1&0\\1&0&0 \end{array}\right]\f$
ECEF2NWU = Tensor(0, 0, 1, 0, -1, 0, 1, 0, 0).versor()

## Passive rotation from ECEF to an LTP defined by:
# \param lat \f$ \phi \f$, latitude.
# \param lon \f$ \lambda \f$, longitude.
# \param hdg \f$ \eta \f$, heading.
# \returns \f$ R(\phi, \lambda, \eta) =
#  \left[\begin{array}{ccc} 0&0&1\\0&-1&0\\1&0&0 \end{array}\right]
#  \left[\begin{array}{ccc} 1&0&0\\0&\cos{\lambda}&-\sin{\lambda}\\
# 0&\sin{\lambda}&\cos{\lambda} \end{array}\right]
# \left[\begin{array}{ccc} \cos{\phi}&0&\sin{\phi}\\0&1&0\\-\sin{\phi}&0&
# \cos{\phi}\end{array}\right]
# \left[\begin{array}{ccc} \cos{\eta}&\sin{\eta}&0\\-\sin{\eta}&\cos{\eta}&0\\
# 0&0&1 \end{array}\right] \f$
def ecef2ltp_rotation(lat, lon, hdg):
    """R = ecef2ltp_rotation(lat, lon, hdg) is the PASSIVE (alias) rotation
    from ECEF to the LTP frame defined by the arguments."""
    from operator import mul
    rot = [ECEF2NWU,
           charts.roll(-lon),
           charts.pitch(-lat),
           charts.yaw(hdg)]

    return reduce(mul, rot)


## Bearing Between two Points specified by latitude and longitude
# (http://mathforum.org/library/drmath/view/55417.html)
# \param lat1 starting latitude
# \param lon1 starting longitude
# \param lat2 ending latitude
# \param lon2 ending longitude
# \returns \f$ b(\phi_1, \lambda_1, \phi_2, \lambda_2)= \tan^{-1}{\frac{\sin{
# (\lambda_2-\lambda_1)}\cos{\phi_2}}{(\cos{\phi_1}\sin{\phi_2}-\sin{\phi_1}
# \cos{\phi_2})\cos{(\phi_2-\phi_1)}}} \f$
# @image html bearing.jpeg
# @image latex bearing.jpeg
def bearing(lat1, lon1, lat2, lon2):
    u"""\u03c8 = bearing(\u03d5_1, \u03bb_1, \u03d5_2, \u03bb_2)

    \u03d5_1, \u03bb_1  are the latitude and longitude of the starting point
    \u03d5_2, \u03bb_2  are the latitude and longitude of the ending point

    \u03c8 is the bearing (heading), in degrees, linking the start to the end.
    """
    from ..utils.trig import sind, cosd, arctand2
    dlon = (lon2 - lon1)
    return arctand2(sind(dlon)*cosd(lat2),
                    cosd(lat1)*sind(lat2) - sind(lat1)*cosd(lat2)*cosd(dlon))


## (lat, lon) tuple
_LatLon = collections.namedtuple("_LatLon", "lat lon")


## (lat, lon, hdg) tuple
_PegPoint = collections.namedtuple("_PrePoint", "lat lon hdg")


## A TBD Mixin for ENU
class _ENU(object):

    """Define ENU heading"""
    hdg = EAST


## A TBD Mixin for NED
class _NED(object):

    """Define NED heading"""
    hdg = NORTH


## Lat by Lon pair
class Graticule(_LatLon):
    """Graticule(lat, lon)"""

    ## Get PegPoint for ENU coordinates at same origin
    def enu(self):
        """get an ENU pegpoint"""
        return PegPoint(self.lat, self.lon, _ENU.hdg)


## The Peg Point
class PegPoint(_PegPoint, Graticule):
    """PegPoint(lat, lon, hdg) all in degrees, of course.

    A peg point need not be a PegPoint, as a tuple works in a pinch.

    Coordinates are designed to use singleton-like peg, but might work
    with array-like pegs (use at you own risk).
    """

    def __str__(self):
        return str(tuple(self))

    def __repr__(self):
        return super(PegPoint, self).__repr__().replace('_Pre', 'Peg')

    ## Write peg to a stream
    def tostream(self, stream, mode='a'):
        """private stream writing function"""
        return stream.write(
            {'b': self._binary, 'a': self._ascii}[mode](*tuple(self))
            )

    ## pack args into a double
    @staticmethod
    def _binary(*args):
        """private stream helper"""
        from struct import pack
        return pack('d'*len(args), *args)

    ## return a line
    @staticmethod
    def _ascii(*args):
        """private ascii writer"""
        return " ".join(map(str, args)) + '\n'

    ## Heritage method.
    # /param ellipsoid A geo.detic.newton.ellipsoid.Ellipsoid
    # /returns
    # geo.detic.newton.ellipsoid.Ellipsoid.local_radius_of_curvature
    def radius_of_curvature(self, ellipsoid):
        """Get the radius of curvature on an ellipsoid"""
        return ellipsoid.local_radius_of_curvature(self.lat, self.hdg)


## Peg is PegPoint
Peg = PegPoint
