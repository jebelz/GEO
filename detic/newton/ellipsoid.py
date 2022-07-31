"""Ellipsoids do a lot. They provide a foundation for c0ordinate classes, they
hold transformation functions, affine computations, and distance computations,
latitude conversions,
"""
## \namespace geo.detic.newton.ellipsoid Ellipsoid with concrete
# geo.detic.coordinates._Coordinate
from .. import coordinates
from .. import origins

from . import _spheroid
from . import _datum

__all__ = ('Ellipsoid',)


CLASSES = {}


## A OblateEllipsoid with coordinate system sub-classes attached
class Ellipsoid(_spheroid.OblateEllipsoidMixin,
                _datum.DatumMixin, _datum.FunctionMixin,
                _datum.SecantTransverseMercatorProject):
    """ellipsoid = Ellipsoid(a, finv [, model="Unknown"])

    a       singleton is the semi-major axis
    finv    is the inverse flattening.
    model   string name of the ellipsoid

    The Ellipsoid is an oblate ellipsoid of revolution, and is a lot more than
    2-parameters

    See __init__.__doc__ for more.
    """
    ## The Ellipsoid needs access to the origins.PegPoint object
    PegPoint = origins.PegPoint

    ## ECEF coordinate instance on Ellispoid
    # \param x
    # \param y
    # \param z
    # \returns geo.detic.coordinates.ECEF
    def ECEF(self, x, y, z):
        """ECEF(x, y, z)"""
        return coordinates.ECEF(x, y, z, ellipsoid=self)

    ## LLH coordinate instance on Ellispoid
    # \param x
    # \param y
    # \param z
    # \returns geo.detic.coordinates.LLH
    def LLH(self, lat, lon, height):
        """LLH(lat, lon, height)"""
        return coordinates.LLH(lat, lon, height, ellipsoid=self)
    
    ## LTP coordinate instance on Ellispoid
    # \param x
    # \param y
    # \param z
    # \returns geo.detic.coordinates.ECEF
    def LTP(self, x, y, z, peg):
        """LTP(x, y, z, peg)"""
        return coordinates.LTP(x, y, z, peg=peg, ellipsoid=self)

    ## ECEF coordinate instance on Ellispoid
    # \param x
    # \param y
    # \param z
    # \returns geo.detic.coordinates.ECEF
    def SCH(self, s, c, h, peg):
        """SCH(s, c, h, peg)"""
        return coordinates.SCH(s, c, h, peg=peg, ellipsoid=self)
    
    def __str__(self):
        return (
            self.model + " ellipsoid\n" +
            "semi-major axis = %f\n" % (self.a,) +
            "inverse flattening = %f\n" % (self.finv,)
            )

    def __repr__(self):
        return "{}(a={}, b={}, model='{}')".format(
            type(self).__name__, self.a, self.b, self.model)


    ## One siderial Day
    day = 86160
    
    ## Angular Velocity for a SIDERIAL day
    # \param peg=None Optional Peg
    # \return \f$ \bf \vec{\Omega} \f$ in ECEF or in peg's LTP.
    def angular_velocity(self, peg=None):
        """Angular Velocity of Earth/planetoid"""
        from math import pi
        from ...metric.euclid.vector import Vector
        omega = Vector(0., 0., 2*pi/self.day)
        if peg:
            ltp = self.create_ltp_origin(*peg)
            omega = (~(ltp.affine2ecef.rotation))(omega)
        return omega

    def coriolis(self, lat, lon, hdg, velocity):
        return -2*(self.angular_velocity([lat, lon, hgt]) ^ velocity)

    def centrifugal(self, lat, lon, hdg, h):
        pass
