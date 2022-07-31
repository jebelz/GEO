"""The ellipsoid as a datim mix-in- this is where Transformations happen."""

## \namespace geo.detic.newton._datum Transformations on an Ellispoid
import abc
import functools

from ...utils.trig import cosd, sind, arcsin, pi, arctan2, arctand2
from ...metric.euclid.kant import Polar
from ...metric.euclid.affine import Affine
from ...metric.euclid.vector import Vector

from .. import origins


## Old School LTP to SCH
# \param  R  Radius of Curavature
# \param  x  LTP x-coordinate
# \param  y  LTP y-coordinate
# \param  z  LTP z-coordinate
# \returns \f$ s = R\arctan{\frac{x}{R+z}}\f$
# \returns \f$ c = R\arcsin{\frac{y}{\sqrt{x^2+y^2+z^2}}}\f$
# \returns \f$ h = \sqrt{x^2+y^2+z^2} - R\f$
def ltp2sch(R, x, y, z):
    """s, c, h = ltp2sch(R, x, y, z)

    Local tangent Cartesian to Local tangent sphere is a function of
    coordinates and the radius of curvature (along x &| s). All inputs
    and outputs have the same units: Length."""
    h = R + z
    rho = abs(Vector(x, y, h)).w #(x**2 + y**2 + h**2)**0.5
    s = R * arctan2(x, h)
    c = R * arcsin(y/rho)
    h = rho - R
    return s, c, h


## Old School SCH to LTP
# \param  R  Radius of Curavature
# \param  s  SCH s-coordinate
# \param  c  SCH c-coordinate
# \param  h  SCH h-coordinate
# \returns \f$ (\vec{r}\cdot\hat{y}, \vec{r}\cdot\hat{z},
# \vec{r}\cdot\hat{x}-R)\f$
def sch2ltp(R, s, c, h):
    """x, y, z = sch2ltp(R, s, c, h)

    Local tangent sphere to Local tangent Cartesian is a function of
    coordinates and the radius of curvature (along x &| s). All inputs
    and outputs have the same units: Length."""
    s_lat = s/R
    c_lat = c/R
    r = h + R
    P = Polar(r, pi/2 - c_lat, s_lat).vector
    return P.y, P.z, P.x-R


## Just a place to put coordinate transforms
class DatumMixin(object):
    """This mixin is a temporary place to put transformations"""

    ## This is a real mix-in.
    __metaclass__ = abc.ABCMeta

    ## Initialized a memoized \f$ \epsilon^2\f$
    e2 = None

    ## Local Radius of Curvature is abstract
    localRad = abc.abstractproperty()

    ## Normal Radius of Curvature is Abstract
    N = abc.abstractproperty()

    ## SCH Helper
    # \param lat Latitude, \f$ \phi \f$.
    # \param lon Longitude, \f$ \lambda \f$
    # \param hdg Heading, \f$ \eta \f$
    # \returns geo.detic.coordinates.SCH Origin
    def create_sch_origin(self, lat, lon, hdg):
        """Return SCH origin at lat, lon, hdg"""
        return self.SCH(0., 0., 0., peg=origins.Peg(lat, lon, hdg))

    ## LTP Helper
    # \param lat Latitude, \f$ \phi \f$.
    # \param lon Longitude, \f$ \lambda \f$
    # \param hdg Heading, \f$ \eta \f$
    # \returns geo.detic.coordinates.LTP Origin
    def create_ltp_origin(self, lat, lon, hdg):
        """Return LTP origin at lat, lon, hdg"""
        return self.create_sch_origin(lat, lon, hdg).ltp()

    # \param lat Latitude, \f$ \phi \f$.
    # \param lon Longitude, \f$ \lambda \f$
    # \returns geo.detic.coordinates.LTP Origin    
    def create_nwu_origin(self, lat, lon):
        return self.create_ltp_origin(lat, lon, hdg=origins.NORTH)

    # \param lat Latitude, \f$ \phi \f$.
    # \param lon Longitude, \f$ \lambda \f$
    # \returns geo.detic.coordinates.LTP Origin    
    def create_enu_origin(self, lat, lon):
        return self.create_ltp_origin(lat, lon, hdg=origins.EAST)
    
    ## Geodetic to ECEF tuple:\n
    # An analytic geodetic coordinates.LLH -> coordinates.ECEF calculation
    # with tuple I/O, \n using method _OblateEllipsoid.N().
    # \param lat Latitude, \f$ \phi \f$.
    # \param lon Longitude, \f$ \lambda \f$
    # \param h Height above ellipsoid, \$ h \f$.
    # \returns ( 
    # \f$ x=(N+h)\cos{\varphi}\cos{\lambda}\f$,\n
    # \f$ y= (N+h)\cos{\varphi}\sin{\lambda}\f$, \n
    # \f$ z= ((1-\epsilon^2)N+h) \sin{\varphi} \f$) \n
    def LatLonHgt2XYZ(self, lat, lon, h):
        """LatLonHgt2XYZ(lat, lon, h) --> (x, y, z)

        lat       is the latitude (deg)
        lon       is the longitude (deg)
        h         is the heigh (m)

        (x, y, z) is a tuple of ECEF coordinates (m)
        """
        N = self.N(lat)
        cos_lat = cosd(lat)
        return (cos_lat * cosd(lon) * (N+h),
                cos_lat * sind(lon) * (N+h),
                sind(lat) * ((1-self.e2)*N + h))

    ## An iterative coordinates.ECEF -> coordinates.LLH  calculation with
    # tuple I/O, \n using ecef2llh()
    # \param x
    # \param y
    # \param z
    # \param iters=10 Number of iterations
    # \returns (\f$ \varphi\f$, \f$\lambda\f$, \f$h\f$)
    def XYZ2LatLonHgt(self, x, y, z, iters=10):
        """XYZ2LatLonHgt(x, y, z {iters=10})--> (lat, lon, hgt)

        calls module function:

        ecef2llh(self, x, y,  [, **kwargs])"""
        return ecef2llh(self, x, y, z, iters=iters)

    ## Get the geo.metric.euclid.vector.Vector from Earf's center to a
    # latitude and longitude on the Ellipsoid
    # \param lat \f$ \phi \f$, degrees OR an orgins.Peg
    # \param lon = None or Peg
    # \param dummy = None
    # \returns Vector from Ceneter to \f$(\phi, \lambda) \f$.
    def center_to_latlon(self, lat, lon=None, dummy=None):
        """center_to_latlon(lat {lon})

        Input:

        lat              is a PegPoint
        lat, lon         are latitude and longitude (degrees)


        Output:

        Vector           instance points from Core() to (lat, lon)
        """
        lat, lon, dummy = self._parse_peg(lat, lon, dummy)
        return self.LLH(lat, lon, 0.).ecef().vector()

    ## TODO: Verify hdg argument changes
    # \param lat \f$ \varphi \f$ or (lat, lon, hdg)
    # \param lon \f$ \lambda \f$
    # \param hdg \f$ \eta \f$
    # \returns affine.Affine from LTP to ECEF
    def ltp2ecef_affine(self, lat, lon=None, hdg=None):
        """"ltp2ecef_affine(lat, lon=None, hdg=None)
        is the Affine transform for LTP to ECEF.

        lat can be a tuple of (lat, lon, hdg), or a signleton."""
        lat, lon, hdg = self._parse_peg(lat, lon, hdg)
        R = origins.ecef2ltp_rotation(lat, lon, -hdg + 180)
        T = self.center_to_latlon(lat, lon)
        return Affine(R, T)

    # \param lat \f$ \varphi \f$ or (lat, lon, hdg)
    # \param lon \f$ \lambda \f$
    # \param hdg \f$ \eta \f$
    # \returns affine.Affine from ECEF to LTP
    def ecef2ltp_affine(self, lat, lon=None, hdg=None):
        """ecef2ltp_affine(lat {lon, hdg})

        Input:

        lat                   is a PegPoint
        lat, lon, hdg         are latitude, longitude, heading (degrees)


        Output:

        Affine  transform from Core to pegged tangent plane.
        """
        lat, lon, hdg = self._parse_peg(lat, lon, hdg)  # is this superfulous?
        return ~(self.ltp2ecef_affine(lat, lon=lon, hdg=hdg))

    ## Compute function from LTP to SCH
    # \param lat \f$ \phi \f$ (or a Peg point).
    # \param lon=None \f$ \phi \f$
    # \returns \f$ _{\rm sch}f_{\rm xyz}(R=R; \ldots) \f$ Partial evaluation
    # of ltp2sch().
    def TangentPlane2TangentSphere(self, lat, lon=None, hdg=None):
        """TangentPlane2TangentSphere(self, lat, lon, hdg)

        Return a function of (x, y, z) that computes transformation from
        tangent plane to (s, c, h) coordiantes for input (lat, lon, hdg).
        """
        lat, lon, hdg = self._parse_peg(lat, lon, hdg)
        return functools.partial(ltp2sch, self.localRad(lat, hdg))

    ## Compute function from SCH to LTP
    # \param lat \f$ \phi \f$ (or a Peg point).
    # \param lon=None \f$ \phi \f$
    # \returns \f$ _{\rm sch}f_{\rm xyz}(R=R; \ldots) \f$ Partial evaluation
    # of sch2ltp().
    def TangentSphere2TangentPlane(self, lat, lon=None, hdg=None):
        """TangentSphere2TangentPlane(self, lat, lon, hdg)

        Return a function of (s, c, h) that computes transformation from
        tangent sphere to (x, y, z) coordinates for input (lat, lon, hdg).
        """
        lat, lon, hdg = self._parse_peg(lat, lon, hdg)
        return functools.partial(sch2ltp, self.localRad(lat, hdg))

    ## Convenience function for parsing arguments that might be a peg point.
    # or a regular tuple
    # \param peg A PegPoint or a latitude
    # \param lon None or a longitude
    # \param hdg None or a heading
    # @retval peg A tuple of peg values
    @staticmethod
    def _parse_peg(peg, lon=None, hdg=None):
        """_parse_peg(peg, lon, hdg) expects that, or:
        _parse_peg(lat, lon, hdg) or
        _parse_peg([lat, lon, hdg], None, None)"""
        return peg if lon is None else (peg, lon, hdg)  # kwd


## 334 Heritage Functions
class FunctionMixin(object):
    """Fortran-style functions attached to an ellipsoid."""

    ## Earth-Centered Earth-Fixed: Change Ellipsoids
    # \param x  array_like x position
    # \param y  array_like y position
    # \param z  array_like z position
    # \kwd ellipsoid None, or another ellipsoid.Ellipsoid
    # \returns list (x, y, z) referenced to ellipsoid or self
    def ecef2ecef(self, x, y, z, ellipsoid=None):
        """ecef2ecef(x, y, z, ellipsoid=None) [changes ellipsoid]

        (of course, (x, y, z) is the same, but the object is not)
        """
        return (self.ECEF(x, y, z) >> (ellipsoid or self)).tolist()

    ## Cartesian to Geodetic
    # \param x  array_like x position
    # \param y  array_like y position
    # \param z  array_like z position
    # \returns list (x, y, z) in LLH
    def ecef2llh(self, x, y, z):
        """ecef2llh(x, y, z)"""
        return self.ECEF(x, y, z).llh().tolist()

    ## Earth-Centered Earth-Fixed to Local Tanget Plane
    # \param x  array_like x position
    # \param y  array_like y position
    # \param z  array_like z position
    # \param peg A peg point tuple
    # \returns list (x, y, z) in LTP
    def ecef2ltp(self, x, y, z, peg):
        """ecef2ltp(x, y, z, peg)
        peg can be a (lat, lon, hdg) tuple or PegPoint"""
        return self.ECEF(x, y, z).ltp(peg=peg).tolist()

    ## Earth-Centered Earth-Fixed to Local Tangent Sphere
    # \param x  array_like x position
    # \param y  array_like y position
    # \param z  array_like z position
    # \param peg A peg point tuple
    # \returns list (x, y, z) in SCH
    def ecef2sch(self, x, y, z, peg):
        """ecef2ltp(x, y, z, peg)
        peg can be a (lat, lon, hdg) tuple or PegPoint"""
        return self.ECEF(x, y, z).sch(peg=peg).tolist()

    ## Geodetic to Earth-Centered Earth-Fixed
    # \param lat  array_like lat position
    # \param lon  array_like lon position
    # \param hgt  array_like hgt position
    # \returns list (lat, lon, hgt) in ECEF
    def llh2ecef(self, lat, lon, hgt):
        """llh2ecef(lat, lon, hgt)"""
        return self.LLH(lat, lon, hgt).ecef().tolist()

    ## Geoetic ellispoid change
    # \param lat  array_like lat position
    # \param lon  array_like lon position
    # \param hgt  array_like hgt position
    # \kwd ellipsoid None, or another ellipsoid.Ellipsoid
    # \returns list (lat, lon, hgt) in LLH, reference to ellipsoid or self
    def llh2llh(self, lat, lon, hgt, ellipsoid=None):
        """llh2llh(lat, lon, hgt, ellipsoid) [changes ellipsoid]"""
        return (self.LLH(lat, lon, hgt) >> (ellipsoid or self)).tolist()

    ## Geodetic to Local Tangent Plane
    # \param lat  array_like lat position
    # \param lon  array_like lon position
    # \param hgt  array_like hgt position
    # \param peg A peg point tuple
    # \returns list (x, y, z) in LTP
    def llh2ltp(self, lat, lon, hgt, peg):
        """llh2ltp(lat, lon, hgt, peg)
        peg can be a (lat, lon, hdg) tuple or PegPoint"""
        return self.LLH(lat, lon, hgt).ltp(peg=peg).tolist()

    ## Geodetic to Local Tangent Sphere
    # \param lat  array_like lat position
    # \param lon  array_like lon position
    # \param hgt  array_like hgt position
    # \param peg A peg point tuple
    # \returns list (x, y, z) in SCH
    def llh2sch(self, lat, lon, hgt, peg):
        """llh2ltp(lat, lon, hgt, peg)
        peg can be a (lat, lon, hdg) tuple or PegPoint"""
        return self.LLH(lat, lon, hgt).sch(peg=peg).tolist()

    ## Local Tangent Plane to EVEF
    # \param x  array_like x position
    # \param y  array_like y position
    # \param z  array_like z position
    # \param peg A peg point tuple
    # \returns list (x, y, z) in ECEF
    def ltp2ecef(self, x, y, z, peg):
        """ltp2ecef(x, y, z, peg)
        peg can be a (lat, lon, hdg) tuple or PegPoint"""
        return self.LTP(x, y, z, peg=peg).ecef().tolist()

    ## Local Tangent Plane to Geodetic
    # \param x  array_like x position
    # \param y  array_like y position
    # \param z  array_like z position
    # \param peg A peg point tuple
    # \returns list (lat, lon, hgt) in LLH
    def ltp2llh(self, x, y, z, peg):
        """ltp2llh(x, y, z, peg)
        peg can be a (lat, lon, hdg) tuple or PegPoint"""
        return self.LTP(x, y, z, peg=peg).llh().tolist()

    ## Local Tangent Plane to another Local Tangent Plane
    # \param x  array_like x position
    # \param y  array_like y position
    # \param z  array_like z position
    # \param peg A peg point tuple
    # \param peg_prime None, or a target peg point tuple
    # \returns list (x, y, z) in ECEF
    def ltp2ltp(self, x, y, z, peg, peg_prime=None):
        """ltp2ltp(x, y, z, peg [, peg_prime=None])
        pegs can be tuples (lat, lon, hdg) or PegPoint objects"""
        #pylint: disable=R0913
        return self.LTP(x, y, z,
                        peg=peg).ltp(peg=peg_prime).tolist()

    ## Local Tangent Plane to (another) Local Tangent Sphere
    # \param x  array_like x position
    # \param y  array_like y position
    # \param z  array_like z position
    # \param peg A peg point tuple
    # \param peg_prime None, or a target peg point tuple
    # \returns list (s, c, h) in SCH
    def ltp2sch(self, x, y, z, peg, peg_prime=None):
        """ltp2sch(x, y, z, peg [, peg_prime=None])
        pegs can be tuples (lat, lon, hdg) or PegPoint objects"""
        #pylint: disable=R0913
        return self.LTP(x, y, z,
                        peg=peg).sch(peg=peg_prime).tolist()

    ## Local Tangent Plane to Earth-Centered Earth-Fixed
    # \param s  array_like s position
    # \param c  array_like c position
    # \param h  array_like h position
    # \param peg A peg point tuple
    # \returns list (x, y, z) in ECEF
    def sch2ecef(self, s, c, h, peg):
        """sch2ecef(s, c, h, peg)
        peg can be a (lat, lon, hdg) tuple or PegPoint"""
        return self.SCH(s, c, h, peg=peg).ecef().tolist()

    ## Local Tangent Plane to Geodetic
    # \param s  array_like s position
    # \param c  array_like c position
    # \param h  array_like h position
    # \param peg A peg point tuple
    # \returns list (lat, lon, hgt) in LLH
    def sch2llh(self, s, c, h, peg):
        """sch2llh(s, c, h, peg)
        peg can be a (lat, lon, hdg) tuple or PegPoint"""
        return self.SCH(s, c, h, peg=peg).llh().tolist()

    ## Local Tangent Sphere to (another) Local Tanget Plane
    # \param s  array_like s position
    # \param c  array_like c position
    # \param h  array_like h position
    # \param peg A peg point tuple
    # \param peg_prime None, or a target peg point tuple
    # \returns list (x, y, z) in LTP
    def sch2ltp(self, s, c, h, peg, peg_prime=None):
        """sch2ltp(s, c, h, peg [, peg_prime=None])
        pegs can be tuples (lat, lon, hdg) or PegPoint objects"""
        #pylint: disable=R0913
        return self.SCH(s, c, h,
                        peg=peg).ltp(peg=peg_prime).tolist()

    ## Local Tangent Sphere to (another) Local Tanget SPhere
    # \param s  array_like s position
    # \param c  array_like c position
    # \param h  array_like h position
    # \param peg A peg point tuple
    # \param peg_prime None, or a target peg point tuple
    # \returns list (s, c, h) in SCH'
    def sch2sch(self, s, c, h, peg, peg_prime=None):
        """sch2ltp(s, c, h, peg [, peg_prime=None])
        pegs can be tuples (lat, lon, hdg) or PegPoint objects"""
        #pylint: disable=R0913
        return self.SCH(s, c, h,
                        peg=peg).sch(peg=peg_prime).tolist()


## An iterative function that converts Earth-Centered Earth-Fixed to LLH,
# given and ellipsoid.Ellipsoid instance \n It's really a method, but is
# broken out, so you can chose a different function without changing the
# class (TBD)
# \param ellipsoid_of_revolution An ellipsoid.Ellipsoid
# \param x ECEF x coordinate (array_like)
# \param y ECEF y coordinate (array_like)
# \param z ECEF z coordinate (array_like)
# \kwd iter=10 Number of iterations
# retval tuple (latitude, longitude, height)
def ecef2llh_iterative(ellipsoid_of_revolution, x, y, z, iters=10):
    """ecef2llh(ellipsoid_of_revolution, x, y, z [, iters=10])-->
    (lat, lon, hgt)

    Input:
    ------
    ellipsoid_of_revolution  an Ellipsoid instance
    x
    y           ECEF coordiantes (singleton, or array, or whatever
    z


    KeyWord:
    --------

    iters      controls the number of iteration in the loop to compute the
                latitude.


    Output:
    -----
    lat       is the latitude (deg)
    lon       is the longitude (deg)
    h         is the heigh (m)
    """

    lon = arctand2(y, x)
    p = (x**2 + y**2)**0.5

    # What is wrong? why is there no 'r'
    #    r = (x**2 + y**2 + z**2)**0.5

    # start with spherical solution..,
    lat = arctand2(p, z)
    #   ##################ITERATION###################################
    while iters:
        RN = ellipsoid_of_revolution.N(lat)
        h = (p/cosd(lat)) - RN
        lat = arctand2(
            z/p,
            1-ellipsoid_of_revolution.e2 * RN / (RN + h))
        iters -= 1
        ##############################################################
    h = p/cosd(lat) - ellipsoid_of_revolution.N(lat)
    return lat, lon, h


## A good to the millimeter level function from Scott Hensely that does not
# use iteration
def ecef2llh_noniterative(ellipsoid_of_revolution, x, y, z):
    """This would be what it is, If the other is too slow."""
    print x, y, z, ellipsoid_of_revolution
    return NotImplemented


## This function gets called by the Ellipsoid's method.
ecef2llh = ecef2llh_iterative


## Mixin for UTM functions
class SecantTransverseMercatorProject(object):
    """Totally untested and in development"""

    ## Function converting latitude, longitude to UTM.
    def toUTM(self):
        """TBD"""
        from ..mercator.kruger import ToUTM
        return ToUTM.fromellipsoid(self)

    ## Function converting UTM to latitude, longitude.
    def fromUTM(self):
        """Inverse of toUTM"""
        return ~(self.toUTM)
