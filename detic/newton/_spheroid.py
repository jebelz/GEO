"""The Ellipsoid takes shape"""
## \namespace geo.detic.newton._spheroid Oblate Spheroid Shape
import abc

from ...utils.trig import cosd, sind

from .. import origins


## The <a href=http://en.wikipedia.org/wiki/Right_angle>Right Angle</a>
RIGHT = 90.


## \f$ \frac{1}{2}\f$
#  <a href="http://en.wikipedia.org/wiki/Half">One Half</a>
HALF = 0.5


## A popular radical function
# @param a semi-major axis
# @param b semi-minor axis
# @param lat Latitude (degrees)
# @retval \f$f(a,b;\phi)\equiv\sqrt{(a\cos{\varphi})^2+(b\sin{\varphi})^2}\f$ 
def radical(a, b, lat):
    """for a, b, lat get radical"""
    return ((a * cosd(lat))**2 + (b * sind(lat))**2)**HALF


## compute bearing and distance on ellipsoid:\n
#  Starting with: \n \f$ \lambda^{(n)} = \lambda_2-\lambda_1 \f$ \n
#  iterate: \n
#  \f$ \sin{\sigma} = \sqrt{(\cos{\beta_2}\sin{\lambda})^2+
#  (\cos{\beta_1}\sin{\beta_2}-\sin{\beta_1}\cos{\beta_2}\cos{\lambda})^2}
#  \f$ \n
#  \f$ \cos{\sigma} =
#  \sin{\beta_1}\sin{\beta_2}+\cos{\beta_1}\cos{\beta_2}\cos{\lambda} \f$ \n
#  \f$ \sin{\alpha} =
#  \frac{\cos{\beta_1}\cos{\beta_2}\sin{\lambda}}{\sin{\sigma}} \f$ \n
#  \f$ \cos{2\sigma_m} =
#  \cos{\sigma}-\frac{2\sin{\beta_1}\sin{\beta_2}}{\cos^2{\alpha}} \f$ \n
#  \f$ C = \frac{f}{16}\cos^2{\alpha}[4+f(4-3\cos^2{\alpha})] \f$ \n
#  \f$ \lambda^{(n+1)} =
#  L + (1-C)f\sin{\alpha}
#  (\sigma+C\sin{\sigma}[\cos{2\sigma_m}+C\cos{\sigma}(-1+2\cos^2{2\sigma_m})])
#  \f$ \n \n
#  Then, with: \n
#  \f$ u^2 = \cos^2{\alpha} \frac{a^2-b^2}{b^2} \f$ \n
#  \f$ A = 1 + \frac{u^2}{16384}(4096+u^2[-768+u^2(320-175u^2)]) \f$ \n
#  \f$ B = \frac{u^2}{1024}(256-u^2[-128+u^2(74-47u^2)]) \f$ \n
#  \f$ s = bA(\sigma-\Delta\sigma) \f$ \n
#  \f$ \Delta\sigma =
#  B\sin{\sigma\big[\cos{2\sigma_m}+
#  \frac{1}{4}B[\cos{\sigma}(-1+2\cos^2{2\sigma_m})-
#  \frac{1}{6}B\cos{2\sigma_m}(-3+4\sin^2{\sigma})
#  (-3+4\cos^2{2\sigma_m})]\big]} \f$ \n
#  see http://www.movable-type.co.uk/scripts/latlong-vincenty.html
#  @param lat1 starting latitude
#  @param lon1 starting longitude
#  @param lat2 ending latitude
#  @param lon2 ending longitude
#  @retval s between start and end, along ellipsoid surface
#  @retval alpha1 heading at start
#  @retval alpha2 heading at end
def _great_circle(self, lat1, lon1, lat2, lon2):
    """s, alpha1, alpha2 = great_circle(lat1, lon1, lat2, lon2)

    (lat1, lon1)    p1's location
    (lat2, lon2)    p2's location

    s               distance along great circle
    alpha1          heading at p1
    alpha2          heading at p2
    """
    from ...utils.trig import arctan2, radians, pi
    #pylint: disable=R0914
    phi1, L1, phi2, L2 = lat1, lon1, lat2, lon2

    a = self.a
    f = self.f
    b = (1-f)*a
    U1 = self.common2reduced(phi1)  # aka beta1
    U2 = self.common2reduced(phi2)
    L = L2-L1

    lam = L

    delta_lam = 100000.

    while abs(delta_lam) > 1.e-10:
        sin_sigma = (
            (cosd(U2)*sind(lam))**2 +
            (cosd(U1)*sind(U2) - sind(U1)*cosd(U2)*cosd(lam))**2
        )**HALF
        cos_sigma = sind(U1)*sind(U2) + cosd(U1)*cosd(U2)*cosd(lam)
        sigma = arctan2(sin_sigma, cos_sigma)

        sin_alpha = cosd(U1)*cosd(U2)*sind(lam)/sin_sigma
        cos2_alpha = 1-sin_alpha**2

        cos_2sigma_m = cos_sigma - 2*sind(U1)*sind(U2)/cos2_alpha

        C = (f/16.) * cos2_alpha*(4.+f*(4-3*cos2_alpha))

        lam_new = (
            radians(L) +
            (1-C)*f*sin_alpha*(
                sigma +
                C*sin_sigma*(
                    cos_2sigma_m +
                    C*cos_sigma*(
                        -1+2*cos_2sigma_m**2
                        )
                    )
                )
            )

        lam_new *= 180/pi

        delta_lam = lam_new-lam
        lam = lam_new
        continue

    u2 = cos2_alpha * (a**2 - b**2)/b**2

    A_ = 1 + u2/16384*(4096+u2*(-768+u2*(320-175*u2)))
    B_ = u2/1024*(256+u2*(-128+u2*(74-47*u2)))

    delta_sigma = B_ * sin_sigma * (
        cos_2sigma_m - (1/4.) * B_ * (cos_sigma * (-1 + 2 * cos_2sigma_m**2)) -
        (1/6.) * B_*cos_2sigma_m * (-3 + 4 * sin_sigma**2) *
        (-3 + 4 * cos_2sigma_m**2)
        )

    s = b*A_*(sigma-delta_sigma)

    alpha_1 = 180*arctan2(
        cosd(U2)*sind(lam),
        cosd(U1)*sind(U2)-sind(U1)*cosd(U2)*cosd(lam)
        )/pi

    alpha_2 = 180*arctan2(
        cosd(U1)*sind(lam),
        -sind(U1)*cosd(U2)+cosd(U1)*sind(U2)*cosd(lam)
        )/pi

    return s, alpha_1, alpha_2


## A 2-parameter Oblate Ellipsoid of Revolution
class OblateEllipsoidMixin(object):
    """This class is a 2 parameter oblate ellipsoid of revolution and serves 2
    purposes:

    __init__  defines the shape

    all the other methods compute the myriad of "other" ellipsoid parameters in
    the literature, and they provide conversion from the common latitude to all
    the other latitudes out there.
    """

    ## This is a mix-in, which could stand on its own, but doesn't.
    __metaclass__ = abc.ABCMeta


    ## Construct from Axes
    # \param a Semi major axis
    # \param b=None a, or semi minor axis
    # \param c=None b, triaxial is not impemented- no asteroids.
    # \return Ellipsoid
    # \throws ValueError if triaxial.
    @classmethod
    def from_axes(cls, a, b=None, c=None, model=""):
        """Build from axes:

        >>>elipsoid = from_axes(a, b=None, c=None, model="")"""
        if c is not None:  # raise
            raise ValueError("Triaxial not Implemented")
        elif b is None:
            from sys import float_info
            finv = float_info.max
        elif bool(a):
            try:
                finv = a / (a - b)
            except ZeroDivisionError:
                return cls.from_axes(a, model=model)
        return cls(a, finv, model=model)

    ## This __init__ defines the shape
    #  @param a semi-major axis
    #  @param finv Inverse flattening
    #  @param model="Unknown" or some name
    def __init__(self, a, finv, model="Unknown"):
        ## The semi-major axes
        self.a = a
        ## The inverse flattening
        self.finv = finv
        ## The semi-minor axis
        self.b = a*(finv-1)/finv
        ## \f$ \epsilon^2 = 1-(\frac{b}{a})^2\f$ \n The first eccentricity
        #  squared.
        self.e2 = (1.-(self.b/self.a)**2)
        ## The model name
        self.model = str(model)

    ## Test equality of parameters-- I use _OblateEllipsoidMixin.a and
    #  _OblateEllipsoidMixin.b for simplicity
    def __eq__(self, other):
        return (self.a == other.a) and (self.b == other.b)

    ## Test inequality of parameters-- I use Ellipsoid.a and Ellipsoid.b for
    #  simplicity
    def __ne__(self, other):
        return (self.a != other.a) or (self.b != other.b)

    ## \f$\cos{oe} = b/a \f$ \n Cosine of the
    #  <a href="http://en.wikipedia.org/wiki/Angular_eccentricity">
    #  angular eccentricity</a>.
    @property
    def cosOE(self):
        """cos(OE) --> b/a"""
        return self.b/self.a

    ## \f$ f=1-\cos{oe} \f$ \n Flattening
    @property
    def f(self):
        """Flattening"""
        return 1./self.finv

    ## \f$ f' \equiv n = \tan^2{\frac{oe}{2}} = \frac{a-b}{a+b} \f$\n
    #  <a href=" http://en.wikipedia.org/wiki/Flattening">
    #  The  Second  Flattening</a> .
    @property
    def f_prime(self):
        """2nd Flattening"""
        return self.f/(1.+self.cosOE)

    ## \f$ n = \frac{a-b}{a+b} \f$ \n Third Flattening
    @property
    def f_double_prime(self):
        """3rd Flattening"""
        return (self.a-self.b)/(self.a+self.b)

    ## \f$n \equiv f'' \f$
    n = f_double_prime

    ## \f$ \epsilon = \sqrt{1-{\frac{b}{a}}^2}  \f$\n First Eccentricity
    @property
    def e(self):
        """1st Eccentricity"""
        return self.e2**HALF

    ## \f$ \epsilon' = \sqrt{{\frac{a}{b}}^2-1}  \f$\n Second Eccentricity
    @property
    def e_prime(self):
        """2nd Eccentricity"""
        return self.e_prime_squared**HALF

    ## \f$ \epsilon'^2 = \cos{OE}^2-1 \f$\n Second Eccentricity Squared
    @property
    def e_prime_squared(self):
        """1st Eccentricity Squared"""
        return self.cosOE**2 - 1.

    ## \f$ e'' = \sqrt{\frac{a^2-b^2}{a^2+b^2}} \f$ \n Third Eccentricity
    @property
    def e_double_prime(self):
        """3rd Eccentricity"""
        return self.e_double_prime_squared**HALF

    ## \f$ e''^2 = \frac{a^2-b^2}{a^2+b^2}  \f$ \n Third Eccentricity Squared
    @property
    def e_double_prime_squared(self):
        """3rd Eccentricity Squared"""
        return (self.a**2-self.b**2) / (self.a**2+self.b**2)

    ##  \f$ R_1 \f$, \n
    #  <a href="http://en.wikipedia.org/wiki/Earth_radius#Mean_radius:_R1">
    #  Mean Radius</a>.
    @property
    def R1(self):
        """Mean Radius"""
        return (2.*self.a + self.b)/3.

    ##  \f$ R_2 \f$, \n
    #  <a href="http://en.wikipedia.org/wiki/Earth_radius#Authalic_radius:_R2">
    #  Authalic Radius</a>
    @property
    def R2(self):
        """Authalic Radius"""
        #pylint: disable=R0201
        return NotImplemented

    ## \f$ R_3 \f$ \n
    #  <a href=
    #  "http://en.wikipedia.org/wiki/Earth_radius#Volumetric_radius:_R3">
    #  Volumetric Radius</a>
    @property
    def R3(self):
        """Volumetric Radius"""
        return (self.b*self.a**2)**(1./3.)

    ## \f$\eta' = 1/\sqrt{1-\epsilon^2\sin^2{\varphi}} \f$;\n
    #  <a href="http://en.wikipedia.org/wiki/Latitude#Elliptic_parameters">
    #  Inverse of the principle elliptic integrand</a>
    #  @param lat Latitude (degrees)
    #  @retval eta_prime Inverse elliptic integrand
    def eta_prime(self, lat):
        """Inverse of the principle elliptic integrand (lat/deg)"""
        return (1 - (self.e/sind(lat))*2)**(-HALF)

    ## \f$ \frac{\pi}{180^{\circ}} M(\phi) \f$ \n
    #  <a href="http://en.wikipedia.org/wiki/Latitude#Degree_length">Latitude
    #  degree length</a>
    #  \param lat Latitude (degrees)
    #  \param h=0 correct for off-shell inquiries
    #  \retval lat-like  meters per degree latitude
    def latitude_degree_length(self, lat, h=0):
        """Length of a degree of latitude (deg-->m) """
        from ...utils.trig import radians
        r = self.meridional_radius_of_curvature(lat)
        return radians(r+h)

    ## \f$ \frac{\pi}{180^{\circ}} \cos{(\varphi)} N(\phi) \f$\n
    #  <a href="http://en.wikipedia.org/wiki/Latitude#Degree_length">
    #  Longitude degree length</a>
    #  @param lat Latitude (degrees)
    #  \param h=0 correct for off-shell inquiries
    #  @retval lat-like  meters per degree longitude
    def longitude_degree_length(self, lat, h=0):
        """Length of a degree of longitude (deg-->m) """
        from ...utils.trig import radians
        r = self.normal_radius_of_curvature(lat)
        return radians(cosd(lat) * (r+h))

    ## \f$M=M(\phi)=
    #  \frac{(ab)^2}{[(a\cos{\varphi})^2+(b\sin{\varphi})^2]^{\frac{3}{2}}}\f$
    #  @param lat Latitude (degrees)
    #  @retval lat-like  meridional radius of curvature
    def meridional_radius_of_curvature(self, lat):
        """North Radius (northRad): Meridional radius of curvature (M), meters
        for latitude in degrees """
        return (self.a * self.b)**2 / radical(self.a, self.b, lat)**3

    ##  \f$N(\phi) = a \eta' \f$, \n
    #  <a href="http://en.wikipedia.org/wiki/Earth_radius#Normal">Normal Radius
    #  of Curvature</a>
    #  @param lat Latitude (degrees)
    #  @retval lat-like  normal radius of curvature
    def normal_radius_of_curvature(self, lat):
        """East Radius (eastRad): Normal radius of curvature (N), meters for
        latitude in degrees """
        return self.a**2 / radical(self.a, self.b, lat)

    ## Synonym for::normal_radius_curvature
    eastRad = normal_radius_of_curvature

    ## \f$\frac{1}{R(\phi, \alpha)}=\frac{\cos^2{\alpha}}{M(\phi)}+
    #  \frac{\sin^2{\alpha}}{N(\phi)}\f$\n
    #  Radius of curvature along a bearing.
    #  @param lat Latitude (degrees)
    #  @param hdg Heading (degrees)
    #  @retval lat-like  radius of curvature along hdg
    def local_radius_of_curvature(self, lat, hdg):
        """local_radius_of_curvature(lat, hdg)"""
        return 1. / (
            cosd(hdg)**2/self.M(lat) +
            sind(hdg)**2/self.N(lat)
            )

    ## camelCase synonym is for compatibility with 334 standards.
    localRad = local_radius_of_curvature

    ## \f$ N(\phi) \f$ \n Normal Radius of Curvature
    N = normal_radius_of_curvature

    ## \f$ M(\phi) \f$ \n Meridional Radius of Curvature
    M = meridional_radius_of_curvature

    ## Synonym for ::meridional_radius_curvature
    northRad = M

    ## \f$ R=R(\phi)=
    #  \sqrt{\frac{(a^2\cos{\varphi})^2+(b^2\sin{\varphi})^2}
    #  {(a\cos{\varphi})^2+
    #  (b\sin{\varphi})^2}}\f$ \n
    #  Radius at a given geodetic latitude.
    #  @param lat Latitude (degrees)
    #  @retval R Radius (m)
    def R(self, lat):
        """R(lat):
        Radius of Curvature at geodetic latitude, lat."""
        return (
            ((self.a**2 * cosd(lat))**2 + (self.b**2 * sind(lat))**2) /
            radical(self.a, self.b, lat)
            )

    ##\f$m(\varphi)= a (1-e^2) \int_0^{\varphi}{(1-e^2\sin^2{x})^
    #  {-\frac{3}{2}}dx}\f$
    #  @param phi The radians latitude phi()
    #  @param ellipeinc=None Elliptic integral per
    #  <a href="http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipeinc.html#scipy.special.ellipeinc">
    #  default's signature</a>.
    #  @retval m The elliptic parameter at \f$ \varphi \f$.
    def m(self, phi, ellipeinc=None):
        """m(phi, [ellipeinc=None]):

        ellipeinc --> scipy.special.ellipeinc unless keyword is set (that is,
        supply your own incomplete elliptic integral."""
        from ...utils.trig import sin, cos
        if ellipeinc is None:  # kwd
            from scipy.special import ellipeinc
        return (
            ellipeinc(phi, self.e2) -
            self.e2*sin(phi)*cos(phi)/(1-self.e2*sin(phi)**2)**HALF
            )

    ## \f$ \chi(\phi) = 2 \tan^{-1} [(\frac{1+\sin{\varphi}}{1-\sin{\varphi}})
    #  (\frac{1-\sin{oe}\sin{\varphi}}{1+\sin{oe}\sin{\varphi}})^{\sin{oe}}]^
    #  {\frac{1}{2}} -\frac{\pi}{2} \f$ \n
    #  @param lat Latitude (degrees)
    #  @retval lat'
    #  <a href="http://en.wikipedia.org/wiki/Latitude#Conformal_latitude">
    #   Conformal latitude, \f$\chi\f$</a>
    def common2conformal(self, lat):
        """Convert common latitude (deg) to conformal latitude (deg) """
        from ...utils.trig import arctand
        sinoe = self.e2**HALF
        sinphi = sind(lat)
        return 2.*arctand(
            (
                (1.+sinphi)/(1.-sinphi) * (
                    (1.-sinphi*sinoe)/(1.+sinphi*sinoe))**sinoe
            ) ** HALF) - RIGHT

    ## \f$ \beta(\phi) = \tan^{-1}{\sqrt{1-e^2}\tan{\varphi}} \f$ \n
    #  <a href="http://en.wikipedia.org/wiki/Latitude#Reduced_latitude">
    #  Reduced Latitude</a>
    #  @param lat Latitude (degrees)
    #  @retval lat' reduced latitude
    def common2reduced(self, lat):
        """Convert common latitude (deg) to reduced latitude (deg) """
        from ...utils.trig import arctand, tand
        return arctand(self.cosOE * tand(lat))

    common2parametric = common2reduced

    ## \f$q(\varphi)=\frac{(1-e^2)\sin{\varphi}}{1-e^2\sin^2{\varphi}}=
    #  \frac{1-e^2}{2e}\log{\frac{1-e\sin{\varphi}}{1+e\sin{\varphi}}} \f$
    #  @param phi Latitude angle (rad)
    #  @retval q Elliptic nome
    def q(self, phi):
        """q(phi) is the elliptic parameter"""
        from ...utils.trig import log, sin
        sinphi = sin(phi)
        return (
            (1 - self.e2) * sinphi / (1 - self.e2*sinphi**2) -
            ((1 - self.e2)/2/self.e) * log(
                (1 - self.e * sinphi)/(1 + self.e * sinphi))
            )

    ## \f$ \varphi \equiv \frac{\pi}{180} \phi \f$ \n
    #  Latitude in radians, __NB__: symbolic variation is used throughout.
    #  @param lat Latitude, \f$\phi\f$  (deg)
    #  @retval phi Latitude angle, \f$ \varphi \f$, (rad)
    @staticmethod
    def phi(lat):
        """from trig import radians as phi
        does exactly wht this function does"""
        from ...utils.trig import radians
        return radians(lat)

    ## \f$ \xi = \sin^{-1}{\frac{q(\varphi)}{q(\pi/2)}} \f$ \n
    #  <a href="http://en.wikipedia.org/wiki/Latitude#Authalic_latitude)">
    #  Authalic Latitude</a>
    #  @param lat Latitude (degrees)
    #  @retval lat' authalic latitude
    def common2authalic(self, lat):
        """Convert common latitude (deg) to authalic latitude (deg) """
        from ...utils.trig import degrees, arcsin, pi
        phi_ = self.phi(lat)
        return degrees(arcsin(self.q(phi_)/self.q(pi/2.)))

    ## \f$ \psi(\phi) = \tan^{-1}{(1-e^2)\tan{\varphi}} \f$ \n
    #  <a href="http://en.wikipedia.org/wiki/Latitude#Geocentric_latitude">
    #  Geocentric Latitude</a>.
    #  @param lat Latitude (degrees)
    #  @retval lat' geocentric latitude
    def common2geocentric(self, lat):
        """Convert common latitude (deg) to geocentric latitude (deg) """
        from ...utils.trig import arctand, tand
        return arctand(tand(lat) * self.cosOE**2)

    ## \f$ \mu{\phi} = \frac{\pi}{2}\frac{m(\varphi)}{m(\varphi/2)} \f$ \n
    #  <a href="http://en.wikipedia.org/wiki/Latitude#Rectifying_latitude">
    #  rectifying latitude</a>
    #  @param lat Latitude (degrees)
    #  @retval lat' rectifying latitude
    def common2rectifying(self, lat, ellipeinc=None):
        """Convert common latitude (deg) to rectifying latitude (deg) """
        from ...utils.trig import radians
        return RIGHT*self.m(radians(lat))/self.m(radians(RIGHT),
                                                 ellipeinc=ellipeinc)

    ## \f$ \psi(\phi)= \sinh^{-1}{(\tan{\varphi}) -
    #  \epsilon\tanh^{-1}{(\epsilon\sin{\varphi})}}\f$\n
    #  <a href=" http://en.wikipedia.org/wiki/Latitude#Isometric_latitude">
    #  isometric latitude </a>
    #  @param lat Latitude (degrees)
    #  @retval lat' isometric latitude
    def common2isometric(self, lat):
        """Convert common latitude to isometric latitude."""
        from ...utils.trig import log, tan, sin, pi, degrees
        phi_ = self.phi(lat)
        sinphi = sin(phi_)
        return degrees(
            log(tan(pi/4. + phi_/2.)) +
            (self.e/2.)*log((1 - self.e*sinphi)/(1 + self.e*sinphi))
            )

    ## Geodetic latitude is the latitude.
    #  @param lat Latitude (degrees)
    #  @retval lat geodetic latitude
    @staticmethod
    def common2geodetic(lat):
        """trivial...."""
        return lat

    ## See origins.bearing()
    @staticmethod
    def bearing(lat1, lon1, lat2, lon2):
        """hdg = bearing(lat1, lon1, lat2, lon2)

        see origins.bearing
        """
        return origins.bearing(lat1, lon1, lat2, lon2)

    ## get c vs hdg, fit it and find correct zero-- make a peg point.
    def sch_from_to(self, lat1, lon1, lat2, lon2):
        """sch_from_to(lat1, lon1, lat2, lon2)

        computes a peg point at (lat1, lon1) such that
        (lat2, lon2) lies on the s-axis (c=0).

        Under Construction
        """
        #pylint: disable=W0612
        import numpy as np
        hdg, c, h = self._make_hdg_c(lat1, lon1, lat2, lon2)
        psi1 = float(np.poly1d(np.polyfit(hdg, c, 1)).roots)
        psi2 = np.poly1d(np.polyfit(hdg, c, 4)).roots

        iarg = np.argmin(abs(psi2-float(psi1)))

        psi = psi2[iarg]

        print psi-psi1
        return origins.PegPoint(lat1, lon1, float(psi))

    ## get cross track coordinate vs. heading for +/-n degrees around
    #  spherical bearing
    def _make_hdg_c(self, lat1, lon1, lat2, lon2, n=1.):
        """private, and under development...."""
        #pylint: disable=R0913
        import numpy as np
        c = []
        h = []

        p2 = self.LLH(lat2, lon2, 0.)

        b = self.bearing(lat1, lon1, lat2, lon2)

        hdg = np.arange(b-n, b+n, 0.001)
        for theta in hdg:
            peg = origins.PegPoint(lat1, lon1, float(theta))
            p = p2.sch(peg)
            c.append(p.c)
            h.append(p.h)
            continue
        c, h = map(np.array, (c, h))
        return hdg, c, h

    ## Distance between 2 points ON A SPHERE.
    ## @param lat1 starting latitude
    ## @param lon1 starting longitude
    ## @param lat2 ending latitude
    ## @param lon2 ending longitude
    ## @retval s arc length on approximating sphere
    def distance_spherical(self, lat1, lon1, lat2, lon2):
        """d = distance(lat1, lon1, lat2, lon2)"""
        from ...utils.trig import arctan2
        llh1 = self.LLH(lat1, lon1, 0.*lat1)
        llh2 = self.LLH(lat2, lon2, 0.*lat2)

        n1 = llh1.n_vector()
        n2 = llh2.n_vector()

        delta_sigma = arctan2((abs(n1 ^ n2)).w, (n1*n2).w)
        return delta_sigma*self.R1

    ## Developmental....
    ## @param lat1 starting latitude
    ## @param lon1 starting longitude
    ## @param lat2 ending latitude
    ## @param lon2 ending longitude
    ## @retval p2  TBS
    def distance_sch(self, lat1, lon1, lat2, lon2):
        """distance_sch(lat1, lon1, lat2, lon2)
        is under development"""
        hdg = self.bearing(lat1, lon1, lat2, lon2)
        peg = origins.PegPoint(lat1, lon1, hdg)
        p2 = self.LLH(lat2, lon2, 0.).sch(peg)
        return p2

    ## Distance on surface
    ## @param lat1 starting latitude
    ## @param lon1 starting longitude
    ## @param lat2 ending latitude
    ## @param lon2 ending longitude
    ## @retval s between start and end, along ellipsoid surface
    ## @retval alpha1 heading at start
    ## @retval alpha2 heading at end
    def great_circle(self, lat1, lon1, lat2, lon2):
        """See ellipsoid._great_circle"""
        return _great_circle(self, lat1, lon1, lat2, lon2)

    ## True distance between 2 points (with ellipsoid confinement)
    ## @param lat1 starting latitude
    ## @param lon1 starting longitude
    ## @param lat2 ending latitude
    ## @param lon2 ending longitude
    ## @retval s between start and end, along ellipsoid surface
    def distance_true(self, lat1, lon1, lat2, lon2):
        """See ellipsoid._great_circle"""
        return self.great_circle(lat1, lon1, lat2, lon2)[0]

    ## Compute the bearing between
    ## @param lat1 starting latitude
    ## @param lon1 starting longitude
    ## @param lat2 ending latitude
    ## @param lon2 ending longitude
    ## @retval alpha1 heading at start
    ## @retval alpha2 heading at end
    def bearings(self, lat1, lon1, lat2, lon2):
        """See ellipsoid._great_circle"""
        return self.great_circle(lat1, lon1, lat2, lon2)[1:]

    ## Default distance function.
    distance = distance_true
