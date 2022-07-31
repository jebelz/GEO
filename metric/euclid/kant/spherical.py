"""Kant is the home Non Cartesian representations of Vectors."""
## \namespace geo.metric.euclid.kant.spherical Spherical Coordinates
import operator

from ....utils.trig import TWO_PI, degrees, radians

from . import _bases


__all__ = ('Polar',)


## Spherical Coordinate Vector
# @image html polar.gif
class Polar(_bases.NonCartesianBase):
    """Polar(radius, theta, phi)"""

    ## \f$ (r, \theta, \phi) \f$
    components = ('radius', 'theta', 'phi')

    ## Cone Angle (deg)
    # \return \f$ \theta \f$, in degrees
    # @image html scn2.gif
    def zenith(self):
        """Polar angle in degrees"""
        return degrees(self.theta)

    ## ENU Azimuth, measured from N, clockwise.
    # \returns \f$ 90 - \phi \f$
    # @image html scn2.gif
    def az(self):
        """Azimuth/Clock angle: degrees clockwise from North"""
        return 90 - degrees(self.phi)

    ## Azimuth is clock
    clock = az

    ## Zenith is cone
    cone = zenith

    ## Elevation Angle
    # \return \f$ 90 - \theta \f$
    def el(self):
        """Elevation Angle: degrees about horizon (XY-plane)"""
        return 90 - self.zenith()

    ## ALT is EL.
    alt = el

    ## Angle Pair.
    # \return az(), el()
    def azel(self):
        """az(), el()"""
        return self.az(), self.el()

    ## Construct from az and el
    # \param az()
    # \param el()
    # \returns cls instance.
    @classmethod
    def fromazel(cls, az, el, az0=90, el0=90):
        """Construct a unit vector from azimuth and elevation angles (deg):

        polar = Polar.fromazel(az, el)"""
        return cls(1.+0*(az+el),  # broadcast or not
                   radians(-(el-el0)),
                   radians(-(az-az0)))

    ## Construct from alt and az.
    # \param alt()
    # \param az()
    # \returns cls instance.
    @classmethod
    def fromaltaz(cls, alt, az):
        """polar = Polar.fromaltaz(alt, az) (in degrees)"""
        return cls.fromazel(az, alt)

    ## Cartesian ro Spherical
    # \param  x
    # \param  y
    # \param  z
    # \returns \f$ r = \sqrt{x^2 + y^2 +z^2} \f$
    # \returns \f$ \theta = \cos^{-1}{z/r} \f$
    # \returns \f$ \phi = \tan^{-1}{y/x} \f$
    @staticmethod
    def _fromcart(x, y, z):
        """Simple Cartesian to spherical, with not-simple routine to
        handle null vectors or array_like with null vectors in them."""
        from ....utils.trig import arccos, arctan2
        radius = (x**2 + y**2 + z**2)**0.5
        phi = arctan2(y, x)
        try:  # array v singleton branch
            # ensure non-divide by 0 for singleton
            if radius:  # divide by zero check
                theta = arccos(z/radius)
            else:
                # radius was 0, hence theta is degenerate
                theta = 0.
        except ValueError:
            # we got here b/c x,y,z are array_like.
            scale = 0.*radius
            for count, item in enumerate(radius):
                # invert each on, with a watch on zeros, a EAFF implementation
                # requires catching a numpy Exception, LBYL is pure python.
                if item:  # divide by zero
                    scale[count] = operator.truediv(z[count], item)
                else:
                    scale[count] = 1.
            theta = arccos(scale)
        return radius, theta, phi

    ## Spherical to Cartesian
    # \param radius \f$ r \f$
    # \param Polar angle \f$ \theta \f$
    # \param Azimuth angle \f$ \phi \f$
    # \returns \f$ x = r \sin{\theta}\cos{\phi} \f$
    # \returns \f$ y = r \sin{\theta}\sin{\phi} \f$
    # \returns \f$ z = r \cos{\theta} \f$
    @staticmethod
    def _tocart(radius, theta, phi):
        from ....utils.trig import sin, cos
        return tuple(radius*item for item in
                     (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)))

    ## Look Angles:
    # \param circumference \f$2\pi\f$
    # \returns \f$ (\theta, \phi) \f$
    def look(self, circumference=TWO_PI):
        """Look angles wrt to Zed"""
        return [item * circumference/TWO_PI for item in (self.theta, self.phi)]

    ## Gradient of a scalar functon at point:
    # \param dfdr \f$ \frac{\partial f}{\partial r} \f$
    # \param dfdt \f$ \frac{\partial f}{\partial \theta} \f$
    # \param dfdp \f$ \frac{\partial f}{\partial \phi} \f$
    # \returns
    # \f$ \nabla f = \frac{\partial f}{\partial r} {\bf \hat{r}} +
    # \frac{1}{r}
    #  \frac{\partial f}{\partial \theta}  {\bf \hat{\theta}} +
    # \frac{1}{r\sin{\theta}}
    # \frac{\partial f}{\partial \phi}  {\bf \hat{\phi}} ) \f$
    def grad(self, dfdr, dfdt, dfdp):
        """Gradient..."""
        from ....utils.trig import sin
        return type(self)(dfdr,
                          dfdt / self.radius,
                          dfdp / (self.radius * sin(self.theta)))

    ## Radial unit vector.
    # \returns \f$ {\bf \hat{r}} = {\bf \vec{r}} / ||r||\f$
    def rhat(self):
        """unit vector along r"""
        return self.vector.hat()

    ## Polar Unit vector
    # \returns  \f$ {\bf \hat{\theta}} =
    # {\bf \hat{\phi}} \times{\bf \hat{r} } \f$
    def thetahat(self):
        """unit vector along theta"""
        from ....utils.trig import sin, cos
        from ..vector import Vector
        return Vector(cos(self.theta)*cos(self.phi),
                      cos(self.theta)*sin(self.phi),
                      -sin(self.theta))

    ## Azimuth Unit vector
    # \returns \f$ {\bf \hat{\phi}} = \frac{1}{\sin{\theta}}
    # {\bf \hat{z}} \times {\bf \hat{r}}    \f$
    def phihat(self):
        """Unit vector along phi"""
        from ....utils.trig import sin, cos
        from ..import Vector
        return Vector(-sin(self.phi),
                      cos(self.phi),
                      0)

    ## \f$ ({\bf \hat{r}}, {\bf \hat{\theta}}, {\bf \hat{\phi}}) \f$
    _hats = (rhat, thetahat, phihat)

    ## Tuple of line-elements
    # \returns tuple \f$(1, r, r\sin{\theta}) \f$
    def line(self):
        """polar.line() --> (1, r, r*sin(theta))"""
        from ....utils.trig import sin
        return (1.,
                self.r,
                self.r * sin(self.theta))

    ## <a href="http://en.wikipedia.org/wiki/Line_element">Line Element</a>
    # \param dr \f$ dr \f$
    # \param dtheta \f$ d\theta \f$
    # \param dphi \f$ d\phi \f$
    # \returns tuple \f$(1, r, r\sin{\theta}) \f$
    # \returns \f$ dr {\bf \hat{r}} + rd\theta {\bf \hat{\theta}} +
    # r \sin{\theta} d\phi {\bf \hat{\phi}} \f$
    def line_element(self, dr, dtheta, dphi):
        """vector = line_element(dr, dtheta, dphi)"""
        line = self.line()
        return (self.rhat() * dr * line[0] +
                self.thetahat() * dtheta * line[1] +
                self.phihat() * dphi * line[2])

    ## Jacobian
    # \returns \f$ \left[\begin{array}{ccc}
    # -\hat{\theta}_z\hat{\phi}_y & -r\hat{\theta}_x &
    # -r\hat{\theta}_z\hat{\phi}_x\\
    # \hat{\theta}_z\hat{\phi}_x & r\hat{\theta}_y &
    # -r\hat{\theta}_z\hat{\phi}_x\\
    # \sqrt{\hat{\theta}_x^2+\hat{\theta}_y^2} &
    # r\hat{\theta}_z & 0 \end{array} \right] \f$
    def jacobian(self):
        """tensor = v.jacobian()"""
        from ..tensor import Tensor
        theta = self.thetahat()
        phi = self.phihat()
        return Tensor(-theta.z * phi.y,
                      self.radius * theta.x,
                      -self.radius * theta.z * phi.x,

                      theta.z * phi.x,
                      self.radius * theta.y,
                      -self.radius * theta.z * phi.x,

                      (theta.x**2 + theta.y**2) **0.5,
                      self.radius * theta.z,
                      0. * self.radius)

## Polar is somewhat misnamed, so Sphere is available-- but it is not to be
# confused with spherical vectors.
Sphere = Polar
