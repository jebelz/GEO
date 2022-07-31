"""See Cylinder.__doc___"""
## \namespace geo.metric.euclid.kant.cylindrical Cylindrical Coordinates

from . import _bases


__all__ = ("Cylinder",)


## Cylindrical Coordinate Vector
# @image html cylinder.gif
class Cylinder(_bases.NonCartesianBase):
    """Cylinder(rho, theta, z)"""

    ## \f$ (\rho, \theta, z) \f$
    components = ('rho', 'theta', 'z')

    ## Don't care what you call \f$\rho\f$
    @property
    def radius(self):
        """radius is rho"""
        return self.rho

    ## for lazy typers
    r = radius

    ## ditto
    p = radius

    ## or \f$\theta\$
    @property
    def phi(self):
        """phi is also theta"""
        return self.theta

    ## Azimuth is phi
    az = phi

    ## _aka_ "z".
    @property
    def zed(self):
        """zed is z."""
        return self.z

    ## Cart to Cyl
    # \param (x, y, z)
    # \returns \f$ \rho = \sqrt{x^2+y^2} \f$
    # \returns \f$ \theta = \tan^{-1}{y/x} \f$
    # \return \f$ z = z\f$
    @staticmethod
    def _fromcart(x, y, z):
        from ....utils.trig import arctan2
        rho = (x**2 + y**2)**0.5
        theta = arctan2(y, x)
        return rho, theta, z

    ## Cylindrical to Cartesian
    # \param \f$ \rho, \theta, z \f$
    # \returns \f$ x = \rho\cos{\theta} \f$
    # \returns \f$ y = \rho\sin{\theta} \f$
    # \return \f$ z = z\f$
    @staticmethod
    def _tocart(rho, theta, z):
        from ....utils.trig import sin, cos
        return rho*cos(theta), rho*sin(theta), z

    ## Contravariant Scaling
    # \param rho
    # \returns \f${\bf \vec{\nabla} r} = \rho {\bf \hat{\rho}} +
    # \frac{\theta}{\rho} {\bf \hat{\theta}} + z {\bf \hat{z}} \f$
    def grad(self, rho):
        """Gradient at rho."""
        return type(self)(self.rho,
                          self.theta/rho,
                          self.z)

    ## Jacobian
    # \returns \f$ \left[\begin{array}{ccc}
    # v_x & -\rho v_y & 0 \\
    # v_y & \rho v_x& 0 \\
    # 0 & 0 & 1 \end{array} \right] \f$
    def jacobian(self):
        """"tensor = v.jacobian()

        is wrt to global bases."""
        from ..tensor import Tensor
        zero = 0 * self.rho  # copy array_like, or not.
        one = 1 + zero
        v = self.vector
        return Tensor(v.x, -self.rho * v.y, zero,
                      v.y, self.rho * v.x, zero,
                      zero, zero, one)

    ## Cartesian Unit Vector
    # \return \f$ \cos{\theta}{\bf \hat{x}} + \sin{\theta}{\bf \hat{y}} \f$
    def rhohat(self):
        """rho unit vector (Cartesian)"""
        from ....utils.trig import sin, cos
        from ..vector import Vector
        return Vector(cos(self.theta),
                      sin(self.theta),
                      0 * self.theta)

    ## Cartesian Unit Vector
    # \return \f$ \hat{\rho}_y {\bf \hat{x}} - \hat{\rho}_x {\bf\hat{y}} \f$
    def thetahat(self):
        """theta unit vector (Cartesian)"""
        from ..vector import Vector
        rho = self.rhohat()
        return Vector(rho.y, -rho.x, rho.z)

    ## Cartesian Unit Vector
    # \return \f$ \hat{\rho}_z {\bf\hat{x}} + \hat{\rho}_z {\bf\hat{y}} +
    # \sqrt{ \hat{\rho}_x^2 +  \hat{\rho}_y^2} {\bf{\hat{z}}} \f$
    def zhat(self):
        """zed unit vector (Cartesian)"""
        from ..vector import Vector
        rho = self.rhohat()
        return Vector(rho.z, rho.z, (rho.x**2 + rho.y**2)**0.5)

    ## <a href="http://en.wikipedia.org/wiki/Line_element">Line Element</a>
    # \param drho
    # \param dtheta
    # \param dz
    # \returns \f$ \delta\rho {\bf\hat{\rho}} +
    # \delta\theta\rho {\bf\hat{\theta}} + \delta z{\bf\hat{z}} \f$
    def line_element(self, drho, dtheta, dz):
        """line element."""
        return (drho * self.rhohat() +
                dtheta * self.rho * self.thetahat() +
                dz * self.zhat())
