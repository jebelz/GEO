"""see Parabola.__doc__"""
## \namespace geo.metric.euclid.kant.parabolic Parabolic coordinates
from . import _bases

__all__ = ("Parabola",)


## <a href="http://mathworld.wolfram.com/ParabolicCoordinates.html">
# Parabolic Coordinate Vector</a>
# @image html parabola.png
class Parabola(_bases.NonCartesianBase):
    """Parabola(u, v, theta)"""

    ## \f$ (u, v, \theta) \f$
    components = ('u', 'v', 'theta')

    ## Conversion from Cartesian
    # \returns \f$ u = \sqrt{x^2+y^2+z^2} + z \f$
    # \returns \f$ v = \sqrt{x^2+y^2+z^2} - z \f$
    # \returns \f$ \theta = \tan^{-1}{y/x} \f$
    @staticmethod
    def _fromcart(x, y, z):
        from operator import add, sub
        from ....utils.trig import arctan2
        radius = (x**2 + y**2 + z**2)**0.5
        u, v = [func(radius, z)**0.5 for func in (add, sub)]
        theta = arctan2(y, x)
        return u, v, theta

    ## Conversion to Cartesian
    # \returns \f$ x = uv\cos{\theta} \f$
    # \returns \f$ y = uv\sin{\theta} \f$
    # \returns \f$ z = \frac{1}{2}(u^2 +v^2) \f$
    @staticmethod
    def _tocart(u, v, theta):
        from ....utils.trig import sin, cos
        z = 0.5*(u**2 - v**2)
        uv = u*v
        return uv*cos(theta), uv*sin(theta), z
