"""Totally in progress, and not checked. Compute VSH from operations
on spherical coordinate tensors.

Is this part of representaiton theory or EM? Both, but it's here
for now.
"""
## \namespace geo.magnetic.legendre <a
# href="http://en.wikipedia.org/wiki/Vector_spherical_harmonics">
# Vector Spherical Harmonics </a>

from scipy.special import sph_harm

from ..utils.trig import tan, exp, sqrt
from ..metric.euclid.kant import Polar


## \f$ f(z) = e^{iz} \f$
def expi(z):
    """exp(iz)"""
    return exp((1j)*z)


## \f$ \cot{\theta} \f$
def cot(theta):
    """Cotangent (rad)"""
    return 1/tan(theta)


## Array_like unity: 1, or ones(len(arg))
def one(*args):
    """Makes ones like arg"""
    result = 1.
    for item in args:
        result = result + 0.*item
        continue
    return result


## (Vector) Spherical Harmonics
class Ylm(object):
    """Ylm(l, m)(theta, phi)

    The derivatives are:

    Ylm.dtheta
    Ylm.phi

    
    Ylm.grad


    The VSH are:
    Ylm.Y
    Ylm.Psi
    Ylm.Phi
    """
    ## \param l degree
    # \param m oreder
    def __init__(self, l, m):
        ## order
        self.l = int(l)
        ## azimuth
        self.m = int(m)

    ## String
    def __str__(self):
        return "  {}\nY\n  {}".format(self.l, self.m)

    ## unicode, with theta and phi.
    def __unicode__(self):
        return u"  {}\nY(\u0398, \u03c6)\n  {}".format(self.l, self.m)
        
    ## \f$ ||m|| \le \f$ checks validity of l, m.
    def __nonzero__(self):
        return abs(self.m) <= self.l

    ## Raising Operator:
    # \param dm \f$ \Delta m \f$
    # \returns \f$ Y_l^{m + \Delta m} \f$
    def __rshift___(self, dm):
        return type(self)(self.l, self.m + dm)

    ## Loweing Operator
    # \param dm \f$ \Delta m \f$
    # \returns \f$ Y_l^{m - \Delta m} \f$
    def __lshift___(self, dm):
        return self >> (-dm)

    ## Straight Function Wwrapper
    # \param theta Polar angle \f$ \theta \f$
    # \param phi Azimuth angle \f$ \phi \f$
    # \returns \f$ Y_l^m(\theta, \phi) \f$
    def __call__(self, theta, phi):
        return sph_harm(self.m,
                        self.l,
                        theta,
                        phi) if self else 0*(theta+phi)

    ## Partial derivative wrt to polar angle
    # \param theta Polar angle \f$ \theta \f$
    # \param phi Azimuth angle \f$ \phi \f$
    # \returns \f$ \frac{\partial Y_l^m(\theta, \phi)}{\partial \theta} =
    # m \cot{(\theta)}Y_l^m(\theta, \phi) + \sqrt{(l-m)(l-m-1)} e^{im\phi}
    # Y_l^{m+1}(\theta, \phi) \f$
    def dtheta(self, theta, phi):
        """dtheta(theta, phi)
        Derivative wrt to theta"""
        n, m = list(self)
        return (
            m * cot(theta) * self(theta, phi) +
            sqrt((n-m)*(n-m-1)) * expi(-phi) * (self >> 1)(theta, phi)
            )

    ## Partial derivative wrt to azimuth angle
    # \param theta Polar angle \f$ \theta \f$
    # \param phi Azimuth angle \f$ \phi \f$
    # \returns \f$ \frac{\partial Y_l^m(\theta, \phi)}{\partial \phi}
    # = imY_l^m(\theta, \phi) \f$
    def dphi(self, theta, phi):
        """dphi(theta, phi)
        Derivative wrt to phi"""
        return (1j) * self.m * self(theta, phi)

    ## Gradient
    # \param r Radial coordiante \f$ r \f$
    # \param theta Polar angle \f$ \theta \f$
    # \param phi Azimuth angle \f$ \phi \f$
    # \returns
    # \f$ {\bf \vec{\nabla}}Y_l^m(\theta, \phi) =
    # \frac{\partial Y_l^m(\theta, \phi)}{\partial \theta}
    # {\bf{\hat{\theta}}} +
    # \frac{\partial Y_l^m(\theta, \phi)}{\partial \phi}
    # {\bf{\hat{\phi}}}  \f$
    def grad(self, r, theta, phi):
        """Returns the gradient in spherical coordinates (Polar) at
        Polar(r, theta, phi)"""
        return Polar(r, theta, phi).grad(0*r,
                                         self.dtheta(theta, phi),
                                         self.dphi(theta, phi))

    ## Tangential Vector Spherical Harmonic
    # \param r Radial coordiante \f$ r \f$
    # \param theta Polar angle \f$ \theta \f$
    # \param phi Azimuth angle \f$ \phi \f$
    # \returns \f$ {\bf Y}^m_l(\theta, \phi) = Y^m_l(\theta, \phi)
    # {\bf \hat{r}} \f$
    def Ylm(self, r, theta, phi):
        """Y vector harmonic = r*Y"""
        return r*Polar(one(theta, phi), theta, phi) * self(theta, phi) / r

    ## Normal Vector Spherical Harmonic
    # \param r Radial coordiante \f$ r \f$
    # \param theta Polar angle \f$ \theta \f$
    # \param phi Azimuth angle \f$ \phi \f$
    # \returns
    # \f$ {\bf \Psi}^m_l = r {\bf \vec{\nabla}} Y^m_l(\theta, \phi)  \f$
    def Psilm(self, r, theta, phi):
        """Psi vector harmonic = r * grad Y"""
        return r * self.grad(r, theta, phi)

    ## Binormal Vector Spherical Harmonic
    # \param r Radial coordiante \f$ r \f$
    # \param theta Polar angle \f$ \theta \f$
    # \param phi Azimuth angle \f$ \phi \f$
    # \returns  \f$ {\bf \Phi}^m_l = {\bf \vec{r}} \times
    # {\bf \vec{\nabla}} Y_l^m(\theta, \phi) \f$
    def Philm(self, r, theta, phi):
        """Phi vector harmonic = r X grad Y"""
        return Polar(r, theta, phi) ^ self.grad(r, theta, phi)
