"""UTM Functions via the Kruger Flattening Series."""

#pylint: disable=R0913,R0914

## \namespace geo.detic.mercator.kruger
# <a href="http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#Simplified_formulas">
#  Mecator flatenning series</a>

import numpy as np


## \f$ \frac{f}{2-f} \f$
def n(f):
    """n(f)"""
    return f / (2 - f)

## \f$ 1 + \frac{n^2}{4} + {n^4}{64} \f$
A = np.poly1d([1/64, 0., 1/4, 0., 1], variable='n')

## \f$ f(j) \equiv 1 \f$
UNITY = np.poly1d([1], variable='j')

## \f$ f(j) \equiv 2j \f$
TWO_J = np.poly1d([2, 0], variable='j')


## The Kruger Series polynimials in n, and their realizations.
class Kruger(object):
    """Kruger is composed of polynomials."""

    ## \f$ \alpha_1 = \frac{1}{2}n-\frac{2}{3}n^2+\frac{5}{16}n^3 \f$
    ALPHA_1 = np.poly1d([5/16, -2/3, 1/2, 0], variable='n')
    ## \f$ \alpha_2 = \frac{132}{48}n^2-\frac{3}{5}n^3 \f$
    ALPHA_2 = np.poly1d([-3/5, 13/48, 0, 0], variable='n')
    ## \f$ \alpha_3 = \frac{16}{240} n^3 \f$
    ALPHA_3 = np.poly1d([61/240, 0, 0, 0], variable='n')
    ## \f$ \beta_1 = \frac{1}{2}n-\frac{2}{3}n^2+\frac{5}{16}n^3 \f$
    BETA_1 = np.poly1d([5/16, -2/3, 1/2, 0], variable='n')
    ## \f$ \beta_2 = \frac{132}{48}n^2-\frac{3}{5}n^3 \f$
    BETA_2 = np.poly1d([-3/5, 13/48, 0, 0], variable='n')
    ## \f$ \beta_3 = \frac{16}{240} n^3 \f$
    BETA_3 = np.poly1d([61/240, 0, 0, 0], variable='n')
    ## \f$ \delta_1 = \frac{1}{2}n-\frac{2}{3}n^2+\frac{5}{16}n^3 \f$
    DELTA_1 = np.poly1d([5/16, -2/3, 1/2, 0], variable='n')
    ## \f$ \delta_2 = \frac{132}{48}n^2-\frac{3}{5}n^3 \f$
    DELTA_2 = np.poly1d([-3/5, 13/48, 0, 0], variable='n')
    ## \f$ \delta_3 = \frac{16}{240} n^3 \f$
    DELTA_3 = np.poly1d([61/240, 0, 0, 0], variable='n')

    ## \f$ \alpha = (\alpha_1, \alpha_2, \alpha_3) \f$
    ALPHA = (ALPHA_1, ALPHA_2, ALPHA_3)
    ## \f$ \beta = (\beta_1, \beta_2, \beta_3) \f$
    BETA = (BETA_1, BETA_2, BETA_3)
    ## \f$ \delta = (\delta_1, \delta_2, \delta_3) \f$
    DELTA = (DELTA_1, DELTA_2, DELTA_3)

    ## Helper sends n to polynomials and on to __init__.
    @classmethod
    def from_n(cls, n_, letter):
        """evaluate functions from 'n'"""
        return cls(*[func(n_) for func in getattr(cls, letter.upper())])

    ## \param *polynomials in 'n'.
    def __init__(self, *polynomials):
        self.polynomials = polynomials

    ## This is a complicated generic series calculator:
    # \param func1 \f$ f_1 \f$ A function
    # \param arg1 \f$ x_1 \f$ A geodetic coordinate
    # \param func2 \f$ f_2 \f$ Another funciton
    # \param arg2 \f$ x_2 \f$ The other geodetic coordinate
    # \param jpoly \f$P(j)\f$ A polynomial in the series sum.
    # \return \f$ \sum_{j=1}^3{ P(j)f_1(2jx_1)f_2(2jx_2)} \f$
    def __call__(self, func1, arg1, func2, arg2, jpoly):
        result = 0
        for count, poly in enumerate(self.polynomials):
            j = count + 1
            result += jpoly(j) * poly * func1(2*j*arg1) * func2(2*j*arg2)
            continue
        return result

## \f$ E_0 = 500,000 \f$
E0 = 500e3
## \f$ k_0 = 0.996 \f$
k0 = 0.9996
## \f$ N_0 = 0 \f$ in the Northren Hemisphere
N0_NORTH = 0.
## \f$ N_0 = 1,000,000 \f$ in the Southern Hemisphere
N0_SOUTH = 10000e3
## The N0 hash table
N0 = {-1: N0_SOUTH, 1: N0_NORTH}

## \f$ \lambda_0(Z) = 6Z-183 \f$
lambda0 = np.poly1d([6, -183], variable='zone')


## Total antipatten class-- this will be refactored
class _Base(object):
    """Base class for coordinate conversion"""

    ## Class helper:
    # \param geo.detic.ellipsoid.Ellipsoid
    @classmethod
    def fromellipsoid(cls, ellipsoid_):
        """Contruct form an ellipsoid."""
        return cls(ellipsoid_.a, ellipsoid_.finv)

    @classmethod
    def alpha(cls, *args, **kwargs):
        """Overidden by contructor."""
        return cls, args, kwargs

    beta = alpha
    delta = alpha

    ## Instaniate from Ellipsoid Parameters
    def __init__(self, a, finv):
        ## Semi-major axis
        self.a = a
        ## Inverse flatening
        self.finv = finv
        f = 1/finv
        ## Flatenning
        self.f = f
        ## ::n
        self.n = n(f)
        ## ::A
        self.A = A(self.n) * self.a / (1 + self.n)
        for letter in ('ALPHA', 'BETA', 'DELTA'):
            setattr(self, letter.lower(), Kruger.from_n(n, letter))

    ## \f$ [_{\phi, \lambda}f_{E, N}]^{-1} = _{E, N}f_{\phi, \lambda} \f$ \n
    # via ::CANNON.
    def __invert__(self):
        return CANNON[type(self)](self.a, self.finv)


## Function Wrapper converts LL to UTM
class ToUTM(_Base):
    """f = ToUTM(a, finv)"""

    ## \f$ t(\phi) = \sinh{(\tanh^{-1}{\sin{\phi}} - \frac{\sqrt{2n}}{1+n}
    # \tanh^{-1}{(\frac{\sqrt{2n}}{1+n} \sin{\phi})})} \f$
    def t(self, phi):
        """t(phi)"""
        c = 2 * np.sqrt(self.n) / (1 + self.n)
        sinphi = np.sin(phi)
        return np.sinh(np.arctanh(sinphi) -
                       c * np.arctanh(c * sinphi))

    ## \f$ \xi'(t, \lambda, \lambda_0) =
    # \tan^{-1}{\frac{t}{\cos{(\lambda-\lambda_0})}} \f$
    @staticmethod
    def csi_prime(t, lam, lam0):
        """xi'"""
        return np.arctan2(t, np.cos(lam-lam0))

    ## \f$ \eta'(t, \lambda, \lambda_0) =
    # \tan^{-1}{\frac{\sin{(\lambda-\lambda_0)}}{\sqrt{1+t^2}}} \f$
    @staticmethod
    def eta_prime(t, lam, lam0):
        """eta'"""
        return np.arctanh(np.sin(lam-lam0) / np.sqrt(1 + t**2))

    ## \f$ \sigma(\xi', \eta') = 1 +
    # {\bf \bar{\alpha}}(\cos, \xi', \cosh, \eta', f(j)=2j) \f$
    def sigma(self, csi_prime, eta_prime):
        """sigma"""
        return 1 + self.alpha(np.cos, csi_prime,
                              np.cosh, eta_prime,
                              jpoly=TWO_J)

    ## \f$ \tau(\xi', \eta') =
    # {\bf \bar{\alpha}}(\sin, \xi', \sinh, \eta', f(j)=2j) \f$
    def tau(self, csi_prime, eta_prime):
        """tau"""
        return self.alpha(np.sin, csi_prime,
                          np.sinh, eta_prime,
                          jpoly=TWO_J)

    ## \f$ E(\xi', \eta') = E_0 + k_0 A [\eta' + {\bf \bar{\alpha}}
    # (\cos, \xi', \sinh, \eta', f(j)=2j)] \f$
    def E(self, csi_prime, eta_prime):
        """E"""
        return E0 + k0 * self.A * (eta_prime + self.alpha(np.cos, csi_prime,
                                                          np.sinh, eta_prime,
                                                          jpoly=TWO_J))

    ## \f$ N(\xi', \eta') = N_0 + k_0 A [\xi' + {\bf \bar{\alpha}}
    # (\sin, \xi', \cosh, \eta', f(j)=2j)] \f$
    def N(self, csi_prime, eta_prime):
        """N"""
        return N0 + k0 * self.A * (csi_prime + self.alpha(np.sin, csi_prime,
                                                          np.cosh, eta_prime,
                                                          jpoly=TWO_J))

    ## \f$ k(\phi, \lambda, \lambda_0, \sigma, \tau, t) =
    # \frac{k_0 A}{a}
    # \sqrt{ (1 + \left(\frac{(1-n)\tan{\phi}}{1+n}\right)^2
    # \frac{ \sigma^2 + \tau^2}{t^2 + \cos{\lambda-\lambda_0}}} \f$
    def k(self, phi, lam, lam0, sig, tau, t):
        """k"""
        return (k0 * self.A / self.a *
                np.sqrt((1 + ((1-self.n)*np.tan(phi)/(1+self.n))**2) *
                        (sig**2 + tau**2) / (t**2 + np.cos(lam-lam0))))

    ## \f$ \gamma(\lambda, \lambda_0, \sigma, \tau, t) =
    # \tan^{-1}{\frac
    #{\tau\sqrt{1+t^2} + \sigma t \tan{\lambda-\lambda_0}}
    #{\sigma\sqrt{1+t^2} - \tau t \tan{\lambda-\lambda_0}}
    # } \f$
    @staticmethod
    def gamma(lam, lam0, sig, tau, t):
        """gamma"""
        return np.arctan2(
            tau * np.sqrt(1 + t**2) + sig * t * np.tan(lam-lam0),
            sig * np.sqrt(1 + t**2) - tau * t * np.tan(lam-lam0)
            )

    ## Convert to UTM
    # \param \phi Latitude
    # \param \lambda Latitude
    def __call__(self, phi, lam, lam0):
        t = self.t(phi)
        csi_prime = self.csi_prime(t, lam, lam0)
        eta_prime = self.eta_prime(t, lam, lam0)
        return self.E(csi_prime, eta_prime), self.N(csi_prime, eta_prime)


## Function Wrapper converts UTM to LL
class FromUTM(_Base):
    """f = ToUTM(a, finv)"""

    ## \f$ \xi(N) = \frac{N-N_0}{k_0 A} \f$
    def csi(self, N):
        """xi"""
        return (N-N0) / (k0 * self.A)

    ## \f$ \eta(E) = \frac{E-E_0}{k_0 A} \f$
    def eta(self, E):
        """eta"""
        return (E-E0) / (k0 * self.A)

    ## \f$ \xi'(\xi, \eta) = \xi - {\bf \bar{\beta}}
    # (\sin, \xi, \cosh, \eta, f(j)=1) \f$
    def csi_prime(self, csi, eta):
        """xi'"""
        return csi - self.beta(np.sin, csi,
                               np.cosh, eta,
                               jpoly=UNITY)

    ## \f$ \eta'(\xi, \eta) = \eta - {\bf \bar{\beta}}
    # (\cos, \xi, \sinh, \eta, f(j)=1) \f$
    def eta_prime(self, csi, eta):
        """eta'"""
        return eta - self.beta(np.cos, csi,
                               np.sinh, eta,
                               jpoly=UNITY)

    ## \f$ \sigma'(\xi, \eta) = 1 - {\bf \bar{\beta}}
    # (\cos, \xi, \cosh, \eta, f(j)=2j) \f$
    def sigma_prime(self, csi, eta):
        """sigma'"""
        return 1 - self.beta(np.cos, csi,
                             np.cosh, eta,
                             jpoly=TWO_J)

    ## \f$ \tau'(\xi, \eta) = 1 - {\bf \bar{\beta}}
    # (\sin, \xi, \sinh, \eta, f(j)=2j) \f$
    def tau_prime(self, csi, eta):
        """tau'"""
        return 1 - self.beta(np.sin, csi,
                             np.sinh, eta,
                             jpoly=TWO_J)

    ## \f$ \chi(\xi', \eta') = \sin^{-1}{\frac{\sin{\xi'}}{\cosh{\eta'}}} \f$
    @staticmethod
    def chi(csi_prime, eta_prime):
        """chi"""
        return np.arcsin(np.sin(csi_prime) / np.cosh(eta_prime))

    ## Conversion Function
    # \param E Easting
    # \param N Northing
    # \param zone Zone
    # \param hemi Hemisphere's primitive (or not)
    # \returns (\f$ \phi=\chi + {\bf \bar{\delta}}(\sin, \chi , 1),
    # \lambda = \lambda_0(Z) + \tan^{-1}{\frac{\sinh{\eta'}}{\cos{\xi'}}} )\f$
    def __call__(self, E, N, zone, hemi):
        csi = self.csi(N)
        eta = self.eta(E)
        csi_prime = self.csi_prime(csi, eta)
        eta_prime = self.eta_prime(csi, eta)
        chi = self.chi(csi_prime, eta_prime)
        phi = chi + self.delta(np.sin, chi, lambda dum: 1, 1)

        lam = lambda0(zone) + np.arctan2(np.sinh(eta_prime), np.cos(csi_prime))

        return phi, lam


## The Cannonical Pairing (for inversion)
CANNON = {ToUTM: FromUTM, FromUTM: ToUTM}
