"""Speed"""
## \namespace geo.politics.units._speed
# <a href="http://en.wikipedia.org/wiki/Speed">Speed</a>
from ._bases import *
from .times import hour2second, hour
from .length import mile2ft, ft2m, nauticalmile2meter, mile

__all__ = ('mph', 'mph2mps', 'mph2kph', 'kph2mph', 'mps2mph', 'kts2mph',
           'mps2mach', 'SaffirSimpson2mph', 'Beaufort2mph', 'mph2Beaufort',
           'Fujita2mph', 'mph2Fujita', 'EF2mph', 'mph2EF', 'TORRO2mph',
           'mph2TORRO', 'Dvorak2mph', 'mph2Dvorak')

mph = mile / hour


## <a href="http://en.wikipedia.org/wiki/Metre_per_second">Meters per seconds
# </a>
mph2mps = (mile2ft*ft2m)/hour2second

## <a href="http://en.wikipedia.org/wiki/Kilometres_per_hour">kph</a>
mph2kph = (mph2mps/(~hour2second)/constants.kilo)

kph2mph = ~mph2kph

## <a href="http://en.wikipedia.org/wiki/Mph">mph</a>
mps2mph = ~mph2mps


## knots: <a href="http://en.wikipedia.org/wiki/Knot_(unit)">Knots</a>.
kts2mph = nauticalmile2meter.__mul__(mps2mph/hour2second,
                                     From="knots",
                                     To="miles_per_hour")


## Mach Speed: <a href="http://en.wikipedia.org/wiki/Mach_number">Mach</a>
# and the
# <a href=
# http://docs.scipy.org/doc/scipy/reference/constants.html?highlight=mach">
# speed of sound</a>.
mps2mach = SI(1./constants.speed_of_sound, 'meters_per_second', 'Mach')


## The <a href="http://en.wikipedia.org/wiki/Saffir-Simpson_Hurricane_Scale">
# Saffir-Simpson Scale</a>.
# \param c SS number
# \returns  \f$ v  = [\alpha^{-1}-100]\cdot 10^{c/15} \f$
# @image html saffir-simpson-sm.gif
# @image latex saffir-simpson-sm.gif
def SaffirSimpson2mph(c):
    return mps2mph(1/(constants.fine_structure) - 100)*10**(c/15.)


## A class for the general form of common wind scales
class _WindScale(object):
    """Class for windscale function."""

    ## Construct from parameters (duh?):
    # \param a scale
    # \param b offset
    # \param c exponent
    def __init__(self, a, b, c):
        ## scale
        self.a = a
        ## offset
        self.b = b
        ## exponent
        self.c = c

    ## Scale to meter per seconds.
    # \param n Scale number
    # \returns  \f$ v(n) = a(n+b)^c \f$
    def __call__(self, n):
        return self.a*(n+self.b)**self.c

    ## Inverse Function:
    # \returns _IWS \f$ f(v) = \sqrt[c]{v/a}-b \f$
    def __invert__(self):
        return _IWS(self.a, self.b, self.c)


## Inverse _WindScale
class _IWS(_WindScale):
    """Inverse Windscale."""

    def __invert__(self):
        return _WindScale(self.a, self.b, self.c)

    ## Inverse Formula:
    # \param v A speed
    # \returns \f$ \sqrt{c}{v/a} - b \f$
    def __call__(self, v):
        return (v/self.a)**(1./self.c)-self.b


##  The <a href="http://en.wikipedia.org/wiki/Beaufort_scale">
# Beaufort Scale</a> to mph:\n \f$
# V = 1.870B^{\frac{3}{2}} \f$
Beaufort2mph = _WindScale(1.87, 0., 1.5)

## mph to Beaufort
mph2Beaufort = ~Beaufort2mph

##  The <a href="http://en.wikipedia.org/wiki/Fujitat_scale">Fujitat
# Scale</a> to mph:\n
# \f$ V =  14.1(F+2)^{\frac{3}{2}} \f$
Fujita2mph = _WindScale(14.1, 2, 1.5)

## mph to Fujita
mph2Fujita = ~Fujita2mph

##  The <a href="http://en.wikipedia.org/wiki/Enhanced_Fujitat_scale">
# Enhanced Fujitat Scale</a> to mph:\n
# \f$ V = 1.32623732(F+6.80405916)^{2.03088117} \f$
EF2mph = _WindScale(1.32623732, 6.80405916, 2.03088117)

## mph to Enhanced Fujita
mph2EF = ~EF2mph

##  The <a href="http://en.wikipedia.org/wiki/TORRO_scale">TORRO
# Scale</a> to mph:\n
# \f$ V = 5.289(T+4)^{\frac{3}{2}} \f$
TORRO2mph = _WindScale(5.289, 4, 1.5)

## mph to Torro
mph2TORRO = ~TORRO2mph

##  The <a href="http://en.wikipedia.org/wiki/Dvorak_technique">
# Dvorak T-Number</a> to mph:\n
# \f$ V = 0.8568(T+6.56137)^{1.82951101} \f$
Dvorak2mph = _WindScale(0.8568568, 6.56137899, 1.82951101)

## mph to Dvorak
mph2Dvorak = ~Dvorak2mph
