"""Intensity"""
## \namespace geo.politics.units._intensity
# <a href="http://en.wikipedia.org/wiki/Intensity_(physics)">Intensity<a>
#
from ._bases import *
from ._logs import _Log

__all__ = ('apparent_magnitude', 'intensity2magnitude', 'magnitude2intensity', 'Langley2Jpersqm', 'solar_constant', 'sqm2watts_solar', 'L_solar', 'jansky')

## TBD, same as intensity2magnitude.
apparent_magnitude = _Log(base=10, prefix=fractions.Fraction(-4, 10))


## <a href="http://en.wikipedia.org/wiki/Magnitude_(astronomy)">Magnitude</a>
# from intensity difference.
# \param i Intensity ratio.
# \returns \f$ m_1-m_2 = -2.5 \log_{10}{\frac{I_1}{I_2}} \f$
def intensity2magnitude(i):
    """m1-m2 = intensity2magnitude(i1/i2)"""
    return -2.5*scipy.log10(i)


## <a href="http://en.wikipedia.org/wiki/Magnitude_(astronomy)">Magnitude</a>
# to intensity ratio.
## \f$ I_1/I_2 = 10^{(m_1-m_2)/(-2.5)} \f$
def magnitude2intensity(m):
    """I1/I2 = magnitude2intensity(m1-m2)"""
    return 10**(m/(-2.5))

## <a href="http://en.wikipedia.org/wiki/Langley_(unit)">Langley</a>
Langley2Jpersqm = SI(418400.00, 'langley', 'J per sqm')


## <a href="http://en.wikipedia.org/wiki/Solar_constant">Solar Constant</a>
solar_constant = 1412.
sqm2watts_solar = SI(solar_constant, "sq-meter (solar)", "W")

## <a href="http://en.wikipedia.org/wiki/Solar_luminosity">Solar Luminosity</a>
L_solar = 3.939e26


## <a href="http://en.wikipedia.org/wiki/Jansky">Jansky</a>.
jansky = SI(1.e-26, "Jy", "W per msq per Hz")
