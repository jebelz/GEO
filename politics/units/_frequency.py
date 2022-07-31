"""Freq."""
## \namespace geo.politics.units._frequency Hz
from ._bases import *
from .times import minute2second
from .times import s2Hz as Hz2s

__all__ = ('Hz2MHz', 'MHz2Hz', 'Hz2GHz', 'GHz2Hz', 'hz2rpm', 'rpm2hz', 'hz2rps', 'RPM', 'A440', 'key2Hz', 'key_from_Hz')


class Frequency(Unit):
    """Convert from/to Hertzs:

    Hertzs = 1 * unit  (convert 1 unit to Hertzs)
    units = 1 / unit   (convert 1 Hertz to units)
    units = unit(1)    (convert 1 Hertz to units)
    """


## Megahertz.
Hz2MHz = SI(1/constants.mega, 'Hz', 'GHz')
MHz2Hz = ~Hz2MHz


## Gigahertz.
Hz2GHz = SI(1/constants.giga, 'Hz', 'GHz')
GHz2Hz = ~Hz2GHz


##<a href="http://en.wikipedia.org/wiki/Revolutions_per_minute">RPMs</a>
hz2rpm = SI(float(minute2second), 'Hz', 'RPM')

# <a href="http://en.wikipedia.org/wiki/Arc_length#Arcs_of_circles">RPM v Hz</a>
rpm2hz = ~hz2rpm

## <a href="http://en.wikipedia.org/wiki/Radians_per_second">Angular Rate</a>
hz2rps = SI(2*constants.pi, 'Hz', 'radians_per_second')

RPM = rpm2hz.unit(Frequency)


## <A href="http://en.wikipedia.org/wiki/A440_(pitch_standard)">Concert A</a>
A440 = Frequency(440)


## Piano key to frequency (H
def key2Hz(n, A=A440):
    """f = key(n, [concert_A=440])
    f (hz) for n = 1 to 88.
    """
    return pow(2, (n-49)/12.) * A


## Convert a frequency to a key number.
def key_from_Hz(f, A=A440):
    """n = key_from_Hz(f, concert_A=440)
    f (Hz)
    n is the piano key for f."""
    return 12 * scipy.log2(float(f)/A) + 49.
