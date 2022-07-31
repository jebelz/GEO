"""Time"""
## \namespace geo.politics.units._times
# <a href="http://en.wikipedia.org/wiki/Time">Time</a>
from ._bases import *

__all__ = ('s2Hz', 's2us', 'us2s', 'us', 'hour2minute', 'minute2second', 'hour2second', 'hour', 'minute', 'shake', 'tropical_year', 'siderial_year', 'SOLAR_DAY', 'SIDERAL_DAY', 'EARTH_ROTATION_PERIOD', 'day', 'year', 'year_leap', 'year_julian', 'year_gregorian', 'year2dogyear', 'planck_t')


class Time(Unit):
    """Convert from/to seconds:

    seconds = 1 * unit  (convert 1 unit to seconds)
    units = 1 / unit   (convert 1 second to units)
    units = unit(1)    (convert 1 second to units)
    """


## <a href="http://en.wikipedia.org/wiki/Hertz">Hz</a>.
s2Hz = functools.partial(operator.div, 1)


## Mircoseconds.
s2us = SI(1/constants.micro, 'second', 'microsecond')
us2s = ~s2us

us = us2s.unit(Time)

## <a href="http://en.wikipedia.org/wiki/Hour
hour2minute = SI(constants.hour/constants.minute, 'hour', 'minute')

## <a href="http://en.wikipedia.org/wiki/Minute
minute2second = SI(constants.minute, 'minute', 'second')
hour2second = (hour2minute * minute2second)

hour = hour2second.unit(Time)
minute = minute2second.unit(Time)


## <a href="http://en.wikipedia.org/wiki/Shake_(unit)">Shake</a>
shake = Time(10*constants.nano)

## <a href="http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_2.html">
# Tropical Year</a>
tropical_year = scipy.poly1d([-0.53/86400., 365.+(5+(48+45.5/60.)/60.)/24.],
                       variable="[Julian Centuries-Y2K]")

siderial_year = scipy.poly1d([-0.01/86400., 365.+(6+(9+9.5/60.)/60.)/24.],
                             variable="[Julian Centuries-Y2K]")

SOLAR_DAY = Time(86400.003)

SIDERAL_DAY = Time(86164.094)

EARTH_ROTATION_PERIOD = Time(86164.102)

day = Time(24 * float(hour))
year = Time(365 * float(day))
year_leap = Time(1+float(year))
year_julian = Time(365.25 * float(day))
year_gregorian = Time(365.2425 * float(day))


## According to Cesar Milan.
year2dogyear = scipy.poly1d([4, 13], variable='years')

planck_t = Time(constants.hbar * constants.G / constants.c**5) ** 0.5
