"""Length"""
## \namespace geo.politics.units._length
# <a href="http://en.wikipedia.org/wiki/Length">length</a>.
from ._bases import *
from .times import year_julian

__all__ = ('SEMIMAJOR_AXIS', 'meridonal_radius', 'ft2m', 'ft', 'm2ft', 'hand2m', 'pica2point', 'inch2pica', 'ft2inch', 'inch2ft', 'ft2yd', 'chain2ft', 'furlong2chain', 'mile2furlong', 'lea2mile', 'yd2ft', 'mile2ft', 'ft2mile', 'mile2km', 'mile', 'rod2yard', 'chain2link', 'nauticalmile2meter', 'geographicmile2meter', 'AU', 'light_year', 'parsec', 'planck_l', 'fm', 'angstrom')


class Length(Unit):
    """Convert from/to meters:

    meters = 1 * unit  (convert 1 unit to meters)
    units = 1 / unit   (convert 1 meter to units)
    units = unit(1)    (convert 1 meter to units)
    """


## Earth's semi-major axis sho
SEMIMAJOR_AXIS = 6378137.0

## Earth's mean meridonal radius.
meridonal_radius = 6367449.1457001315


## To <a href="http://en.wikipedia.org/wiki/Foot_(length)">The Foot</a>
ft2m = SI(constants.foot, 'ft', 'meters')

ft = ft2m.unit(Length)

## From <a href="http://en.wikipedia.org/wiki/Foot_(length)">The Foot</a>
m2ft = ~ft2m

## Hands.
hand2m = SI(10.16e-2, "hand", "m")


## <a href="http://en.wikipedia.org/wiki/Pica_(typography)">Pica</a> tp
# <a href="http://en.wikipedia.org/wiki/Point_(typography)">Point</a>.
pica2point = SI(12., 'pica', 'point')
inch2pica = SI(6., 'inch', 'pica')

## <a href="http://en.wikipedia.org/wiki/Inches">Inches</a> to Feet.
ft2inch = SI(constants.foot/constants.inch, 'ft', 'inch')
inch2ft = ~ft2inch

## Feet to <a href="http://en.wikipedia.org/wiki/Yard">Yards</a>
ft2yd = SI(constants.foot/constants.yard, 'ft', 'yard')

## <a href="http://en.wikipedia.org/wiki/Chain_(length)">Chains</a> to Feet.
chain2ft = SI(66., 'chain', 'ft')

## <a href="http://en.wikipedia.org/wiki/Furlong">Furlongs</a> to Chains
furlong2chain = SI(constants.deka, 'furlong', 'chain')
## <a href="http://en.wikipedia.org/wiki/Mile
mile2furlong = SI(8., 'mile', 'furlong')
lea2mile = SI(3., 'lea', 'mile')


yd2ft = ~ft2yd


mile2ft = (mile2furlong * furlong2chain * chain2ft)
ft2mile = ~mile2ft

# <a href="http://en.wikipedia.org/wiki/Kilometer">Klics</a> and miles.
mile2km = ((mile2ft * ft2m).__div__(constants.kilo, 'km'))

mile = (mile2km/(1/constants.kilo)).unit(Length)

## <a href="http://en.wikipedia.org/wiki/Rod_(length)">Rod</a>.
rod2yard = SI(5.5, 'rod', 'yard')

## <a href="http://en.wikipedia.org/wiki/Link_(unit)">Link</a>.
chain2link = SI(constants.hecto, 'chain', 'link')

## <a href="http://en.wikipedia.org/wiki/Link_(unit)">Link</a>.
chain2link = SI(constants.hecto, 'chain', 'link')

## <a href="http://en.wikipedia.org/wiki/Nautical_mile">Nautical Mile</a>
nauticalmile2meter = SI(constants.nautical_mile, 'nautical mile', 'meter')

## <a href="http://en.wikipedia.org/wiki/Geographical_mile">Geographoc Mile</a>
geographicmile2meter = SI(SEMIMAJOR_AXIS*2*constants.pi/360./60.,
                          From='geographic mile',
                          To='meter')


## <a href="http://en.wikipedia.org/wiki/Astronomical_unit">Astronomical
# Unit</a>
AU = Length(149597870700)

## <a href="http://en.wikipedia.org/wiki/Light-year">Light Year</a>
light_year = Length(constants.speed_of_light * float(year_julian))

## <a href="http://en.wikipedia.org/wiki/Parsec">Parsec</a>, with 1st order
# Taylor expansion of arc-tangent function.
parsec = Length(1/sum([f(scipy.radians(1./3600.) / float(AU)) for f in map(scipy.poly1d, ([1, 0], [-1./3, 0, 0, 0]))]))


planck_l = Length(constants.hbar * constants.G / constants.c**3) ** 0.5


## <a href="http://en.wikipedia.org/wiki/Femtometre">Fermi<a/>
fm = Length(1e-15)

angstrom = Length(1e-10) # Fail
