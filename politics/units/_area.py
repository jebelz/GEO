"""Area"""
## \namespace geo.politics.units._area
# <a href="http://en.wikipedia.org/wiki/Area">area</a>
from ._bases import *
from .length import ft2yd, m2ft

__all__ = ('barn', 'hectare2sqmeter', 'm2ft', 'acre', 'acre2hectare', 'acre2sqft', 'ft2yd', 'section2acre', 'acre2sqmeter', 'twp2acre', 'acre2sqyd', 'ha', 'twp2section', 'hectare2acre')


class Area(Unit):
    """Convert from/to square meters:

    sqaure meters = 1 * unit  (convert 1 unit to meters**2)
    units = 1 / unit   (convert 1 meter**2 to units)
    units = unit(1)    (convert 1 meter**2 to units)
    """


## <a href="http://en.wikipedia.org/wiki/Acre">Acre</a>
acre2sqyd = SI(constants.acre/constants.yard**2, 'acre', 'yard**2')
acre2sqft = (acre2sqyd * ((~ft2yd)**2))

## <a href="http://en.wikipedia.org/wiki/Hectare">Hectare</a>
hectare2sqmeter = SI(constants.hecto**2, 'hectare', 'meters**2')

hectare2acre = (hectare2sqmeter * (m2ft**2) * (~acre2sqft))
acre2hectare = ~hectare2acre

acre2sqmeter = (acre2hectare * hectare2sqmeter)

section2acre = SI(640, 'section', 'acre')

## Township
twp2section = SI(36, 'twp', 'section')
twp2acre = (twp2section * section2acre)


acre = acre2sqmeter.unit(Area)

ha = hectare2sqmeter.unit(Area)

## <a href="http://en.wikipedia.org/wiki/Barn_(unit)">barn</a>.
barn = Area(1e-28)
