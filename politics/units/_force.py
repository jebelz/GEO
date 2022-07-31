"""Force"""
## \namespace geo.politics.units._force
#<a href="http://en.wikipedia.org/wiki/Force">Force</a>
from ._bases import *

__all__ = ('kgf2N', 'N2kgf', 'N2dyn', 'N2lbF', 'N2pdl', 'lbs2N')

kgf2N = SI(constants.g, "kgF", "N")
N2kgf = ~kgf2N
N2dyn = SI(constants.kilo/constants.centi, "N", "dyn")
N2lbF = SI(0.22481, "N", "lbF")
N2pdl = SI(7.2330, "N", "pdl")


## <a href="http://en.wikipedia.org/wiki/Pound-force">Pound Force</a>
lbs2N = SI(constants.pound_force, 'lbs', 'N')
