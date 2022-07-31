"""The Euclidean subpackage:
Core:

scalar.py   Rank 0
vector.py   Rank 1
tensor.py   Rank 2
three.py    Rank 3 and 4

affine.py   General Transformations of Euclidean Space."""
## \namespace geo::metric::euclid
# <a href="http://en.wikipedia.org/wiki/Euclidean_space">\f$E^3\f$</a> Animals.
# @image html euclid.jpg
# @image latex euclid.jpg
from .scalar import *
from .vector import *
from .tensor import *
from .affine import *
from .kant import *
from .three import *
# this must be loaded now, or never- thanks to ranked
from .four import *


## <a href="http://en.wikipedia.org/wiki/Parity_(physics)">Parity</a> operator.
# \param T A tensor
# \returns \f$ P{\bf T} \f$
def parity(T):
    """P(T) = -DELTA_i'i...-DELTA_m'm * T_i...k"""
    return T * pow(-1, T.rank)  # cart before horse implementation
