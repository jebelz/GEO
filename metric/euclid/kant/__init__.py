"""Kant is the home Noncartesian representations of Vectors.

Contents:

========================================================================
Class                  |  Constructor       |        Description       |
                       | (from Cartesian)   |                          |
-----------------------+--------------------+--------------------------|
Polar(r, theta, phi)   |polar(x, y, z)      | Spherical[1] Coordinates |
.......................|....................|..........................|
Cylinder(rho, theta, z)|cylinder(x, y, z)   | Cylindrical Coordinates  |
.......................|....................|..........................|
Parabola(u, v, theta)  |parabola(x, y, z)   | Parabolic Coordinates    |
========================================================================

And the CoordinateFactory metaclass to brew you own (at runtime).

These object can play amongst themselves, and with Cartesian Vectors,
e.g.:

>>>Vector(3, 4, 12) - Polar(13, 0.394791197, 0.927295218002)

is NULL.


[1] The term Polar is used, b/c a Spherical vector is something entirely
different (see wigner/eckart.py)"""

## \namespace geo.metric.euclid.kant Non-cartesian wrappers.
# @image html t_polar.png
# @image latex t_polar.png
import functools

from . spherical import *
from . cylindrical import *
from . parabolic import *


## Non-Cartesian Coordinate Factory (It is this easy).
class CoordinateFactory(type):
    """C = CoordinateFactory(fromcart,
    tocart,
    components=('a', 'b', 'c'),
    name='Noncartesian')

    where the conversion functions are given:
    a, b, c = fromcart(x, y, z)
    x, y, z = tocart(a, b, c)
    """

    ## Class Factory Method
    # \param fromcart \f$(a, b, c) = F(x, y, z)\f$
    # \param tocart   \f$(x, y, z) = F^{-1}(a, b, c)\f$
    # \param components=('a', 'b', 'c') Names of the new coordinates
    # \param name='Noncartesian'  Name of the new coordinate system
    # \returns type Class derived from _bases.NCE
    def __new__(mcs, fromcart, tocart,
                components=('a', 'b', 'c'), name='Noncartesian'):
        from . import _bases
        return type(name,
                    (_bases.NonCartesianBase,),
                    {'_fromcart': staticmethod(fromcart),
                     '_tocart': staticmethod(tocart),
                     'components': components})


## Chain vector to cls via fromvector.
def _helper(cls):
    """Help make vectors"""

    def decorator(func):
        """Decoator for cls"""

        @functools.wraps(func)
        def helper(x, y, z):
            """convert x, y, z into instance of cls"""
            from .... import Vector
            return cls.fromvector(Vector(x, y, z))
        return helper
    return decorator


@_helper(Polar)
def polar():
    """Polar = polar(x, y, z)"""


@_helper(Cylinder)
def cylinder():
    """Cylinder = cylinder(x, y, z)"""


@_helper(Parabola)
def parabola():
    """Parabola = parabola(x, y, z)"""
