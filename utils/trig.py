"""A place to store trig functions using degrees-- so If you don't have numpy
you can use math-- but just have numpy."""
#pylint: disable=W0611, C0111, C0301

## \namespace geo.utils.trig Trig functions in degrees from
#  <a href="http://www.numpy.org/">numpy</a> or
#  <a href="http://docs.python.org/2/library/math.html">math</a>
import operator
import functools

# Here we completely violate ths style guide so that this module
# olds the funtion, regardless of its source
try:
    from numpy import (degrees, radians, sin, cos, tan, arccos, arctan2,
                       arctan, pi, arcsin, log, sqrt, sign, exp)
except ImportError as err:
    from math import degrees, radians, sin, cos, tan, pi, sqrt, log, exp
    # importh python with numpy name convensions
    from math import acos as arccos
    from math import asin as arcsin
    from math import atan as arctan
    from math import atan2 as arctan2
    from math import copysign as _copysign
    sign = functools.partial(_copysign, 1)
    print "numpy import failed: using math built-in"
    NUMPY = False
else:
    NUMPY = True

## \f$2\pi\f$
TWO_PI = 2*pi

## \f$ \frac{2\pi}{360} \f$
DEGRAD = (TWO_PI/360)

## A collection of angle measuring names, and their range
CIRCUMFERENCE = {'degree': 360,
                 'radian': TWO_PI,
                 'arcsecond': 1296000,
                 'turn': 1.,
                 'quadrant': 4,
                 'sextant': 6,
                 'hour': 24,
                 'point': 32,
                 'grad': 400,
                 'arcminute': 21600}


## decorator converts input radians to degrees
def from_degrees(func):
    """func(radians) --> func(degrees) """
    @functools.wraps(func)
    def funcd(x):
        return func(radians(x))
    return funcd


## decorator converts input radians to degrees
def to_degrees(func):
    """radian = f(x) --> degrees = f(x)"""
    @functools.wraps(func)
    def dfunc(*args):
        return degrees(func(*args))
    return dfunc


## \f$ \cos{\frac{2\pi}{360}\theta^{\circ}} \f$ \n
#  cosine (degrees) using from_degrees() decorator
@from_degrees
def cosd(x):
    """cosine of argument in degrees"""
    return cos(x)


## \f$ \sin{\frac{2\pi}{360}\theta^{\circ}} \f$ \n
#  sine (degrees)
@from_degrees
def sind(x):
    """sine of argument in degrees"""
    return sin(x)


## \f$ \tan{\frac{2\pi}{360}\theta^{\circ}} \f$ \n
#  tangent (degrees)
@from_degrees
def tand(x):
    """tangent of argument in degrees"""
    return tan(x)


## \f$ f(y, x) = \frac{360}{2\pi}\arctan_2(y, x) \f$ \n
#  2-argument inverse tangent, in degrees via to_degrees() decorator
@to_degrees
def arctand2(y, x):
    """arctangent in degrees (2-argument form)"""
    return arctan2(y, x)


## \f$ f(x) = \frac{360}{2\pi}\arctan(x) \f$ \n
#  1-argument inverse tangent, in degrees via to_degrees() decorator
@to_degrees
def arctand(x):
    """arctangent in degrees (1-argument form)"""
    return arctan(x)


## \f$ f(x) = \frac{360}{2\pi}\arccos(x) \f$ \n
#  1-argument inverse tangent, in degrees via to_degrees() decorator
@to_degrees
def arccosd(x):
    """arccose in degrees"""
    return arccos(x)


## \f$ f(x) = \frac{360}{2\pi}\arcsin(x) \f$ \n
#  1-argument inverse tangent, in degrees via to_degrees() decorator
@to_degrees
def arcsind(x):
    """arcsine in degrees"""
    return arcsin(x)


## \f$ \Re{z} \f$
Re = operator.attrgetter('real')


## \f$ \Im{z} \f$
Im = operator.attrgetter('imag')


## Complex Conjugate by any means necessary.
# \param z
# \returns \f$ z^* \equiv \Re{z} - i \Im{z} \f$
def conj(z):
    """complex conjugate of z-- by any means possible."""
    try:
        return z.conjugate()
    except AttributeError:
        try:
            return z.conj()
        except AttributeError:
            try:
                return z.real - z.imag * (1j)
            except AttributeError:
                x = (1 + 1j) * z
                y = (1 - 1j) * z
                return (
                    (x + y.conjugate()) / (1 + 1j)  -
                    (x - y.conjugate()) / (1 - 1j)
                ) / 2
