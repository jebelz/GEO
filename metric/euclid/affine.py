"""Affine Transformation: A general linear function from R3 --> R3:

A(x) = M*x + b

Where M is either a Tensor of any kind, or a rotation object.
b is a Vector.

Helmert:
With the helmert() function, you get a standard geodesy affine transformation
that includes a scaling factor."""
## \namespace geo::metric::euclid::affine
# <a href="http://en.wikipedia.org/wiki/Affine_transformation">
# Affine transformations</a>
import fractions
from . import euclid

DIMENSIONS = 3

__all__ = ('Affine', 'helmert', 'wgs84_to_mgi')


## Limited <a href="http://en.wikipedia.org/wiki/Affine_transformation">Affine
# Transformations</a>.
class Affine(euclid.LinearMap):
    """Affine(linear_map, translation)

    linear_map: preferably a euclid.chart on SO3
    translation: euclid.Vector in E3

    Methods:
    == == == ==

    A(v)    applies transformation
    A*A'    composition: (A*A')(x) = A(A'(x))
    ~A      returns the inverse transformation

    A ** (n/m) computes any rational power of the thing.
    """
    ## Init: a callable linear_map and translation
    # \param linear_map Anything that rotates
    # \param translation Anything that translates
    def __init__(self, linear_map, translation):
        """see class docstring for signature"""
        ## Alias linear_map (or not, you could make it a shear, dilation, ...)
        self.linear_map = linear_map
        ## Translation
        self.translation = translation

    ## The vector is the translation
    @property
    def vector(self):
        """get the translation vector"""
        return self.translation

    ## Rotation portion
    @property
    def rotation(self):
        """Rotation portion of transform."""
        return self.linear_map

    ## Get Tensor version of the linear_map
    # \returns tensor A tensor.Tensor, regardless.
    @property
    def tensor(self):
        """get the linear_map as a tensor"""
        return self.linear_map.tensor

    def __str__(self):
        return "A(v) = %s(v) + %s" % (str(self.linear_map),
                                      str(self.translation))

    ## Affine transform:\n
    # \f$ A(\vec{v}) \rightarrow \vec{v}' = R(\vec{v}) + \vec{T} \f$
    # \param vector_ A vector
    # \returns Affine transformed vector
    def alibi_transform(self, vector_):
        """vector = affine(vector) is
        affine.translation + affine.linear_map(vector)"""
        return self.translation + (self.linear_map)(vector_)


    ## Convolution is left to right (for alibi)\n
    #  \f$ AA' =(R, T)(R', T') \rightarrow (R'(R), R'(T) + T) \f$ \n
    def __mul__(self, other):
        return type(self)(self.linear_map * other.linear_map,
                          self.linear_map(other.translation) +
                          self.translation)


    ## See root().
    # \return \f$ A^{\frac{1}{2}} \f$
    def sqrt(self):
        """Square root of transform."""
        return self.root(2)

    ## 1/n-th fractional transform.
    # \param n Root
    # \return \f$ A^{\frac{1}{n}} \f$
    def root(self, n):
        """n-th root of transform."""
        from operator import add
        M = (self.linear_map.versor() ** (1./n)).tensor()  # move to tensor
        T = self.translation * reduce(add, [M**p for p in range(n)], 0*M).I
        return type(self)(M, T)

    ## Rational powers are fractional transforms.
    # \returns \f$ A^{\frac{n}{m}} =
    # A^{\frac{1}{m}}\cdot A^{\frac{1}{m}}\cdots A^{\frac{1}{m}} \f$ n times
    def __pow__(self, n):
        """Usage:

        A ** 0         Affine(DELTA, NULL)
        A ** n         A*A*A...n times *A    n = 1, 2, 3
        A ** (-n)      (~A)**n
        A ** (1/m)     A.root(m) for fractions.Fraction(1, m) (see stdlib)
        A ** (n/m)     (A ** (1/m)) ** n
        A ** flt       A ** (n/m) with n, m from Fraction.from_float(flt)
                                                 witn m <= 10000
        """
        from operator import mul
        from itertools import repeat
        if isinstance(n, fractions.Fraction):  # guard on fractions
            result = self.root(n.denominator) ** n.numerator
        elif isinstance(n, float):
            result = self.__pow__(
                fractions.Fraction.from_float(n).limit_denominator(10000)
            )
        else:
            if n < 0:  # power branch
                result = (~self) ** (-n)
            elif n == 0:
                result = type(self)(self.linear_map ** 0, self.translation * 0)
            elif n == 1:
                result = self
            else:
                i, f = divmod(n, 1)
                result = reduce(mul, repeat(self, int(i)), self**f)
        return result

    ## Inversion
    #  \returns
    #  \f$ AA^{-1} = ({\bf 1}, 0) \rightarrow A^{-1} =
    # (R^{-1}, -R^{-1}(T)) \f$
    #  \throws err.ProjectionInversionError Non-invertible projection
    #  \throws err.RejectionInversionError  Non-invertible rejection
    #  \throws ZeroDivisionError If it could not be ID's as priors.
    def __invert__(self):
        """inverse transformation-- or raise an error."""
        try:
            # If the matrix part is not a linear_map, this could fail
            inv_rot = ~(self.linear_map)
        except ZeroDivisionError as err:
            from ...utils import exceptions as err
            # get trace to try and ID the problem
            trace = round(self.linear_map.trace())
            # A lil' dictionary of known traces...
            AFFINE_ERROR = {1: err.ProjectionInversionError,
                            DIMENSIONS-1: err.RejectionInversionError}
            # get correct error, or None
            affine_error = AFFINE_ERROR.get(trace)
            # raise special error or regular error
            raise affine_error or err
        return type(self)(inv_rot, -(inv_rot(self.translation)))

    ## Inverse
    @property
    def T(self):
        """Inverse..."""
        raise TypeError("Affine Transform cannot be tranposed.")

    ## Sub is a debugging overload that compares frames of 2 transforms.
    def __sub__(self, other):
        return self.frame() - other.frame() # DBG


## \f$ \vec{v}' = \vec{C} + [\mu {\bf I} + \vec{r} {\bf \times}]\vec{v}\f$\n
#  A <a href="http://en.wikipedia.org/wiki/Helmert_transformation">
# Helmert</a> transformation:
# @image html helmert.jpg
# @image latex helmert.jpg
def helmert(c_tuple, s, r_tuple):
    """
    affine = Helmert((cx, cy, cz), s, (rx, ry, rz))
    cx, cy, cz  in meters (a Vector)
    mu in ppm
    rx, ry, rz in arc-seconds (*r as a Vector--
                               since it is a small linear_map)"""
    from math import pi
    from .tensor import DELTA
    from .vector import Vector

    C = Vector(*map(float, c_tuple))
    # note: sign is correct, since R will take an anterior product,
    # matrix representing the cross product.
    R = (-Vector(*r_tuple)*pi/180./3600.).dual()
    # mu is defined in ppm
    mu = 1.+s/1.e6

    return Affine(mu*DELTA + R, C)


## An example of a Helmert transform.
def wgs84_to_mgi():
    """return the helmert transform for MGI - untested
    and probably inverted."""
    return helmert((-577.326, -90.129, -463.920),
                   -2.423,
                   (5.137, 1.474, 5.297))
