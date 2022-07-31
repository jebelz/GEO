"""A base class for unit-converting functions."""
##\namespace geo.politics.units._bases The base unit-conversion function.
import functools
import fractions
import operator
import datetime
import math

import scipy
from scipy import constants


## A Unit wrapper?
class Unit(object):
    """unit = Unit(scale [,bias=0])

    si = unit(value)
    value = si * unit

    value = (~unit)(si)
    sq-si = (unit**2)(sq-value)
    etc..."""

    ## Construct from a scale factor and a mostly 0 bias.
    # \param scale The conversion factor
    # \param bias Usually 0, or maybe 273.15K
    def __init__(self, scale, bias=0):
        ## scale factor to SI
        self.scale = float(scale)
        ## offset of zeros (for Temp); consider pushing down to temp,
        # though that makes a circle-ellipse problem
        self.bias = float(bias)

    ## float is the scale
    # \returns Unit.scale
    def __float__(self):
        return self.scale

    ## Conversion functions
    # \return \f$ (x-b)/m \f$
    def __call__(self, x):
        return (x-self.bias) / self.scale


    def __mul__(self, other):
        if self.bias or other.bias:
            raise ValueError("Can't combine units w/ biases")
        return Unit(self.scale * other.scale)

    def __div__(self, other):
        if self.bias or other.bias:
            raise ValueError
        return Unit(self.scale / other.scale)

    def __truediv__(self, other):
        return self.__div__(other)

    def __rtruediv__(self, other):
        return self.__rdiv__(other)

    def __rmul__(self, other):
        return self.bias + other * self.scale

    def __rdiv__(self, other):
        return self(other)

    def __pow__(self, n):
        if self.bias:
            raise ValueError
        return type(self)(self.scale ** n)


## A class for coverting between two units systems:
# <a href="http://en.wikipedia.org/wiki/Fundamental_unit">Fundamental</a>
# <a href="http://en.wikipedia.org/wiki/SI">SI</a> and
# <a href="http://en.wikipedia.org/wiki/U.S._customary_units">US Customary</a>
# units based on scipy's
# <a href="http://docs.scipy.org/doc/scipy/reference/constants.html#constants-database">constants</a> package.
class SI(object):
    """SI(constant, From="a non SI unit", To="an SI unit") """

    ## Init Takes a scale factor and optional "From" and "To" keywords.
    def __init__(self, constant, From="a non SI unit", To="an SI unit"):

        ## The scale factor
        self.constant = constant
        ## The name of the "from" unit
        self.From = From
        ## The name of the
        self.To = To

    ## This method converts the argument from SI.From to SI.To
    def __call__(self, x):
        """__call__(x) will apply conversion to x,

        __call__(<SI instance">) --"> __mul__(<SI instance">) so that:

        (A(B))(x) --"> (A*B)(x) == A(B(x))
        """
        return self.__mul__(x) if isinstance(x, (SI,)) else self.constant*x

    ## This method returns an SI instance that will convert from
    # SI.To to SI.From .
    def __invert__(self):
        """self.__class__(1./float(self), self.To, self.From)"""
        return self.__class__(1./float(self), self.To, self.From)

    ## float(SI()) is SI.constant.
    def __float__(self):
        """float(self.constant)"""
        return float(self.constant)

    ## "Convert from SI.From to SI.To with scaling :SI.__float__().
    def __str__(self):
        return ("Convert from " +
                self.From+ " to " +
                self.To + " with scaling: "
                +str(float(self)))

    ## Negative Power invokes inverse, positive power does something too.
    def __pow__(self, n):
        if n == (-1):
            return ~self
        else:
            return self.__class__(float(self)**n,
                                  self.From+"**"+str(n),
                                  self.To+"**"+str(n))

    ## Stack classes.
    def __mul__(self, other, From = None, To = None):
        """(a/b)*(b/c) --"> a/c, with kwargs From=None, To=None """
        result = self.__class__(float(self)*float(other), self.From, other.To)
        if From is None and To is None:
            if 0 and self.To != other.From:
                print "Warning, possible mismath: " + self.To+ "<-->" + other.From
            return result
        result.From = From
        result.To = To
        return result

    ## Scale (or not)
    def __div__(self, other, string = ''):
        if isinstance(other, (self.__class__,)):
            return self.__class__(float(self)/float(other),
                                  self.From+"_per_"+other.From,
                                  self.To+"_per_"+other.To)
        else:
            return self.__class__(float(self)/other, self.From, string)

    def unit(self, cls=Unit):
        return cls(float(self))
