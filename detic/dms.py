"""Two functions for converting dd.dddddd <--> DD MM' SS.ssss"

degree2dms
dms2degree


The DMS class does arithmetic with DDMMSS tuples.


There is a major architecural problem:
For a DMS object x:

cos(x) works, but cosd(x) does not, as it does cos(radians(x))

"""
## \namespace geo::detic::dms Degree, Minutes, Seconds <--> Decimal conversion
import functools
import math

import numpy as np

__all__ = ('ddmmss2degree', 'degree2ddmmss', 'DMS')


## Convert method argument to degrees
def convert(method):
    """Decorator converts arguments to degrees, e.g:

    >>>cosd = convert(cos)
    """
    #pylint: disable=C0111
    @functools.wraps(method)
    def converted_method(self, other):
        try:
            other = other.degrees()
        except AttributeError:
            pass
        return method(self, other)
    return converted_method


## Revert method results back to class to which method is bound
def revert(method):
    """decorator restore convert in arc functions"""
    #pylint: disable=C0111
    @functools.wraps(method)
    def reverted_method(self, *args):
        return type(self)(*degree2ddmmss(method(self, *args)).tolist())
    return reverted_method


## A Degree, Minute, Seconds object that permits +, -, *, /
class DMS(object):
    """DMS(dd, [mm=0, ss=0.])

    dd   degrees
    mm   minutes
    ss   seconds
    """
    def __init__(self, dd, mm=None, ss=None):
        ## Floor of degrees
        self.dd = dd
        ## Floor of minutes
        self.mm = 0*self.dd if mm is None else mm  # kwd
        ## Decimal seconds
        self.ss = 0*self.mm if ss is None else ss  # kwd

    def __len__(self):
        try:
            result = len(self.dd)
        except TypeError:
            try:
                result = len(self.mm)
            except TypeError:
                try:
                    result = len(self.ss)
                except TypeError:
                    result = NotImplemented
        return result

    def __str__(self):
        return """(%s, %s', %s")""" % tuple(map(str, self.tolist()))

    __repr__ = __str__

    def tolist(self):
        """[dd, mm, ss.ssssss]"""
        return [self.dd, self.mm, self.ss]

    def __getitem__(self, index):
        return type(self)(*[item[index] for item in self.tolist()])

    def __setitem__(self, index, value):
        self.dd[index] = value.dd
        self.mm[index] = value.mm
        self.ss[index] = value.ss

    ## not valid for array_like data;
    def __float__(self):
        return self.degrees()

    ## Convert to degrees
    def degrees(self):
        """convert to decimal degrees"""
        print "deg"
        return ddmmss2degree(*self.tolist())

    deg = degrees

    ## radians, If you must
    def radians(self):
        """radians...."""
        print "rad"
        return math.pi*self.degrees()/180.

    @convert
    def __eq__(self, other):
        return self.degrees() == other

    @convert
    def __gt__(self, other):
        return self.degrees() > other

    @convert
    def __lt__(self, other):
        return self.degrees() < other

    @convert
    def __ne__(self, other):
        return self.degrees() != other

    @convert
    def __ge__(self, other):
        return self.degrees() >= other

    @convert
    def __le__(self, other):
        return self.degrees() <= other

    @revert
    def __abs__(self):
        return abs(self.degrees())

    @revert
    def __neg__(self):
        return -self.degrees()

    ## force positive minutes and seconds
    @revert
    def __pos__(self):
        return self.degrees()

    @revert
    @convert
    def __add__(self, other):
        return self.degrees() + other

    def __sub__(self, other):
        return self.__add__(-other)

    @revert
    def __mul__(self, x):
        return x*self.degrees()

    def __div__(self, x):
        return self.__mul__(1./x)

    ## This overload allows np.cos, np.sin, etc.. to work automatically.
    def __array__(self):
        try:
            len(self)
        except TypeError:
            return np.array(self.radians())
        return self.radians()

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __rmul__(self, x):
        return self.__mul__(x)


## Basic singelton function
def _deg(d, m=0, s=0.):
    """dd.ddd = ddmmss2degree(d [, m=0 [, s=0.]])"""
    f = math.copysign(1, d)
    return d + f*(m + s/60.)/60.


## Basic singelton function
def _dms(deg):
    """dd, mm, ss = degree2ddmmss(degree)"""

    f = math.copysign(1, deg)

    dd, m_ = divmod(abs(deg), 1)
    dd = int(dd)

    mm, ss = divmod(60*m_, 1)
    mm = int(mm)

    ss *= 60.

    return f*dd, mm, ss


## Vectorized private functions
_ddmmss2degree = np.vectorize(_deg)
## Vectorized private functions
_degree2ddmmss = np.vectorize(_dms)


## ufunc style function
def ddmmss2degree(dd, mm, ss):
    """(ddd, mm, ss.ssss) --> decimal degrees"""
    deg = _ddmmss2degree(dd, mm, ss)
    try:
        deg = float(deg)
    except TypeError:
        pass
    return deg


## ufunc style function
def degree2ddmmss(deg):
    """decimal degrees --> DMS"""
    dd, mm, ss = _degree2ddmmss(deg)
    try:
        dd = int(dd)
        mm = int(mm)
        ss = float(ss)
    except TypeError:
        dd = dd.astype(int)
        mm = mm.astype(int)

    return DMS(dd, mm, ss)
