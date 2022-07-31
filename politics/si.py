"""objects with units. There are 7 module constants, 1 for each of the 7
base dimensions in the SI system. Each one is an instance of the corresponding


class, instansiated with a value of 1. (That is, they are unit units). They
are the base units, but they are not base classes. Each base unit inherits
from the same super classes, and have the same metaclass. Here is a table:

SI Base Units:
|------------|------------|-----------------------------------------|
|instance[1] |  class[2]  |   Dimensionality  (a named tuple)[3,4]  |
|------------|------------|-----------------------------------------|
|   m        |   Meter    | (L=1, M=0, T=0, I=0, Theta=0, N=0, J=0) |
|   kg       |  Kilogram  | (L=0, M=1, T=0, I=0, Theta=0, N=0, J=0) |
|   s        |   Second   | (L=0, M=0, T=1, I=0, Theta=0, N=0, J=0) |
|   A        |   Ampere   | (L=0, M=0, T=0, I=1, Theta=0, N=0, J=0) |
|   K        |   Kelvin   | (L=0, M=0, T=0, I=0, Theta=1, N=0, J=0) |
|   mol      |    Mole    | (L=0, M=0, T=0, I=0, Theta=0, N=1, J=0) |
|   cd       |   Candela  | (L=0, M=0, T=0, I=0, Theta=0, N=0, J=1) |
|------------|------------|-----------------------------------------|
[1] _BASEUNIT is a nameddtuple with the instance assigned to the correct
field, e.g

_BASEUNIT.L is m

[2] _BASECLASS is a nameddtuple with the instance assigned to the correct
field, e.g

_BASECLASS.N is Mole

[3, 4] SYMBOL and _NAME are namedtuples with string fields assigned
       representing the symbol (default instance name) and name (default
       class name).

Why: Because D.   on't      is beter than W. rite
             R.   epeat                   E. verything
             Y.   ourself                 T. twice...to... ten thousand time

That's where the metaclass comes in:

Unit instances are made by their classes, and unit classes are made by
their metaclass:

|-------------|----------------|---|----------------|-------------------|
|  Unit Type  | Module Constant| # | Base Classes   |    metaclass      |
|=============|================|===|================|===================|
| Base Unit   |   _BASECLASSES  | 7 | -> Base ->_SI  | _base_metaclass   |
|-------------|----------------|---|----------------|-------------------|
|Derived Unit |   None         | 0 |-> Derived ->_SI|_derived_metaclass |
|=============|================|===|================|===================|

Say what? Why metaclasses: because the make the classes right. When you
load a module, the 7 base-unit classes are built by type.__new__ and then
extended by _base_metaclass.__new__ as follows:

Their dimensionality is computed from the global variables describing
the base dimensions. The root class _SI holds the global dictionary
of units (_SI.units), and that is updated whenever a class is made.

The key is the dimensionality tuple, and the value is the class.

So now you see why there ARE NO DERIVED CLASSES. They get made when you
need them, and that is done in __invert__ or __mul__ (that is, unit
inversion or unit multiplication)--all other algebraic operations that
create new units farm the work out to these methods. So for example:

g = 9.8*m/s**2

creates new classes in the dictionary, via:
s**2 --> s*s --> s.__mul__(s) ==> Second2
/(s**2) -> ~(s2) ==> s2.__invert__() ==> perSecond2
m(/s**2) --> m*perSecond2 --> m.__mul__(perSecond2) --> MeterperSecond2

In [9]: for key, value in m.units.iteritems():
   ...:     print key, value
   ...:
Dimensionality(L=0, M=0, T=0, I=0, Theta=0, N=0, J=1) <class 'si.Candela'>
Dimensionality(L=0, M=0, T=-2, I=0, Theta=0, N=0, J=0) <class 'si.perSecond2'>
Dimensionality(L=0, M=0, T=1, I=0, Theta=0, N=0, J=0) <class 'si.Second'>
Dimensionality(L=0, M=0, T=2, I=0, Theta=0, N=0, J=0) <class 'si.Second2'>
Dimensionality(L=1, M=0, T=0, I=0, Theta=0, N=0, J=0) <class 'si.Meter'>
Dimensionality(L=0, M=0, T=0, I=0, Theta=1, N=0, J=0) <class 'si.Kelvin'>
Dimensionality(L=0, M=1, T=0, I=0, Theta=0, N=0, J=0) <class 'si.Kilogram'>
Dimensionality(L=0, M=0, T=0, I=1, Theta=0, N=0, J=0) <class 'si.Ampere'>
Dimensionality(L=0, M=0, T=0, I=0, Theta=0, N=1, J=0) <class 'si.Mole'>
Dimensionality(L=1, M=0, T=-2, I=0, Theta=0, N=0, J=0) <class 'si.MeterperSecond2>

These new classes are derived classed, and are created and put into the
dictionary by the _derived_metaclass meta-class.

Unit classes have one dynamic attribute: "value" -- it can be anything the
represents one or many numbers (e.g. a scalar or an array). They have all the
basic funcitons overloaded:

complex
float
int
long
trunc
len
iter

Finally: all classes have a "to_imperial()" method that will convert their
value into a dimensionless value that is/are correct in imperial units. Base
units are converted with the static class attribute "correction" -- which is a
correction factor to the units into American. For derived units, the
correction factor is computed on the fly. E.g:

In [14]: print g.to_imperial()
32.1522309711

(that's in ft/sec**2)

In [15]: print g.correction
3.28083989501

(i.e. ft/s**2 /  (m/s**2) )


Finally: don't expressed dimensionally silly equations, that *WILL* be
rejected, and a UnitsError will be thrown."""

## \namespace geo.politics.si
# <a href="http://en.wikipedia.org/wiki/International_System_of_Units">
# Dimensionally Aware </a> numbers.

import collections
import functools
import itertools
import operator

## All the base units.
__all__ = ('m', 'kg', 's', 'A', 'K', 'mol', 'cd')


## Unit Dimensions Error
class UnitError(UserWarning, TypeError):
    """Apples per seconds versus Oranges per Kelvin"""


## A collection of the 7
# <a href="http://en.wikipedia.org/wiki/International_System_of_Units#Base_units>fundamental dimensions</a>
_Dim = collections.namedtuple('Dimensionality', 'L M T I Theta N J')


## Todo: promote all dimension math to this class
class Dim(_Dim):
    """Dim as a dimension."""

    def __invert__(self):
        return type(self)(*map(operator.neg, self))

    def __mul__(self, other):
        return type(self)(*itertools.starmap(operator.add, itertools.izip(
            self, other)))

    def __div__(self, other):
        return self + ~other

    def __pow__(self, n):
        return type(self)(*[n*item for item in self])

## Symbols of Base Units
_SYMBOL = _Dim('m', 'kg', 's', 'A', 'K', 'mol', 'cd')

## Name of Base Units
_NAME = _Dim('Meter',
             'Kilogram',
             'Second',
             'Ampere',
             'Kelvin',
             'Mole',
             'Candela')

## 7.
NUMBER_OF_BASE_DIMENSIONS = len(_SYMBOL)


## Do meters plus kilograms? You get this:
class UnitsError(Exception):
    """raised when units don't match dimensons"""


#pylint: disable=R0201
## decorator to apply method to value
def _value(func):
    """_value(func) make a method return:
    func(self.value)
    """

    @functools.wraps(func)
    def wrapped_method(self):
        return func(self)(self.value)

    return wrapped_method


## decorator for methods that return something with dimensions
def _physical(func):
    """_physical(func) makes a method return
    type(self)(func(self.value(*args)))
    """

    @functools.wraps(func)
    def wrapped_method(self, *args):
        return type(self)(func(self)(self.value, *args))

    return wrapped_method


## decorator for rich comparison
def _cmp(method):
    func = method(None, None)

    def new_method(self, other):
        if not isinstance(self, type(other)):
            msg = "Can't compare {} with {}."
            raise UnitError(msg.format(type(self), type(other)))
        return func(self.value, other.value)

    return new_method


## The base class of all united objects
class _SI(object):
    """base class for all SI object"""

    _dimensions = None

    name = None

    symbol = None

    ## default correction factor for imperial units
    correction = 1.

    ## A dictionary of known units, gets updated when new
    # units are created.
    units = {}

    ## An instance is a number, or list of-, or array-of, with units--
    # hence the generic name "value"
    def __init__(self, value):
        ## The quantity of what you have
        self.value = value

    ## if value has len, then so does self
    @_value
    def __len__(self):
        return len

    ## if value has float, then so does self
    @_value
    def __float__(self):
        return float

    ## if value has complex, then so does self
    @_value
    def __complex__(self):
        return complex

    ## if value has int, then so does self
    @_value
    def __int__(self):
        return int

    ## if value has long, then so does self
    @_value
    def __long__(self):
        return int

    ## if value has trunc, then so does self
    @_value
    def __trunc__(self):
        from math import trunc
        return trunc

    ## if value is an interator, then so is self
    @_value
    def __iter__(self):
        return iter

    ## if, then....
    @_value
    def __abs__(self):
        return abs

    ## If value can be indexed, then so call self,
    # and it returns another instance of type(self)
    @_physical
    def __getitem__(self):
        from operator import getitem as result
        return result

    ## Container behavior continues
    @_physical
    def __setitem__(self):
        from operator import setitem as result
        return result

    ## Container behavior continues
    @_physical
    def __delitem__(self):
        from operator import delitem as result
        return result

    ## +x makes a (shallow) copy (see neg)
    @_physical
    def __pos__(self):
        from operator import pos as result
        return result

    ## -x --> x.__class__(-x.item)
    @_physical
    def __neg__(self):
        from operator import neg as result
        return result

    ## x + y better make sense! kg+kg, A+A, etc...
    def __add__(self, other):
        try:
            if self._dimensions != other._dimensions:
                message = "cannot add (sub) %s to (from) %s" % (other.name,
                                                                self.name)
                raise UnitsError(message)
        except AttributeError:
            message = (
                "cannot do arithmitic in %ss without another %s!" %
                (2 * (self.name,))
                )
            raise UnitsError(message)
        return type(self)(self.value + other.value)

    ## x - y --> x + (-y)
    def __sub__(self, other):
        return self.__add__(-other)

    ## a-x --> (-x) + a
    def __rsub__(self, other):
        return self.__neg__() + other

    ## This should only be called with a nonsense equation
    def __radd__(self, other):
        message = (
            "cannot do arithmitic in %ss without another %s!"
            %
            (2 * (self.name,))
            )
        raise UnitsError(message)

    ## x*alpha or x*t:
    # you can scale via __rmul__
    # or make new derived units in the
    # module function: multiply()
    def __mul__(self, other):
        if isinstance(other, _SI):
            dds = tuple([sum(item) for item in zip(self._dimensions,
                                                   other._dimensions)])
            # chose existing or make a new cls on the fly
            cls = _SI.units[dds] if dds in _SI.units else _derived_factory(dds)
            return cls(self.value*other.value)
        return self.__rmul__(other)

    ## alpha*x --> x.__class__(alpha*self.value)
    def __rmul__(self, other):
        return type(self)(other*self.value)

    ## x/alpha or x/y --> x*(~y):
    # you can scale, or make new derived units, the latter lets
    # __mul__ and __invert__ do the work
    def __div__(self, other):
        if isinstance(other, _SI):
            if self._dimensions == other._dimensions:
                result = self.value/other.value
            else:
                result = self*(~other)
        else:
            result = type(self)(self.value/other)
        return result

    ## alpha/x --> alpha*(~x)
    def __rdiv__(self, other):
        return self.__invert__()*other

    ## Invert an _SI arguemnt in 3 steps:                          \n
    # (1) get inverse-dimension tuple                            \n
    # (2) get a class for that tuple (from COTS or from factory) \n
    # (3) instansiate the class with the multiplicative invese   \n
    def __invert__(self):
        from operator import neg
        dds = tuple(map(neg, self._dimensions))
        # note: should move this choice into the metaclass--as it is what
        # knows whcih units exists.
        cls = _SI.units[dds] if dds in _SI.units else _derived_factory(dds)
        return cls(1./self.value)

    ## pow is for integer arguments only, for n>0 (n<0):       \n
    #  x**n --> reduce(mul, [x]*n) --> x*x*...n-times...*x*x  \n
    # (x**n) --> (~x)**(-n)
    def __pow__(self, n):
        if n == 0:
            result = self.__mul__(self.__invert__())
        elif n < 0:
            result = self.__invert__().__pow__(-n)
        else:
            from itertools import repeat
            from operator import mul
            result = reduce(mul, repeat(self, n))
        return result

    ## String includes a symbol
    def __str__(self):
        return str(self.value) + self.symbol

    ## repr is set here for development
    def __repr__(self):
        return str(self)

    def __nonzero__(self):
        return bool(self.value)

    @_cmp
    def __eq__(self, other):
        from operator import eq as result
        return result

    @_cmp
    def __ne__(self, other):
        from operator import ne as result
        return result

    @_cmp
    def __gt__(self, other):
        from operator import gt as result
        return result
    @_cmp
    def __ge__(self, other):
        from operator import ge as result
        return result
    @_cmp
    def __lt__(self, other):
        from operator import lt as result
        return result

    @_cmp
    def __le__(self, other):
        from operator import le as result
        return result

    @_cmp
    def __contains__(self, other):
        from operator import contains as result
        return result

    ## Convert to imperial units
    def to_imperial(self):
        """convert value to imperial units (magnitude)"""
        return self.value*self.correction

    ## New classes call this to notify -_SI.units about their existence
    # though it should be in metaclass, since metaclass knows all.
    # \param cls
    # \sideffect modifies cls.units
    @staticmethod
    def update(cls):
        """update(cls) is a static method"""
        if cls._dimensions not in cls.units:
            cls.units.update({cls._dimensions: cls})
        return cls


## A class factor for derived units
def _derived_factory(dds):
    """cls = _derived_factoy(dds)

    dds  is a sequence of 7 dimensional exponenets defining the dimensions
         of the derive unit-- see global _Dim for ordering.

    cls  is a new sub-class _DerivedClass-- it is put in the _SI.units
    dictionary."""
    new_name = ''
    per_name = ''
    New_name = ''
    Per_name = ''

    # compute class's name
    for d, b, n in zip(dds, _SYMBOL, _NAME):
        if not d:
            continue
        if d > 0:
            new_name += n
            New_name += n
            if d > 1:
                new_name += str(d)
                New_name += str(d)
        else:
            per_name += n
            Per_name += n
            if d < -1:
                print per_name, -d
                per_name += str(-d)
                Per_name += str(-d)
    if per_name:
        new_name = new_name+"_per_"+per_name
        New_name = New_name+"_per_"+Per_name
    ## OK: name is computed, now make a new class for it:
    return _derived_metaclass(New_name,
                              (_Derived,),
                              {'_dimensions': dds,
                               'name': New_name,
                               'symbol': new_name})


## This metaclass updates the _SI.units dictionary whenever
# a new units class is made
class _units_metaclass(type):

    ## this throw a type error-I don't get it.
    def update(self, cls):
        return _SI.update(cls)


## Metaclass for derived methods from the factory:
# it updates the dimension dictionary with a
# _Dim tuple key.
class _derived_metaclass(_units_metaclass):

    def __new__(mcs, *args, **kwargs):
        cls = type.__new__(mcs, *args, **kwargs)
        cls._dimensions = _Dim(*(cls._dimensions))
        return _SI.update(cls)


## THis makes code run like a script --_BASECLASSES starts empty, and is filled
# with base unit classes as the metaclass makes them--weird.
_BASECLASSES = []


## This metaclass uses the class name to figure out the Dimensionality-
# it also puts the class in _SI.units dictionary of units; the alternative is
# a decorator.
class _base_metaclass(_units_metaclass):

    def __new__(self, *args, **kwargs):
        cls = type.__new__(self, *args, **kwargs)
        cls.name = cls.__name__.lstrip('_')
        try:
            idx = _NAME.index(cls.name)
            cls.symbol = _SYMBOL[idx]
            dimensions = [0]*NUMBER_OF_BASE_DIMENSIONS
            dimensions[idx] = 1
            cls._dimensions = _Dim(*dimensions)
            _BASECLASSES.append(_SI.update(cls))
        except ValueError:
            pass
        return cls


## A base class for derived units-- doesn't do anything yet,
# and it may go away; original purpose was to have the
# metaclass for derived classes, but that did not work.
class _Derived(_SI):

    ## A computed read-only attribute.
    @property
    def correction(self):
        """Unknown"""
        from operator import mul, attrgetter
        return reduce(
            mul,
            [b**a for a, b in
             zip(self._dimensions, map(attrgetter('correction'), _BASECLASSES))]
            )


## Base class for base units- it assigns their
# meta class: _base_metaclass
class _Base(_SI):

    ## All base units use this metaclass
    __metaclass__ = _base_metaclass


class _Meter(_Base):

    ## m -> ft
    correction = (100./2.54/12.)


class _Kilogram(_Base):

    ## kg -> lbs
    correction = 2.20462


class _Second(_Base):
    pass


class _Ampere(_Base):
    pass


class _Kelvin(_Base):

    def to_imperial(self):
        return 32 + 1.8*(self.value-273.15)


class _Mole(_Base):
    pass


class _Candela(_Base):

    ## cd -> cp (Candelpower)
    # http://en.wikipedia.org/wiki/Candlepower
    correction = 1/0.981


## collect the classes in a _Dim tuple-- recall _BASECLASSES was built
# by the _base_metaclass when it made the above 7 classes.
_BASECLASSES = _Dim(*_BASECLASSES)

## intansiate _BASECLASSES with unit argument and put them in a Dum tuple
_BASEUNIT = _Dim(*(cls(1.) for cls in _BASECLASSES))

## Unpack individual units
m, kg, s, A, K, mol, cd = _BASEUNIT
