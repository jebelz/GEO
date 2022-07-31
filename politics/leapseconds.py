"""The Leap Seconds module. It has one class:

>>>datetime(2015, 6, 30, 23, 59, 60, 0)

which can be used EXACTLY like datetime.datetime, except that it knows about
leapseconds, and includes them when:

(1) Subtracting 2 leapseconds.datetime instances (do NOT mix with standard
                  datetime.datetime instances).

(2) Add/subtract-ing a standard datetime.timedelta instance.

Leap seconds cannot be predicted beyond 6 months, so y'all have to maintain
the internal tuple datetime.datetime instances at which leapseconds are
created:

PRELEAP

From this, POSTLEAP - the datetime.datetime instance immediately following a
leap second (an used internally).

Finally a tuple of leapseconds.datetime instances represented the actual
moment of a leapsecond is created:

LEAPS

You cannot create a leap second instance that is not in this tuple.
"""
## \namespace geo.politics.leapseconds
# <a href="http://en.wikipedia.org/wiki/Leap_second">leap seconds</a> aware
# <a href="http://docs.python.org/library/datetime.html#datetime-objects">
# datetime.datetime</a> class.
import functools
import datetime as DT
import operator
import itertools
from math import trunc
try:
    from pytz import timezone
except ImportError:
    pass
else:
    # UTC Time-zone info is mandatory
    UTC = timezone('UTC')


## make this available to module users
timedelta = DT.timedelta
    

## We need to extract the datetime.datetime error message for putting
# leap second data into one of those, so we can catch it an make a leap
# second.
try:
    _dummy = DT.datetime(1965, 12, 21, 1, 3, 60)
except ValueError as _err:
    MESSAGE = _err.message


## Zero time difference
ZERO = DT.timedelta(0)


## One second
ONE = DT.timedelta(0, 1)


## One Day
DAY = DT.timedelta(1)


## One Microseconds
US = DT.timedelta(0, 0, 1)


## Components of a datetime
COMPS = ('year', 'month', 'day', 'hour', 'month', 'minute', 'second',
         'microsecond')


## Leap Fail
class LeapError(ValueError, UserWarning):
    """When a leap is created that is NOT a known leap"""


## A function that is called only if the 2 arguments are equal
def _sub0(self, other):
    #pylint: disable=W0613
    return ZERO


## subtract datetimes who's difference is negative
def _subn(self, other):
    return -_subp(other, self)


## subtract 2 well ordered datetimes: let datetime to the 1st order
# calculation, and then add corrections b/c of intervening leap seconds.
def _subp(self, other):
    # farm out to builtin.
    base = self.datetime - other.datetime
    # now find leaps between endpoints for d in LEAPS (with a silly loop):
    base += len(list(itertools.dropwhile(functools.partial(
        operator.gt, other.datetime), itertools.takewhile(
            functools.partial(operator.gt, self.datetime),
            LEAPS)))) * ONE
    # return corrected time delta.
    return base + (self.leap - other.leap)*ONE


## Hash Table for subtraction operation
_SUB = {-1: _subn, 0: _sub0, 1: _subp}


## The tuple of
# <a href="http://docs.python.org/library/datetime.html#datetime-objects">
# WHEN</a> leap-seconds were added--this controls other records of leap,
# these are built-in datatime.datetime instance just before the leap.
PRELEAP = (
    DT.datetime(1972, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1972, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1973, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1974, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1975, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1976, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1977, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1978, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1979, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1981, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1982, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1983, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1985, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1987, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1989, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1990, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1992, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1993, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1994, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1995, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1997, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(1998, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(2005, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(2008, 12, 31, 23, 59, 59, tzinfo=UTC),
    DT.datetime(2012, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(2015, 6, 30, 23, 59, 59, tzinfo=UTC),
    DT.datetime(2016, 12, 31, 23, 59, 59, tzinfo=UTC))


## Coount Leap Seconds.
NLEAP = len(PRELEAP)


## This is the bultin datetime second AFTER the leap.
POSTLEAP = tuple(item + ONE for item in PRELEAP)


def _float2dt(x):
    a, b = map(int, divmod(round(1e6 * x), 1000000))
    return a * ONE + b * US


#def deb(method):
#    def dmethod(self, *args):
#        res = method(self, *args)
#        return res
#    return dmethod


## Parse Valid Leap Second Inputs.
# \param args Like a builtin datetime
# \param kwargs Like a builtin datetime
# \returns datetime for the builtin
# \throws ValueError for non-UTC time
# \throws LeapError for a LEAP that is not recognized via ::PRELEAP.
def _parse60(*args, **kwargs):
    if len(args) == 8 or 'tzinfo' in kwargs:  # raise
        raise ValueError("UTC Only")
    kwargs['tzinfo'] = UTC
    try:
        x = DT.datetime(*args, **kwargs)
    except ValueError as err:
        if err.message != MESSAGE:  # raise
            raise err
        if kwargs.get('second') is not None and kwargs['second'] == 60:
            kwargs['seconds'] -= 1
        elif len(args) >= 6 and args[5] == 60:
            args = list(args)
            args[5] -= 1
        else:
            raise err
        x = DT.datetime(*args, **kwargs)
        if DT.datetime(x.year, x.month, x.day,  # raise
                       x.hour, x.minute, x.second,
                       tzinfo=UTC) not in PRELEAP:
            raise LeapError("Not a Leap Second!")
        y = DT.datetime(*args, **kwargs), True
    else:
        y = x, False
    return y


### Count (signed/ordered) Leap Seconds between 2 endpoint (exclusive).
#def delta(a, b):
#    y, x = sorted([a, b])
#    return sum(x > lo > y for lo in LEAPS) * cmp(b, a)


## The Leap Second Aware UTC-only datetime.datetime wrapper.
@functools.total_ordering
class datetime(object):
    """See datetime.datetime documentation.

    This is that, with the allowance of leap seconds."""
    
    ## Most instances are not leap seconds.
    leap = False

    @classmethod
    def fromdatetime(cls, x):
        """construct from a normal datetime.datetime"""
        return cls(x.year, x.month, x.day, x.hour, x.minute, x.second,
                   x.microsecond)

    ## Parse standard datetime args and kwargs into:\n
    def __init__(self, *args, **kwargs):
        # farm out parsing
        dt, lp = _parse60(*args, **kwargs)

        ## Standard datetime portion of instance
        self.datetime = dt
        ## This IS a Leap Second
        self.leap = bool(lp)

    ## Truncated to seconds resolution (leap included).
    def __trunc__(self):
        return type(self)(self.year, self.month, self.day,
                          self.hour, self.minute, self.second)

    ## Leaps before self
    def __int__(self):
        return sum(item <= self for item in LEAPS)

    ## Could it be leap?
    def __nonzero__(self):
        return bool(self.datetime)

    ## Forward requests to datetime.datetime.
    def __getattr__(self, attr):
        return getattr(self.datetime, attr)

    def __eq__(self, other):
        for attr in COMPS:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def __ne__(self, other):
        for attr in COMPS:
            if getattr(self, attr) != getattr(other, attr):
                return True
        return False

    def __gt__(self, other):
        for attr in COMPS:
            x, y = getattr(self, attr), getattr(other, attr)
            if x != y:
                return x > y
        return False

    def __radd__(self, other):
        return self + other

    def __str__(self):
        result = str(self.datetime)
        if self.leap:
            result = "{}{}{}".format(result[:17], '60', result[19:])
        return result

    def __repr__(self):
        return "{}.{}({}, {}, {}, {}, {}, {}, {})".format(
            self.__module__, type(self).__name__,
            self.year, self.month, self.day, self.hour, self.minute,
            self.second, self.microsecond)

    ## The second field is either normal, or 60.
    @property
    def second(self):
        """Seconds allows 60"""
        return 60*self.leap or self.datetime.second

    ## now is UTC aware.
    @classmethod
    def now(cls):
        """Now is UTC aware"""
        return cls.fromdatetime(DT.datetime.now(UTC))

    ## add a datetime.timedelta
    def __add__(self, other):
        return next(_add4(self, other))

    ## Count Leap Seconds between 2 endpoint (exclusive).
    def __mod__(self, other):
        """n = a % b

        is the number of leapsecond that are after a and before b.

        TODO:
        {1: __rmod__,
        0: lambda x,y: 0,
        -1: sum(other > lp > self for lp in LEAPS)}cmp(self, other)
        """
        # sort arguments
        right, left = sorted([self, other])
        # count leaps inbetween and apply correct sign based on order.
        return sum(left > lp > right for lp in LEAPS) * cmp(other, self)

    ## This is invoked by datetime.datetime % leapsecond.datetime (sans if)
    def __rmod__(self, other):
        return self.fromdatetime(other) % self

    ## Subtract leapseacond from another (no datetime allowed)
    def __sub__(self, other):
        # Pass this instance over to __add__
        if not isinstance(other, type(self)):
            return self + (-other)
        # Subtract at second resolution and micro separately.
        return (
            trunc(self)._sub(trunc(other)) +
            US * (self.microsecond - other.microsecond)
        )

    ## Subtraction of 2 integer second datetimes.
    def _sub(self, other):
        assert not (self.microsecond or other.microsecond)
        return _SUB[cmp(self, other)](self, other)

    def _epoch(self, epoch, ilb=0):
        return (self.datetime +
                ONE * (int(self)-int(epoch)) +
                _float2dt(ilb * (self.datetime -
                                 epoch.datetime).total_seconds()))

    ## Convert to GPS clock
    def gps(self):
        """Convert to GPS Time."""
        return self._epoch(GPS_EPOCH)

    ## Convert to TAI
    def tai(self, _dt=10*ONE):
        """Convert to Internatinal Atomic Time."""
        return self._epoch(TAI_EPOCH) + _dt

    ## Convert to Terrestial Time
    def tt(self, _dt=(32*ONE + 184*US)):
        """Convert to Terrestial Time."""
        return self.tai() + _dt


## A tuple of actual leap seconds (created manually)
LEAPS = tuple(datetime(item.year, item.month, item.day, 23, 59, 60)
              for item in map(datetime.fromdatetime, PRELEAP))

## GPS Epoch
GPS_EPOCH = datetime(1980, 1, 6)

## TAI Epoch
TAI_EPOCH = datetime(1958, 1, 1)

## TT Epoch.
TCB_EPOCH = datetime(1977, 1, 1)
TCB_ilB = 1.5505e-8


## When solving addition problem, these are the pairs of seconds and
# 'leap' by which you might err.
_OFFSETS = zip((ZERO, ONE, ZERO), (False, False, True))


## Generate correct addition
def _add4(self, other):
    """generate correct result for datetime + timedelta. (don't worry
    about multiple solutions, as the second is the same as the 1st, and
    the 1st is always right)."""

    # leap-free answer
    base = self.datetime + other

    # nominal leap correction
    base -= ONE * (self % base)

    # Now search for a leap correction, a correction larger than 1 requires
    # for atleast 86400 leap seconds to exist-- also, you may land during
    # leap second, and that mandates a post facto correction.
    for offset, leap in _OFFSETS:
        # correct guess for offsets
        guess = self.fromdatetime(base + offset)
        # check if the result *could* be a leapsecond
        if leap and trunc(guess).datetime in PRELEAP:
            # and test both case: seconds=59 or 60.
            guess.leap = leap
        # If guess is good, kick it back.
        if  guess - self == other:
            yield guess


## dev debug func
#def ddd():
#    Left = LEAPS[:]
#    Right = LEAPS[:]
#
#    for k in (-40176001,-40176000,-40175999, -2, -1, 0, 1, 2, 40175999,40176000,40176001):
#        dd = ONE * k
#        for i, left in enumerate(Left):
#            for j, right in enumerate(Right):
#                base = right - left + dd
#                __, offset, leap = add(left, base)
#                if offset:
#                    dt = __._sus().datetime
#                    print offset, leap, left.leap, base, dt in PRELEAP, dt in POSTLEAP, dt-ONE in POSTLEAP, dt


## dev debug func
#def db1(nn):
#    t0 = LEAPS[nn]
#    for n in (-2,-1,0,1,2):
#        t = t0 + ONE*n
#        print n, int(t), t, t.datetime, t.gps(), int((t.gps()-t.datetime).total_seconds())
