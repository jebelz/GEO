"""interval arithmetic:

A = Interval(start, stop)

who's single responsibility is to represent a closed interval on an ordered
field.



Basically arithmetic is defined according standard definitions
(http://en.wikipedia.org/wiki/Interval_arithmetic)

A + B
A - B
A * B
A / B


Rich comparison For 2 intervals A, B is defined as follows:

A < B          A is totally earlier/less than B
A <= B         A is less/earlier than, but overlaps B from the left
A == B         total equality
A !=B          not A == B
A >=B          A is larger/later than B, but overlaps it form the right
A > B          A is totally larger/later than B

For field elements, x, y:
A < x          x is above A
A <= x         A < x or x is at A's top limit
A > x          x is below A
A >= x         A > x or x is at A's lower limit


See operator docstring For definition of arithmetic operations:

+A             indentity
-A             reflect boundaries, preseving order
~A             flip order
B in A, x in A ....interval is contain in larger one (or not)
A & B          min interval overlapping A and B
A | B          min interval containing all if A and B
A + B
A - B
A >> x         shift right boundary
A << x         shift left boundary

A ^ B
A % B

Standard Functions:

int(A)     --> +/- is the direction of the interval TODO: make int(float)
float(A)
abs(A)
bool(A)

Iterator/Container behavior:

iter(A)
A[index]

Non-Elementary Function:
------------------------
Methods:

support, diameter, center, radius              # properties of the interval

intersects (overlaps), is_covered_by, filter   # relation to other interval

Module Functions:
length                is NOT len(). (len undefined, even though it's 2).

sort_start            # sorting lists of intervals
sort_stop
sort_length

coverage              # coverage of a list of intervals
"""
## \namespace camp2ex.utils.interval Intervals over
# \f$\mathbb R\f$ and UTC

#pylint: disable=C0103
import abc
import operator
import itertools

## Docstring getter
docstr = operator.attrgetter("__doc__")

## Get an intervals length-- less mathy than support, negative is allowed.
length = operator.methodcaller("support")

minmax = (min, max)


## Sort some number of list by their left endpoint
#  \param list_ A list of Interval objects
#  \retval list_ sorted by Interval.start
def sort_start(*args):
    """sort_start(i1, i2, ...) sorts on i_n.start value."""
    result = list(args)
    result.sort(lambda i1, i2: cmp(i1.start, i2.start))
    return result


## Sort some number of list by their right endpoint.
#  \param list_ A list of Interval objects
#  \retval list_ sorted by Interval.stop
def sort_stop(*args):
    """sort_start(i1, i2, ...) sorts on i_n.stop value."""
    result = list(args)
    result.sort(lambda i1, i2: cmp(i1.stop, i2.stop))
    return result


## Sort some number of list by their length (rule of 3?)
#  \param list_ A list of Interval objects
#  \retval list_ sorted by ::length
def sort_length(*args):
    """sort_start(i1, i2, ...) sorts on length(i_n)"""
    result = list(args)
    result.sort(lambda i1, i2: cmp(length(i1), length(i2)))
    return result


## is this coverage?
def coverage_dev1(*args):
    """reduce on xor"""
    return reduce(operator.xor, sort_start(*args))


## Minimal continuous coverage of args, from the left.
#  \param *args Any number of Interval objects.
#  \retval interval of args continuous coverage from the left.
def coverage(*args):
    """interval = coverage(interval1, ..., intervalN)

    is the continuous coverage (from the left).
    """
    result = ()
    for n, item in enumerate(sort_start(*args)):
        if not n:  # initializer
            result = item
        else:
            new_result = result ^ item
            if new_result is None:  # guard condition
                return result
            else:
                result = new_result

    # result or empty tuple
    return result


## Warning the Interval is reversed, which is ok in topological algebra
#  but here, not so much.
class ReversedIntervalWarning(UserWarning):
    """TBD"""


## Standard Python Behaviors
class _Ob(object):
    """Define pythonic stuff:
    __init__
    __iter__
    __str__
    __getitem__

    """

    _fields = ('start', 'stop')

    ## just a mix-in
    __metaclass__ = abc.ABCMeta

    ## \param start starting point
    #  \param stop ending point
    def __init__(self, start, stop):
        ## Left end point
        self.start = start
        ## Right end point
        self.stop = stop

    ## Intervals are 2ples
    def __len__(self):
        return 2

    ## Closed interval string, with normal ordering and arrow reversing.
    def __str__(self):
        if int(self) >= 0:
            result = "[" + " =>> ".join(map(str, self)) + "]"
        else:
            result = "[" + " <<== ".join(map(str, ~self)) + "]"
        return result

    __repr__ = __str__

    ## \returns generator yields self[0], self[1]
    def __iter__(self):
        """iterate endpoints in order"""
        yield self[0]
        yield self[1]

    ## Container of endpoints
    #  \param index in (-1, 0, 1)
    #  \retval element of the field over which Interval is defined
    def __getitem__(self, index):
        """[start, stop]"""
        return getattr(self, self._fields[index])


## Diagnostics Mix-in
class _Dx(object):
    """Interval Diagnostics:
    int(), float(), abs(), bool()

    and:
    support, radius, center, ...
    """

   ## just a mix-in
    __metaclass__ = abc.ABCMeta

    ## sgn(interval)
    def __int__(self):
        """-1  Reversed ordered
        0   Degenerate
        1  Proper, Normal ordered"""
        import math
        return int(math.copysign(1, float(self))) if self else 0  # guard

    ## float of support() or _to_float()
    #  \retval float
    def __float__(self):
        """Try to convert support to a float"""
        try:
            return float(self.support())
        except TypeError:
            return float(self._to_float())

    ## | float |
    #  \retval float
    def __abs__(self):
        """abs(float(self))"""
        return abs(float(self))

    # pylint: disable=R0201
    ## TBD float (iF, e.g., support is a datetime.timedetla)
    #  \throws TypeError until user defines it.
    def _to_float(self):
        """_to_float is NotImplemented. By default, you don't need it.
        But: If your interval is defined over a field that does not
        support "float", the you can override this method to define
        some sort of python-float to the elements of the field."""
        raise TypeError("Interval's field doesn't have a float defined")

    ## Size of Interval's support
    #  \retval element diFference
    def support(self):
        """DiFference of endpoints"""
        return self[1]-self[0]

    ## \retval a positive (wrt to field) diFference
    #  \throws TypeError If the field doesn't support the notion of positive
    #  and negative dIfferences.
    def diameter(self):
        """ diameter is a positive measure of the interval"""
        # get the +/- support
        result = self.support()
        # Try to decide sign- this formulation works on 'datetime' objects.
        is_pos = result >= result*0
        # general center finder -not just real numbers
        return result if is_pos else -result  # guard

    ## Degeneracy
    #  \retval bool If non-degenerate
    def __nonzero__(self):
        """non zero support"""
        return bool(float(self))

    ## Center
    #  \retval Element at center of Interval
    def center(self):
        """center = start + radius()

        Note: start + (start-stop)/2 is computable For datetimes
        while (start+stop)/2 is not- hence, the former.

        Note: use min(self), not self.start, so it works For negative
        intervals. radius is always positive."""
        return min(self) + self.radius()

    ## Radius
    #  \retval radius diFference()/2
    #  \throws TypeError If field dIfference doesn't support __div__,
    #  __mul__(float) or __rmul__(float)
    def radius(self):
        """radius is half of the diameter"""
        # Try __div__(2)
        diameter = self.diameter()
        try:
            return diameter/2
        except TypeError:
            pass
        # Try __mul__(0.5)
        try:
            return diameter*0.5
        except TypeError:
            pass
        # Force __rmul__(0.5).
        return 0.5*diameter

    ## Divvy up interval into equal boundaries
    # \param n Integer number of slices
    # \throws TypeError
    # \yields points
    def ilinspace(self, n):
        try:
            xgen = xrange(n+1)
        except TypeError as err:
            raise err("Expected int, got: {}".format(type(n)))
        start = self[0]
        end = self[1]
        da = (end - start) / float(n)
        return (count * da for count in xgen)

    def divvy(self, n):
        from sys import maxint
        points = self.ilinspace(n)
        g1, g2 = itertools.tee(points, 2)  # does this burn memory?
        s1 = itertools.islice(g1, 0, maxint)
        s2 = itertools.islice(g2, 1, maxint)
        return (type(self)(a, b) for a, b in
                itertools.izip(s1, s2))


## <a href="http://en.wikipedia.org/wiki/Interval_arithmetic">Standard
#  Interval Arithmetic</a>
class _Arithmetic(object):
    """Official Interval Arithmetic
    C = A + B
    C = A - B
    C = A * B
    C = A / B

    These are mostly used For fuzzy work with a measurement and it's
    error-bar defined interval.
    """

    ## just a mix-in
    __metaclass__ = abc.ABCMeta

    ## Identity
    #  \retval interval \f$ +[x, y] \equiv [+x, +y] \f$
    def __pos__(self):
        """+[x, y] --> [+x, +y]"""
        return self.unary(operator.pos)

    ## Reflection about 0
    #  \retval interval  \f$ -[x, y] \equiv [-x, -y] \f$
    def __neg__(self):
        """-[x, y] --> [-x, -y]"""
        return self.unary(operator.neg)

    ## Add endpoints, pairwise (non-simplIfied computation)
    #  \retval interval
    #  \f$ [x, y] + [x', y'] \equiv [x+x', y+y'] \f$
    def __add__(self, other):
        """[x, y] + [x', y'] --> [x+x', y+y']"""
        return self.binary(other, operator.add)

    ## Subtract Endpoints, pairwise
    #  \param other an Interval
    #  \retval interval
    #  \f$ [x, y] - [x', y'] \equiv [x-x', y-y'] \f$
    def __sub__(self, other):
        """[x, y] - [x', y'] --> [x-x', y-y']"""
        return self.binary(other, operator.sub)

    ## Interval Multiply
    #  \param other an Interval
    #  \retval interval
    #  \f$ [x, y] * [x', y'] \equiv
    #  [\min{(xx', xy', yx', yy')}, \max{(xx', xy', yx', yy')}] \f$
    def __mul__(self, other):
        """[a, b] * [c, d] =
        [min (a*c, a*d, b*c, b*d), max(a*c, a*d, b*c, b*d)]"""
        return self.binary(other, operator.mul)

    ## Interval Divide
    #  \param other an Interval
    #  \retval interval \f$
    #  [\min{(x/x', x/y', y/x', y/y')}, \max{(x/x', x/y', y/x', y/y')}] \f$
    def __div__(self, other):
        """[a, b] / [c, d] =
        [min (a/c, a/d, b/c, b/d), max (a/c, a/d, b/c, b/d)]"""
        return self.binary(other, operator.div)

    ## TBD- there are nuances that are not treated here.
    def __pow__(self, n):
        return reduce(operator.mul, itertools.repeat(self, n))

    ## Apply arbitrary unary function
    #  \param func a unary operator
    #  \retval interval \f$ [f(x), f(y)] \f$
    def unary(self, func):
        """[x, y] --> [func(x), func(y)]"""
        return type(self)(*map(func, self))

    ## Binary function with [min(F), max(F)]
    #  \param other an Interval
    #  \param func a binary operator
    #  \retval interval
    #   \f$ [\min{(f(x,x'), f(x,y'), f(y,x'), f(y,y'))},
    #  \max{(f(x,x'), f(x,y'), f(y,x'), f(y,y'))}] \f$
    def binary(self, other, func):
        """f([x, y], [x', y'] -->
        [min(f(x, x'), f(x, y'), f(y, x'), f(y, y'),
        max(f(x, x'), f(x, y'), f(y, x'), f(y, y'),]
        """
        return type(self)(*[
            ext(itertools.starmap(func, itertools.product(self, other)))
            for ext in (min, max)])


## Boolean functions comparing intervals and/or elements.
class _Comparator(object):
    """Interval-to-Interval comparison (with Boolean results)"""

    ## just a mix-in
    __metaclass__ = abc.ABCMeta

    ## Set-like "in"
    #  \param other Interval or element
    #  \retval bool \f$ A \in B \iff  x \in B \ \forall \ x \in A \f$
    def __contains__(self, other):
        """a in [x, y] --> x <= a <= y \n
        [x', y'] in [x, y] --> x' in [x, y] and y' in [x, y]"""
        if isinstance(other, type(self)):  # parse inputs
            return all(map(self.__contains__, other))
        return min(self) <= other <= max(self)

    ## [x, y] == [x', y'] --> x == x' and y == y'
    #  \param other Interval or element
    #  \retval bool \f$ x \equiv x'\wedge y \equiv y' \f$
    def __eq__(self, other):
        """total equality"""
        return self[0] == other[0] and self[1] == other[1]

    ## [x, y] != [x', y'] --> not [x, y] == [x', y']
    #  \param other Interval or element
    #  \retval bool \f$ x \ne x' \vee    y \ne y' \f$
    def __ne__(self, other):
        """not total equality """
        return self[0] != other[0] or self[1] != other[1]

    ## Less Than
    #  \param other Interval or element
    #  \retval bool \f$ y < x' \f$
    def __lt__(self, other):
        """[x, y] >= [x', y'] --> y < x'
        [x, y] >= a --> x > a and y > a"""
        if isinstance(other, type(self)):  # parse inputs
            return self[1] < other[0]
        return max(self) < other

    ## Overlap Left
    #  \param other Interval or element
    #  \retval bool \f$ x < x' \land y \in [x', y'] \f$
    def __le__(self, other):
        """[x, y] <= [x', y'] --> y in [x', y'] ' and x < x'
        [x, y] >= a --> x > a and y > a"""
        if isinstance(other, type(self)):  # parse inputs
            return self[0] < other[0] and self[1] in other
        return max(self) <= other

    ## Greater Than
    #  \param other Interval or element
    #  \retval bool \f$ x > y' \f$
    def __gt__(self, other):
        """[x, y] > [x', y'] --> x > y'
        [x, y] > a --> x > a and y > a"""
        if isinstance(other, type(self)):  # parse inputs
            return self[0] > other[1]
        return min(self) > other

    ## Overlaps Right
    #  \param other Interval or element
    #  \retval bool \f$ y > y' \land y \in [x', y'] \f$
    def __ge__(self, other):
        """[x, y] >= [x', y'] --> x in [x', y'] ' and y > y'
        [x, y] >= a --> x > a and y > a"""
        if isinstance(other, type(self)):  # parse inputs
            return self[0] in other and self[1] > other[1]
        return min(self) > other

    ## Any overlap at all?
    #  \param other Interval
    #  \bool \f$ (A \le B) \vee (A \ge B) \vee (A\in B) \vee (B \in A) \f$
    def intersects(self, other):
        """True IFF  self intersects other in anyway"""
        return (
            any(map(self.__contains__, other)) or
            any(map(other.__contains__, self))
            )

    ## maybe a better name.
    overlaps = intersects

    ## Need a mathy name.
    #  \param *args variable number of Interval items.
    #  \retval bool If args cover all of self
    def is_covered_by(self, *args):
        """is_covered_by(i1, i2, ...) IFF arguments cover self"""
        return self in coverage(*args)

    __call__ = is_covered_by  # what was I thinking?

    ## filter(intersects(), ...) - like the built-in
    #  \param list_ A list of Interval objects
    #  \retval list_ of intersecting intervals.
    def filter(self, list_):
        """filter-like Format: filter in self.intersects"""
        return filter(self.intersects, list_)


## Combine 2 Interval to make a 3rd.
class _Combine(object):
    """Non arithmetic combos:
    C = A&B       maximal C such that C in A and C in B
    C = A|B       minimum C such that A in C and B in C
    C = A^B       minimum C such that A in C and B in C IFF A&C >= 0.
                  (that is, A|B IF there are no gaps)
    C = A%B       maximal C such that A^C^B == A|B and C not in (A, B).
    """
    ## just a mix-in
    __metaclass__ = abc.ABCMeta

    ## Minimal Interval containing both
    #  \param other Interval
    #  \retval interval Interval
    def __or__(self, other):
        """[x, y] | [x', y'] --> [min(x, x'), max(y, y')]"""
        return type(self)(min(self[0], other[0]),
                          max(self[1], other[1]))

    ## Maximal Interval contained by both
    #  \param other Interval
    #  \retval interval Interval
    def __and__(self, other):
        """[x, y] & [x', y'] --> [max(x, x'), min(y, y')]"""
        return type(self)(max(self[0], other[0]),
                          min(self[1], other[1]))

    ## Join IFF there are no gaps (unlike |, which fills in the gap)
    #  \param other Interval
    #  \retval interval Interval or None.
    def __xor__(self, other):
        """a ^ b: a | b If there are no gaps, None otherwise. """
        return None if int(self % other) > 0 else self | other  # guard

    ## Find the gap
    #  \param other Interval
    #  \retval interval Interval
    def __mod__(self, other):
        """a % b --> {a[1], b[0]} Gap"""
        return type(self)(self[1], other[0])


## Unary modiFiers
class _ModIfy(object):
    """Interval ModiFers"""

    ## just a mix-in
    __metaclass__ = abc.ABCMeta

    ## Transposition
    #  \retval interval \f$ [y, x] \f$
    def __invert__(self):
        """~[x, y] --> [y, x]"""
        return type(self)(*reversed(self))

    ## Increment right endpoint
    #  \param delta element
    #  \retval interval \f$ [x, y + \delta] \f$
    def __rshift__(self, delta):
        """[x, y] >> a --> [x, y+a]"""
        return type(self)(self[0], self[1] + delta)

    ## Decrement left endpoint
    #  \param delta element
    #  \retval interval  \f$ [x+\delta, y] \f$
    def __lshift__(self, delta):
        """[x, y] << a --> [x-a, y]"""
        return type(self)(self[0] - delta, self[1])

    ## Normal order, in-place
    def sort(self):
        """INPLACE normal ordering of endpoints"""
        if self[0] > self[1]:  # filter
            self = ~self


## A finite interval over an orderable field
class Interval(_Ob, _Comparator, _Combine, _Dx, _ModIfy, _Arithmetic):
    """interval = Interval(start, stop)

    start, stop:  almost anything that can be ordered, especially iF
    start-stop    makes sense, and For sure iF
    stop+stop     makes sense too.

    See Interval.__bases__.__doc__ For more.
    """


## A class just from datetime intervals (so NO _Arithmetic)
class TimeInterval(_Ob, _Comparator, _Dx, _ModIfy, _Combine):
    """TimeInterval(datetime.datetime, datetime.datetime')

    See TimeInterval.__bases__.__doc__ For more.
    """
    ## Convert support() to seconds
    def _to_float(self):
        return self.support().total_seconds()
