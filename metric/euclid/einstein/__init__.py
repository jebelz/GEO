"""Thus sub-packageimplements lexigraphic processing for Einstein summation
convention.

This init module has some wrappers and a hook that allows tensors to use
it. The heavy lifting is done in albert.py

## Normal Gymnastics:

>>>T.ii    # is a trace.
>>>T.iijk  # is just that (a partial trace)
>>>T.ijkl  # is a transpose relative to T.ijkl.

## Albert (via getattr, of course)
>>>T.[ij]k = (T.ijk - T.jik)/2
>>>T.{ij}k = (T.ijk + T.jik)/2


## Fixing a coordinate:

>>>T.ixj = T.e(Ellipsis, 0, Ellipsis) etc...

Any of the above can be compounded, so for example:

nine.ijmxyzijk

transposes m and k (indices 2 and 8), and then contracts on i (0, 6) and
j (1, 7).. the contraction order is irrelevant, finally, the remaining rank 5

tmp.kxyzm

has indices 1, 2, 3 fixed on x, y, z respectively, the resulting 9 values
are stuffed into a rank-2 Tensor.

The catch: running indices must run from i, j, .., s, t; however, they
need not be consecutive.
"""

__BUGS__ = "DD.iljk is not right--or is it?"

## \namespace geo.metric.euclid.einstein
# <a href="http://en.wikipedia.org/wiki/Einstein_notation">Einstein Notation
# </a>.
import operator
import itertools

__all__ = ('AXES', 'GymClass')

## Names of fixed axes.
AXES = 'xyz'


## Sign convention of symmetrizing symbols
SIGN = {'[': operator.neg, '{': operator.pos}


## Open: Close pairs of symbols-> [] and {}.
PAIR = {key: chr(ord(key)+2) for key in SIGN}


## List of 'extra' symbols that do not refer to indices.
XTRA = list(itertools.chain(*PAIR.iteritems()))


## The names of running indices.
RUN = 'ijklmnopqrstuv'


## Failed Gymnastic Requests.
class GymnasticsError(AttributeError):
    """Raised in __getattr__ If your request is not right."""


## Parser makes albert.Albert or not: this metaclass exists to do work
# before instantiation-- b/c I don't do work in __init__.
class Parser(type):
    """Parser(tensor, attr)

    Either makes a Albert object for parsing, or throws
    and AttributeError, which pushes things to super.__getattr__
    (some one else's super---very inappropriate intimacy."""

    ## Choose class or fail
    # \param tensor
    # \param attr
    # \returns albert.Albert()
    # \throws GymnasticsError If attr is not an index request on T.
    def __new__(mcs, tensor, attr):
        # 1st scan for errors.
        map(mcs._check_index, attr)
        if len(attr) - sum(map(attr.count, XTRA)) != tensor.rank:  # raise
            msg = "Wrong Number of Indices for rank {}: {}"
            raise GymnasticsError(msg.format(tensor.rank, attr))
        from .albert import Albert
        return Albert(tensor, attr)

    ## This is non-local control flow: If the attr request does not
    # belong here, then raise AttributeError, which kicks the request
    # up to ArrayLike, which checks for numpy requests.
    # \param x A character
    # \returns None
    # \throws GymnasticsError on junk x.
    @staticmethod
    def _check_index(x):
        if x in AXES or x in RUN or x in PAIR or x in PAIR.values():  # raise
            pass
        else:
            raise GymnasticsError("Can't understand: {}".format(x))


## Mixin inserts Einstein summation at just the right level.
class GymClass(object):
    """Tensor mixing interprets:

    T.xjkl  Slice
    T.ikj   Transpose
    T.iij   Contraction

    getattr(T, '[ij]kl')    antisym on indices 0,1
    getattr(T, 'i{jk}l')    sym on indices 1, 2

    It *must* throw AttributeError on a numpy request-- it is
    very much a fancy GOTO statement.
    """

    ## Allow direct index getting -note this throws errors on
    # numpy stuff and kicks it up to super-- not transparent.
    def __getattr__(self, attr):
        """self.e(*[AXES.index(item) If item in AXES else Ellipsis
        for item in attr])"""
        return Parser(self, attr)()
    
