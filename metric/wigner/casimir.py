"""Utilities for dealing with (half)-integer tensor product
representations of SU(2) and their projections onto Irreps.

Note on fractions: python's fraction library is ideal for dealing with
half integer math, but the math.factorial library does not accept them.
schur_weyl.pascal.factorial does--the code is in transition from the former
to the latter.
"""

## \namespace geo.metric.wigner.casimir \f$\frac{1}{2}\f$ and whole integer
# utilities.

import itertools
import fractions
import functools
import math


## \f$ \frac{1}{2} \f$
HALF = fractions.Fraction(1, 2)

## Unity, as a fractions.Fraction.
WHOLE = fractions.Fraction(1, 1)



## Natural number Test
# \param j degree
# \returns bool \f$j \in \mathbb{N}\f$
def is_in_N(j):
    """bool is natural number"""
    return math.floor(j) == float(j) and j >=0


## 1/2 Integer Test
# \param j degree
# \returns bool \f$j\f$ is 1/2-integer.
def ishalfint(j):
    """bool is half-int, not whole int."""
    return isrep(j) and not is_in_N(j)


## Valid representation degree
# \param j degree
# \returns bool \f$j\f$ is valid rep.
def isrep(j):
    """Is valid rep (1/2) or whole integer."""
    return (j >= 0) and (is_in_N(j) or is_in_N(2*j))


## Integer nint function
# \param j degree
# \returns int nearest.
def nint(j):
    """Nearest int."""
    return int(round(j))


## Dimension of rep
# \param j degree
# \return int \f$2j+1\f$
def dim(j):
    """Dimension of rep j = 2j+1"""
    return nint(2*j + 1)


## Casimir Invariant Eigenvalue
# \param j degree
# \return float \f$ j(j+1) \f$
def invariant(j):
    """Casimir invariant(j) = j(j+1) (as a float)."""
    return float(j*(j+1))


## Wrapped order iterator (see msort)
# \param j (degree)
# \yields \f$ [0, 1, \ldots, j, -j, -j+1, \ldots, -1] \f$ (or with 1/2 int)
def miter(j):
    """Iterates If valid m for a j."""
    # start at -j and pull out 2j+1 items, sort them.
    return iter(msort(itertools.islice(itertools.count(-j), dim(j))))


## Compare wrt to [0,1..,m,-m, ...-1] ordering
def mcmp(a, b):
    """mcmp(m1, m2) is a cmp function wrt to the python index ordering:
    [0, 1, ..., m, -m, -m+1, ..., -1]"""
    # zero is always 1st
    if b == 0:
        result = int(a != 0)
    elif (a >= 0 and b >= 0) or (a < 0 and b < 0):
        result = cmp(a, b)
    elif a >= 0:
        result = -1
    else:
        result = int(b != 0)
    return result


## Sort wrt to mcmp().
msort = functools.partial(sorted, cmp=mcmp)


## Degree iterator, for degrees less than or equal to an int or 1/2 int.
# \param j (degree)
# \yields \f$ [(j - \lfloor j \rfloor), \ldots, (j-1), j] \f$
def jiter(j):
    """iterate of j's from smallest sub-rep to j (in integer steps.)
    e.g:

    2 --> 0.0, 1.0, 2.0

    2.5 --> 0.5, 1.5, 2.5
    """
    return itertools.takewhile(functools.partial(cmp, j+1),
                               itertools.count(j-math.floor(j)))
